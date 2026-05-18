#' Fit a sparse GLM with PDB-selected lambda.
#'
#' Top-level entry point. Dispatches on the `family` argument,
#' runs PDB to choose \eqn{\lambda}, and fits a
#' warm-started 4-step regularisation path with FISTA.
#'
#'
#' @param X Numeric design matrix of shape `(n, p)`.
#' @param y Response. For `family = "cox"`, a 2-column matrix `(time, event)`.
#' @param family A character name in `c("gaussian","binomial","poisson",
#'   "exponential","gumbel","cox")` (default `"gaussian"`).
#' @param penalty Penalty name: one of `"lasso"`, `"scad"`,
#'   `"mcp"` (lowercase). Default `"lasso"`.
#' @param standardize Whether to standardise columns of `X` to zero mean and
#'   unit variance. Strongly recommended for PDB to be calibrated.
#'   Default is TRUE.
#' @param intercept Whether to fit an unpenalised intercept. Forced to
#'   `FALSE` for the Cox family.
#' @param lambda Optional fixed regularisation parameter; when `NULL`, PDB
#'   chooses it automatically.
#' @param lambda_method One of `"mc_exact"`, `"mc_gaussian"`, `"analytical"`.
#' @param relax If `TRUE`, run an unpenalised refit on the selected support
#'   after the regularisation path (debiasing).
#' @param mcp_gamma MCP concavity parameter when `penalty = "mcp"`. Default 3.0.
#' @param scad_a SCAD concavity parameter when `penalty = "scad"`. Default 3.7.
#' @param lambda_alpha Nominal level for PDB.
#' @param lambda_n_simu Number of Monte Carlo draws for PDB.
#' @param tol FISTA convergence tolerance at the final path step.
#' @param maxit Maximum number of iterations for FISTA if convergence is not yet reached.
#'
#' @return An object of class `c("pic.<family>", "pic")` with components:
#'   \describe{
#'     \item{beta}{Fitted coefficients (length `p`).}
#'     \item{intercept}{Fitted intercept or `NULL`.}
#'     \item{df}{The number of selected variables.}
#'     \item{selected}{Indices of the selected variables.}
#'     \item{lambda}{Selected (or fixed) \eqn{\lambda}.}
#'     \item{family}{The family object used.}
#'     \item{penalty}{The penalty object used.}
#'     \item{lambda_pdb}{Result from [lambda_pdb()] when `lambda = NULL`.}
#'     \item{n_samples}{The number of observations in the training set.}
#'     \item{n_features}{The number of variables in the training set.}
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 100
#' p <- 100
#' s <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- numeric(p)
#' beta[sample.int(p, s)] <- 3
#' y <- as.numeric(X %*% beta + rnorm(n))
#' fit <- pic(X, y, family = "gaussian", penalty = "scad")
#' coef(fit)
#' }
#' @export
pic <- function(
  X, y,
  family = c("gaussian", "binomial", "poisson", "exponential", "gumbel", "cox"),
  penalty = "lasso",
  standardize = TRUE,
  intercept = TRUE,
  lambda_method = c("mc_exact", "mc_gaussian", "analytical"),
  relax = FALSE,
  mcp_gamma = 3.0,
  scad_a = 3.7,
  lambda = NULL,
  lambda_alpha = 0.05,
  lambda_n_simu = 5000,
  tol = 1e-5,
  maxit = 10000
) {
  lambda_method <- match.arg(lambda_method)
  family_obj  <- get_family(match.arg(family))
  penalty_obj <- get_penalty(penalty, scad_a = scad_a, mcp_gamma = mcp_gamma)
  family_name <- family_obj$name

  # Response kind drives validation in check_y / sorting for survival.
  y_kind <- switch(family_name,
    gaussian    =, gumbel = "continuous",
    exponential =, poisson = "positive",
    binomial    = "binary",
    cox         = "survival"
  )
  
  if (family_name == "cox") intercept <- FALSE

  prep <- check_Xy(X, y, y_kind = y_kind, standardize_X = standardize)
  X_std <- prep$X
  y_use <- prep$y

  if (is.null(lambda)) {
    pdb_obj <- lambda_pdb(X_std, family_obj,
      n_simu = lambda_n_simu, alpha = lambda_alpha,
      method = lambda_method
    )
    lambda_val <- pdb_obj$value
  } else {
    pdb_obj <- NULL
    lambda_val <- as.numeric(lambda)
  }

  # 4-step warm-started path: lambda fractions ((i/4)^0.8), L1 warm-up,
  # final step uses the requested penalty with the tight tol.
  fracs <- (seq_len(4L) / 4)^0.8
  lam_path <- lambda_val * fracs
  tol_path <- c(rep(1e-3, 3L), tol)
  pen_path <- c(replicate(3L, .lasso(), simplify = FALSE), list(penalty_obj))

  coef_prev <- NULL
  int_prev <- NULL
  coef_ <- numeric(ncol(X_std))
  int_ <- 0.0   # placeholder; overwritten by the path loop below

  for (idx in seq_along(lam_path)) {
    out <- fista(
      X = X_std, y = y_use, family = family_obj,
      penalty = pen_path[[idx]], lambda_reg = lam_path[[idx]],
      fit_intercept = intercept, rel_tol = tol_path[[idx]],
      coef_init = coef_prev, intercept_init = int_prev, max_iter = maxit
    )
    coef_ <- out$coef
    int_ <- out$intercept
    coef_prev <- coef_
    int_prev <- int_
  }

  if (relax) {
    support <- which(coef_ != 0)
    if (length(support) > 0L) {
      out <- fista(
        X = X_std[, support, drop = FALSE], y = y_use, family = family_obj,
        penalty = penalty_obj, lambda_reg = 0.0,
        fit_intercept = intercept, rel_tol = 1e-3,
        coef_init = coef_[support], intercept_init = int_, max_iter = maxit
      )
      coef_full <- numeric(length(coef_))
      coef_full[support] <- out$coef
      coef_ <- coef_full
      int_ <- out$intercept
    }
  }

  sel_idx <- which(coef_ != 0)
  fit <- list(
    beta          = coef_,
    intercept     = if (intercept) int_ else NULL,
    df            = length(sel_idx),
    selected      = sel_idx,
    lambda        = lambda_val,
    family        = family_obj,
    penalty       = penalty_obj,
    lambda_pdb    = pdb_obj,
    n_samples     = nrow(X_std),
    n_features    = ncol(X_std)
  )

  class(fit) <- c(paste0("pic.", family_name), "pic")

  attr(fit, "preproc") <- list(X_mean      = prep$X_mean,
                               X_std       = prep$X_std,
                               standardize = standardize)

  # Cox-specific post-fit: baseline cumulative hazard / survival.
  if (family_name == "cox") {
    eta <- as.numeric(X_std %*% fit$beta)
    bf  <- baseline_functions(y_use[, 1L], y_use[, 2L], eta)
    fit$baseline_cumulative_hazard <- bf[, c("time", "cumulative_hazard")]
    fit$baseline_survival          <- bf[, c("time", "survival")]
    fit$unique_times               <- bf$time
    fit$censoring_rate             <- mean(y_use[, 2L] == 0) * 100
  }

  fit
}
