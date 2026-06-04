#' Sparse linear regression using the Pivotal Information Criterion.
#'
#' Fits sparse regression, classification, and survival models built on a
#' linear predictor, using a family-specific loss whose gradient at the null
#' is (asymptotically) pivotal, combined with a sparsity-inducing penalty.
#' Supported families are Gaussian, binomial, Poisson, exponential, Gumbel,
#' and Cox proportional hazards.
#'
#' Minimises the objective function
#' \deqn{\hat\beta = \arg\min_\beta\; L(\beta)
#'        \;+\; \lambda_\alpha^{\rm PDB}\,\mathrm{pen}(\beta)
#'        \quad\text{with}\quad
#'        L(\beta) = \phi\!\left(\ell_n\!\left(g(\beta)\right)\right),}
#' where the regularization parameter is selected automatically using
#' the Pivotal Detection Boundary (PDB) principle implemented in
#' [lambda_pdb()]. The PDB choice is not just a substitute for
#' cross-validation: by calibrating \eqn{\lambda} against a pivotal
#' null-distribution statistic rather than out-of-sample prediction
#' error, it yields sharper support recovery - fewer false positives
#' on noise variables and more accurate identification of the true
#' non-zero coefficients. Here, \eqn{\phi} and \eqn{g} are
#' family-specific transformations chosen so that the gradient at the
#' null has a distribution free of the nuisance parameter, which is
#' precisely what makes the use of \eqn{\lambda^{\rm PDB}_\alpha}
#' valid. Optimization is performed using a warm-started FISTA
#' regularization path.
#'
#' ## Details on family-specific losses
#'
#' \describe{
#'   \item{`"gaussian"`}{
#'     Uses the mean squared error for \eqn{\ell_n}, with
#'     \eqn{\phi(\cdot) = \sqrt{\cdot}} and \eqn{g(\cdot) = \mathrm{Id}},
#'     giving the square-root lasso loss
#'     \deqn{L(\beta) = \sqrt{\frac{1}{n}\sum_{i=1}^n
#'       \left(y_i - (\beta_0 + \mathbf{x}_i^\top\beta)\right)^2}.}
#'     This makes the gradient at \eqn{\beta = 0} scale-free, so the PDB
#'     threshold needs no estimate of the noise standard deviation.
#'   }
#'   \item{`"binomial"`}{
#'     A variance-stabilized transformation of the Bernoulli likelihood.
#'     With \eqn{\phi(\cdot) = \mathrm{Id}} and the logistic link
#'     \eqn{\theta_i = g(\eta_i) = (1 + e^{-\eta_i})^{-1}}, the loss is
#'     \deqn{L(\beta) = \frac{1}{n}\sum_{i=1}^n
#'       \left(2 y_i \sqrt{\frac{1 - \theta_i}{\theta_i}}
#'             + 2 (1 - y_i) \sqrt{\frac{\theta_i}{1 - \theta_i}}\right).}
#'     The classical logistic link is preserved; only the loss itself is
#'     modified to obtain a pivotal gradient at the null.
#'   }
#'   \item{`"poisson"`}{
#'     The same pivotalization applied to the Poisson likelihood. With
#'     \eqn{\phi(\cdot) = \mathrm{Id}} and the canonical log link
#'     \eqn{\theta_i = g(\eta_i) = e^{\eta_i}}, the loss is
#'     \deqn{L(\beta) = \frac{1}{n}\sum_{i=1}^n
#'       \left(\frac{2 y_i}{\sqrt{\theta_i}} + 2 \sqrt{\theta_i}\right).}
#'     The canonical log link is kept unchanged; pivotality is obtained
#'     through the loss rather than through the link.
#'   }
#'   \item{`"exponential"`}{
#'     Uses the standard Exponential negative log-likelihood directly (no
#'     transformation). With \eqn{\phi(\cdot) = \mathrm{Id}} and
#'     \eqn{\theta_i = g(\eta_i) = e^{\eta_i}}, the loss is
#'     \deqn{L(\beta) = \frac{1}{n}\sum_{i=1}^n
#'       \left(\log\theta_i + \frac{y_i}{\theta_i}\right).}
#'   }
#'   \item{`"gumbel"`}{
#'     A location-scale model with extreme-value noise. The base
#'     log-likelihood is
#'     \deqn{\ell_n(\theta, \sigma) = \log(\sigma)
#'       + \frac{1}{n}\sum_{i=1}^n \left(z_i + e^{-z_i}\right),
#'       \qquad z_i = \frac{y_i - \theta_i}{\sigma}.}
#'     With \eqn{\phi(\cdot) = \exp(\cdot)} and the identity link
#'     \eqn{\theta_i = g(\eta_i) = \eta_i}, the objective is
#'     \eqn{L(\beta, \sigma) = \exp(\ell_n(\theta, \sigma))}. The scale
#'     \eqn{\sigma} is re-estimated internally by maximum likelihood at
#'     every iteration; the user does not supply it.
#'   }
#'   \item{`"cox"`}{
#'     A square-root-transformed Cox partial log-likelihood. With
#'     \eqn{\phi(\cdot) = \sqrt{\cdot}} and the identity link
#'     \eqn{\theta_i = g(\eta_i) = \eta_i}, the objective is
#'     \deqn{L(\beta) = \sqrt{-\frac{1}{n}\left(
#'       \sum_{i=1}^n \delta_i \eta_i
#'       - \sum_{i=1}^n \delta_i
#'         \log\left(\sum_{j \in R_i} e^{\eta_j}\right)\right)},}
#'     where \eqn{\delta_i} is the event indicator and \eqn{R_i} the risk
#'     set at event time \eqn{t_i}. The Breslow approximation is used for
#'     tied event times, and the model is always fitted without an intercept.
#'   }
#' }
#'
#' @param X Numeric design matrix of shape `(n, p)`, where rows are
#'   observations and columns are candidate predictors. `X` may be a matrix
#'   or a data frame coercible to a numeric matrix. It must contain at least
#'   two columns and no `NA`, `NaN`, or infinite values. 
#' @param y Response variable with length `n`. For Gaussian and Gumbel models,
#'   `y` must be numeric. For Binomial models, `y` must contain only `0` and
#'   `1`. For Poisson and Exponential models, `y` must be non-negative. For
#'   `family = "cox"`, `y` must be a two-column numeric matrix `(time, event)`,
#'   where `time` is non-negative and `event` is coded as `0` or `1`.
#' @param family A character name determining the response distribution and
#'   the associated loss function used during optimization. See
#'   \strong{Details on family-specific losses} below for the exact
#'   objective functions associated with each family. Default `"gaussian"`.
#' @param penalty Penalty name: one of `"lasso"`, `"scad"`, or `"mcp"`
#'   (lowercase). Default `"lasso"`. See [pic_penalties] for the precise
#'   form of each penalty and the role of `scad_a` / `mcp_gamma`.
#' @param standardize Whether to standardize columns of `X` to zero mean and
#'   unit variance prior to computing the PDB regularization parameter and
#'   fitting the model. Strongly recommended for proper PDB calibration.
#'   When `TRUE`, the optimization is performed on the standardized design
#'   matrix and the stored fitted coefficients therefore correspond to the
#'   standardized scale. Use [coef.pic()] to recover coefficients on the
#'   original scale of `X`. Default is `TRUE`. Standardization centers `X`;
#'   for non-Cox families the resulting column-mean shift is absorbed by the
#'   intercept, so combining `standardize = TRUE` with `intercept = FALSE`
#'   assumes the data are already centered (see `intercept`).
#' @param intercept Whether to fit an unpenalized intercept. Forced to
#'   `FALSE` for the Cox family. We assume the data has been centered if intercept = FALSE.
#' @param lambda Optional user-supplied regularization parameter. When `NULL`,
#'   [lambda_pdb()] is used to select it automatically. Providing `lambda`
#'   avoids recomputing the PDB calibration, which is useful for example when
#'   fitting several penalties (`"lasso"`, `"scad"`, `"mcp"`) on the same
#'   dataset, since the PDB choice of \eqn{\lambda} does not depend on the
#'   penalty itself.
#' @param lambda_method One of `"mc_exact"`, `"mc_gaussian"`, `"analytical"`. See
#'   [lambda_pdb()] for details. The default `"mc_exact"` is a safe choice
#'   for small to moderate designs; for very large \eqn{n} or \eqn{p},
#'   `"mc_gaussian"` matches `"mc_exact"` closely at a fraction of the
#'   cost (it amortises only one BLAS \eqn{\tt gemm}, no family-specific
#'   residual draws), and `"analytical"` is closed-form.
#' @param relax If `TRUE`, run an unpenalized refit on the selected support
#'   after the regularization path (debiasing).
#' @param mcp_gamma MCP concavity parameter when `penalty = "mcp"`. Default 3.0.
#' @param scad_a SCAD concavity parameter when `penalty = "scad"`. Default 3.7.
#' @param lambda_alpha Nominal level for PDB. See [lambda_pdb()].
#' @param lambda_n_simu Number of Monte Carlo draws for PDB. See [lambda_pdb()].
#' @param tol FISTA convergence tolerance at the final path step.
#' @param path_length Number of points in the warm-start regularization
#'   path running from \eqn{\lambda_{\max}} down to \eqn{\lambda^{\rm PDB}}.
#'   Default 10.
#' @param maxit Maximum number of iterations for FISTA if convergence is not yet reached.
#'
#' @return An object of class `c("pic.<family>", "pic")`.
#'     \item{beta}{Fitted coefficients (length `p`).}
#'     \item{intercept}{Fitted intercept or `NULL`.}
#'     \item{df}{The number of nonzero coefficients.}
#'     \item{selected}{Indices of the nonzero coefficients.}
#'     \item{lambda}{\eqn{\lambda} used during training (PDB or user-supplied).}
#'     \item{family}{The family object used.}
#'     \item{penalty}{The penalty object used.}
#'     \item{lambda_pdb}{Result from [lambda_pdb()] if \eqn{\lambda} is not user-supplied.}
#'     \item{n_samples}{The number of observations in the training set.}
#'     \item{n_features}{The number of variables in the training set.}
#' For `family = "cox"`, the fit additionally carries:
#'     \item{baseline_cumulative_hazard}{Estimated Breslow baseline cumulative hazard function.}
#'     \item{baseline_survival}{Estimated baseline survival function.}
#'     \item{unique_times}{Sorted unique event times used in the baseline estimation.}
#'     \item{censoring_rate}{Percentage of censored observations in the training set.}
#'     
#' @examples
#' data(QuickStartExample)
#' X <- QuickStartExample$X
#' y <- QuickStartExample$y
#' fit <- pic(X, y, family = "gaussian", penalty = "scad")
#' fit$beta
#' fit$selected
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
  lambda_n_simu = 2000,
  tol = 1e-5,
  path_length = 10L,
  maxit = 1000
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

  lambda_max_val <- .lambda_max(X_std, y_use, family_obj, intercept)
  if (lambda_val == 0) {
    lam_path <- 0
    pen_path <- list(penalty_obj)
    tol_path <- tol
  } else if (lambda_val >= lambda_max_val) {
    lam_path <- lambda_val
    pen_path <- list(penalty_obj)
    tol_path <- tol
  } else {
    K <- max(1L, as.integer(path_length))
    if (K == 1L) {
      lam_path <- lambda_val
      pen_path <- list(penalty_obj)
      tol_path <- tol
    } else {
      lam_path <- exp(seq(log(lambda_max_val), log(lambda_val),
                          length.out = K))
      pen_path <- c(replicate(K - 1L, .lasso(), simplify = FALSE),
                    list(penalty_obj))
      tol_path <- c(rep(1e-3, K - 1L), tol)
    }
  }

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
  feature_names <- prep$feature_names
  selected_out <- if (!is.null(feature_names)) feature_names[sel_idx] else sel_idx

  fit <- list(
    beta          = coef_,
    intercept     = if (intercept) int_ else NULL,
    df            = length(sel_idx),
    selected      = selected_out,
    lambda        = lambda_val,
    family        = family_obj,
    penalty       = penalty_obj,
    lambda_pdb    = pdb_obj,
    n_samples     = nrow(X_std),
    n_features    = ncol(X_std)
  )

  class(fit) <- c(paste0("pic.", family_name), "pic")

  # Per-column default value grids are only retained for Cox fits, where
  # `feature_effects_on_survival()` uses them to spare the user from
  # passing the training X back in. Kept NULL for the other families to
  # avoid useless memory bloat.
  attr(fit, "preproc") <- list(
    X_mean         = prep$X_mean,
    X_std          = prep$X_std,
    standardize    = standardize,
    feature_names  = feature_names,
    feature_values = if (family_name == "cox") prep$feature_values else NULL
  )

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
