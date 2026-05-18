#' Phase-transition analysis of support recovery.
#'
#' For a fixed `(n, p)` configuration, varies the sparsity level `s`
#' from `0` to `s_max` and, for each `s`, estimates by Monte Carlo
#' (`m` replications) the probability that `pic()` recovers exactly
#' the support of the true coefficient vector.
#'
#' At each Monte Carlo replicate a design matrix is drawn from a
#' standard Gaussian, `s` features are sampled uniformly at random,
#' the response is generated under the chosen `type`, and one `pic`
#' fit is run for each requested penalty. The selected support is
#' compared to the truth and three metrics are stored:
#'
#' \describe{
#'   \item{exact_recovery}{1 if the selected set equals the true set, 0 otherwise.}
#'   \item{tpr}{|selected \eqn{\cap} true| / s — true positive rate. For `s = 0`, this is set to 1 when no true feature exists.}
#'   \item{fdr}{|selected \\ true| / max(|selected|, 1) — false discovery rate.}
#' }
#'
#' `n` and `p` may be scalars (one configuration) or vectors of equal length
#' (one configuration per `(n[k], p[k])` pair). If several penalties are
#' supplied, one curve is produced for each `(n, p, penalty)` combination.
#'
#' @param n Integer or integer vector — number of observations per configuration.
#' @param p Integer or integer vector of the same length as `n` — number of features.
#' @param type Family name: `"gaussian"`, `"binomial"`,
#'   `"poisson"`, `"exponential"`, `"gumbel"`, or `"cox"`.
#' @param s_max Largest sparsity level evaluated. Must satisfy `s_max < min(p)`.
#' @param m Number of Monte Carlo replications per `s`.
#' @param penalty One or more penalties for `pic()`: `"lasso"`, `"scad"`, or `"mcp"`.
#' @param beta_value Magnitude of the non-zero coefficients used to generate
#'   the response. The sign is fixed to `+`.
#' @param lambda_method Passed to `pic()`. Default `"mc_exact"`.
#' @param lambda_alpha Nominal level for the PDB selector.
#' @param lambda_n_simu Monte Carlo size for the PDB selector.
#' @param verbose Logical; if `TRUE`, prints a one-line progress message per `(n, p, penalty, s)`.
#'
#' @return An object of class `c("pic.phase_transition", "pic.diagnostic")`
#'   with components:
#'   \describe{
#'     \item{s_grid}{`0:s_max`.}
#'     \item{exact_recovery, tpr, fdr}{matrices of shape `(length(n) * length(penalty), s_max + 1)` —
#'       one row per `(n, p, penalty)` curve, one column per sparsity level.}
#'     \item{curve_n, curve_p, curve_penalty}{Curve descriptors aligned with the rows of the metric matrices.}
#'     \item{config}{A list of all configuration arguments for downstream plotting / reporting.}
#'     \item{call}{The call.}
#'   }
#'
#' @seealso [plot.pic.phase_transition()] for visualisation.
#'
#' @examples
#' \dontrun{
#' # Single configuration, multiple penalties
#' pt <- phase_transition(n = 200, p = 500, type = "gaussian",
#'                        s_max = 20, m = 50,
#'                        penalty = c("lasso", "scad", "mcp"))
#' plot(pt)
#'
#' # Two configurations and two penalties on the same plot
#' pt <- phase_transition(n = c(200, 400), p = c(500, 500),
#'                        type = "binomial", s_max = 15, m = 50,
#'                        penalty = c("lasso", "scad"))
#' plot(pt)
#' }
#' @export
phase_transition <- function(
    n, p,
    type          = c("gaussian", "binomial", "poisson",
                      "exponential", "gumbel", "cox"),
    s_max,
    m             = 100,
    penalty       = c("lasso", "scad", "mcp"),
    beta_value    = 3,
    lambda_method = "mc_exact",
    lambda_alpha  = 0.05,
    lambda_n_simu = 5000L,
    verbose       = TRUE
) {
  type <- match.arg(type)
  penalty <- match.arg(penalty, several.ok = TRUE)
  
  n <- as.integer(n)
  p <- as.integer(p)
  
  if (length(n) != length(p))
    stop("`n` and `p` must have the same length.")
  if (any(n <= 0L) || any(p <= 0L))
    stop("`n` and `p` must be positive.")
  
  s_max <- as.integer(s_max)
  if (length(s_max) != 1L || s_max < 0L)
    stop("`s_max` must be a non-negative scalar integer.")
  if (s_max >= min(p))
    stop("`s_max` must be strictly less than min(p).")
  
  m <- as.integer(m)
  if (length(m) != 1L || m < 1L)
    stop("`m` must be a positive scalar integer.")
  
  s_grid <- 0:s_max
  K <- length(n)
  
  results <- vector("list", K * length(penalty) * length(s_grid))
  row_id <- 0L
  
  for (k in seq_len(K)) {
    nk <- n[k]
    pk <- p[k]
    
    for (pen in penalty) {
      for (s in s_grid) {
        exact <- numeric(m)
        tpr_s <- numeric(m)
        fdr_s <- numeric(m)
        
        for (rep in seq_len(m)) {
          sim <- .generate_recovery_data(nk, pk, s, type, beta_value)
          
          fit <- tryCatch(
            pic(
              sim$X,
              sim$y,
              family = type,
              penalty = pen,
              lambda_method = lambda_method,
              lambda_alpha = lambda_alpha,
              lambda_n_simu = lambda_n_simu
            ),
            error = function(e) NULL
          )
          
          if (is.null(fit)) {
            exact[rep] <- 0
            tpr_s[rep] <- 0
            fdr_s[rep] <- 1
            next
          }
          
          sel <- sort(fit$selected)
          true <- sort(sim$true_features)
          
          exact[rep] <- as.integer(identical(sel, true))
          
          tpr_s[rep] <- if (length(true) == 0L) {
            1
          } else {
            length(intersect(sel, true)) / length(true)
          }
          
          fdr_s[rep] <- length(setdiff(sel, true)) / max(length(sel), 1L)
        }
        
        row_id <- row_id + 1L
        
        results[[row_id]] <- data.frame(
          n = nk,
          p = pk,
          penalty = pen,
          s = s,
          exact_recovery = mean(exact),
          tpr = mean(tpr_s),
          fdr = mean(fdr_s),
          stringsAsFactors = FALSE
        )
        
        if (verbose) {
          message(sprintf(
            "(n=%d, p=%d, penalty=%s, s=%d) recovery=%.2f  tpr=%.2f  fdr=%.2f",
            nk, pk, pen, s,
            mean(exact), mean(tpr_s), mean(fdr_s)
          ))
        }
      }
    }
  }
  
  results <- do.call(rbind, results)
  
  config <- list(
    n = n,
    p = p,
    type = type,
    s_max = s_max,
    m = m,
    penalty = penalty,
    beta_value = beta_value,
    lambda_method = lambda_method,
    lambda_alpha = lambda_alpha,
    lambda_n_simu = lambda_n_simu
  )
  
  out <- list(
    results = results,
    exact_recovery = .phase_metric_table(results, "exact_recovery"),
    tpr = .phase_metric_table(results, "tpr"),
    fdr = .phase_metric_table(results, "fdr"),
    s_grid = s_grid,
    config = config,
    call = match.call()
  )
  
  class(out) <- c("pic.phase_transition", "pic.diagnostic")
  out
}


.phase_metric_table <- function(results, metric) {
  curves <- unique(results[c("n", "p", "penalty")])
  s_vals <- sort(unique(results$s))
  
  out <- curves
  
  for (s in s_vals) {
    col <- paste0("s_", s)
    out[[col]] <- NA_real_
    
    for (i in seq_len(nrow(curves))) {
      idx <- results$n == curves$n[i] &
        results$p == curves$p[i] &
        results$penalty == curves$penalty[i] &
        results$s == s
      
      out[[col]][i] <- results[[metric]][idx]
    }
  }
  
  out
}


#' Print phase-transition analysis.
#'
#' @param x A `pic.phase_transition` object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#' @export
print.pic.phase_transition <- function(x, ...) {
  cfg <- x$config
  K <- length(x$curve_n)
  
  cat("pic phase-transition analysis\n")
  cat("  type      : ", cfg$type, "\n", sep = "")
  cat("  penalties : ", paste(cfg$penalty, collapse = ", "), "\n", sep = "")
  cat("  s grid    : 0..", cfg$s_max, "\n", sep = "")
  cat("  reps      : ", cfg$m, "\n", sep = "")
  cat("  curves    : ", K, "\n", sep = "")
  
  for (k in seq_len(K)) {
    cat(sprintf(
      "    [%d] n = %d, p = %d, penalty = %s\n",
      k,
      x$curve_n[k],
      x$curve_p[k],
      x$curve_penalty[k]
    ))
  }
  
  invisible(x)
}


# ---- internal: data generation for recovery experiments -----------------
#
# For each call returns a list with `X`, `y`, and the indices of the
# non-zero true coefficients (`true_features`). All families use a
# standard Gaussian design and a sparse beta with `s` entries set to
# `beta_value` at uniformly-random positions.
#
# Survival case (Cox): event times follow Exp(rate = exp(X %*% beta)),
# censoring times follow Exp(1); event indicator is min(T, C).
.generate_recovery_data <- function(n, p, s, type, beta_value) {
  X <- matrix(stats::rnorm(n * p), n, p)
  
  true_features <- if (s == 0L) {
    integer(0L)
  } else {
    sort(sample.int(p, s))
  }
  
  beta <- numeric(p)
  beta[true_features] <- beta_value
  eta <- as.numeric(X %*% beta)
  
  y <- switch(
    type,
    gaussian    = eta + stats::rnorm(n),
    binomial    = stats::rbinom(n, size = 1, prob = 1 / (1 + exp(-eta))),
    poisson     = stats::rpois(n, lambda = exp(eta)),
    exponential = stats::rexp(n, rate = exp(eta)),
    gumbel      = eta - log(-log(stats::runif(n))),
    cox         = {
      t_event <- stats::rexp(n, rate = exp(eta))
      t_cens  <- stats::rexp(n, rate = 1)
      time    <- pmin(t_event, t_cens)
      event   <- as.integer(t_event <= t_cens)
      cbind(time = time, event = event)
    }
  )
  
  list(X = X, y = y, true_features = true_features)
}


#' Asymptotic behaviour of the PDB null distribution.
#'
#' For each `n` in `n_grid`, draws a standardised Gaussian design matrix
#' of shape `(n, p)` and computes the null gradient-norm statistic via
#' the three available selectors: `"mc_exact"`, `"mc_gaussian"`, and
#' `"analytical"`. Stores the simulated Monte Carlo statistics and the
#' three resulting \eqn{\hat\lambda} values per `n`.
#'
#' The intended use is to **visualise the convergence** of the exact
#' family-specific null distribution to the Gaussian approximation as
#' `n` grows — i.e., to check empirically that `mc_gaussian` is a valid
#' (and much faster) substitute for `mc_exact` in the asymptotic
#' regime.
#'
#' @param n_grid Integer vector of sample sizes to evaluate.
#' @param p Number of features (scalar integer).
#' @param type Family name: `"gaussian"`, `"binomial"`,
#'   `"poisson"`, `"exponential"`, `"gumbel"`, or `"cox"`.
#' @param alpha Nominal level used for the (1 - alpha) quantile.
#' @param n_simu Monte Carlo size for each selector.
#' @param verbose Logical; if `TRUE`, prints a one-line progress
#'   message per `n`.
#'
#' @return An object of class `c("pic.pdb_asymptotic", "pic.diagnostic")`
#'   with components:
#'   \describe{
#'     \item{n_grid, p, type, alpha, n_simu}{Configuration.}
#'     \item{stats_exact, stats_gaussian}{Lists of length `length(n_grid)`;
#'       each element is a numeric vector of length `n_simu` containing
#'       the simulated null statistics from the corresponding selector.}
#'     \item{lambda_exact, lambda_gaussian, lambda_analytical}{Numeric
#'       vectors of length `length(n_grid)` — the (1 - alpha) quantile
#'       under each selector at each `n`.}
#'     \item{call}{The call.}
#'   }
#'
#' @seealso [plot.pic.pdb_asymptotic()] for visualisation.
#'
#' @examples
#' \dontrun{
#' as_ <- pdb_asymptotic(n_grid = c(50, 200, 1000),
#'                       p = 200, type = "poisson")
#' plot(as_)
#' }
#' @export
pdb_asymptotic <- function(
    n_grid,
    p,
    type    = c("gaussian", "binomial", "poisson",
                "exponential", "gumbel", "cox"),
    alpha   = 0.05,
    n_simu  = 5000L,
    verbose = TRUE
) {
  type <- match.arg(type)

  n_grid <- as.integer(n_grid)
  if (length(n_grid) < 1L) stop("`n_grid` must be non-empty.")
  if (any(n_grid <= 0L))   stop("`n_grid` must contain positive integers.")

  p <- as.integer(p)
  if (length(p) != 1L || p <= 0L)
    stop("`p` must be a positive scalar integer.")

  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || alpha <= 0 || alpha >= 1)
    stop("`alpha` must lie in (0, 1).")

  n_simu <- as.integer(n_simu)
  if (length(n_simu) != 1L || n_simu < 1L)
    stop("`n_simu` must be a positive scalar integer.")

  fam <- get_family(type)
  K   <- length(n_grid)

  stats_exact       <- vector("list", K)
  stats_gaussian    <- vector("list", K)
  lambda_exact      <- numeric(K)
  lambda_gaussian   <- numeric(K)
  lambda_analytical <- numeric(K)

  for (k in seq_len(K)) {
    nk <- n_grid[k]
    if (verbose)
      message(sprintf("[n = %d, p = %d] running pdb_asymptotic...", nk, p))

    # Standardised Gaussian design (n divisor, matching check_X).
    X  <- matrix(stats::rnorm(nk * p), nk, p)
    mu <- colMeans(X)
    X  <- sweep(X, 2L, mu, "-", check.margin = FALSE)
    sv <- sqrt(colMeans(X^2))
    sv[sv == 0] <- 1
    X  <- sweep(X, 2L, sv, "/", check.margin = FALSE)

    pdb_e <- lambda_pdb(X, fam, n_simu = n_simu, alpha = alpha,
                        method = "mc_exact")
    pdb_g <- lambda_pdb(X, fam, n_simu = n_simu, alpha = alpha,
                        method = "mc_gaussian")
    pdb_a <- lambda_pdb(X, fam, n_simu = n_simu, alpha = alpha,
                        method = "analytical")

    stats_exact[[k]]     <- pdb_e$statistics
    stats_gaussian[[k]]  <- pdb_g$statistics
    lambda_exact[k]      <- pdb_e$value
    lambda_gaussian[k]   <- pdb_g$value
    lambda_analytical[k] <- pdb_a$value
  }

  structure(
    list(
      n_grid            = n_grid,
      p                 = p,
      type              = type,
      alpha             = alpha,
      n_simu            = n_simu,
      stats_exact       = stats_exact,
      stats_gaussian    = stats_gaussian,
      lambda_exact      = lambda_exact,
      lambda_gaussian   = lambda_gaussian,
      lambda_analytical = lambda_analytical,
      call              = match.call()
    ),
    class = c("pic.pdb_asymptotic", "pic.diagnostic")
  )
}


#' Print PDB asymptotic diagnostic.
#'
#' @param x A `pic.pdb_asymptotic` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.pic.pdb_asymptotic <- function(x, ...) {
  K <- length(x$n_grid)
  cat("pic PDB asymptotic diagnostic\n")
  cat("  type      : ", x$type, "\n", sep = "")
  cat("  p         : ", x$p, "\n", sep = "")
  cat("  alpha     : ", format(x$alpha, digits = 3), "\n", sep = "")
  cat("  n_simu    : ", x$n_simu, "\n", sep = "")
  cat("  n_grid    : ", paste(x$n_grid, collapse = ", "), "\n", sep = "")
  cat("\n  lambda values per n:\n")
  df <- data.frame(
    n          = x$n_grid,
    mc_exact   = signif(x$lambda_exact, 4),
    mc_gauss   = signif(x$lambda_gaussian, 4),
    analytical = signif(x$lambda_analytical, 4)
  )
  print(df, row.names = FALSE)
  invisible(x)
}