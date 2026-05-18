#' Pivotal Detection Boundary regularisation parameter selector.
#'
#' Computes \eqn{\hat\lambda} as the \eqn{(1-\alpha)} empirical quantile of a
#' null-distribution gradient-norm statistic. Three methods:
#'
#' * `"mc_exact"` (default) — Monte Carlo simulation under the family's null
#'   distribution. Implemented in C++ (`src/lambda_pdb.cpp`).
#' * `"mc_gaussian"` — Gaussian CLT approximation. Family-agnostic, C++.
#' * `"analytical"` — Bonferroni bound. Closed form.
#'
#' Cox uses a dedicated permutation-based estimator
#' (`cox_null_grad_norms_cpp`) for `"mc_exact"`; `"mc_gaussian"` and
#' `"analytical"` fall back to the family-agnostic paths.
#'
#' @param X Standardised design matrix.
#' @param family A pic family object (see [pic_families]).
#' @param n_simu Number of Monte Carlo draws (default 5000).
#' @param alpha Nominal level (default 0.05).
#' @param method One of `"mc_exact"`, `"mc_gaussian"`, `"analytical"`.
#' @return A list with components `value` (the scalar \eqn{\hat\lambda}),
#'   `statistics` (the simulated null statistics, or `NULL` for `"analytical"`),
#'   `method`, `alpha`, `n_simu`.
#' @export
lambda_pdb <- function(X, family,
                       n_simu = 5000L, alpha = 0.05,
                       method = c("mc_exact", "mc_gaussian", "analytical")) {
  method <- match.arg(method)
  family <- get_family(family)
  if (!is.matrix(X)) stop("`X` must be a numeric matrix.")

  mu <- colMeans(X)
  sd <- sqrt(colMeans(sweep(X, 2L, mu, "-") ^ 2))
  if (max(abs(mu)) >= 1e-4 || max(abs(sd - 1.0)) >= 1e-4) {
    warning("lambda_pdb: design matrix X does not appear standardised; ",
            "results may be unreliable.")
  }

  # Variance scaling c(n) used by the family-agnostic Gaussian / analytical
  # paths. Hardcoded here as a single source of truth: 6 cases, no abstraction.
  n <- nrow(X)
  c_n <- switch(family$name,
                gumbel = exp(2 * (-digamma(1) + 1)),
                cox    = 1.0 / (4.0 * log(n)),
                1.0)   # gaussian, binomial, poisson, exponential

  # Cox uses a custom permutation-based pivotal estimator for mc_exact.
  # mc_gaussian and analytical fall back to the family-agnostic paths.
  if (family$name == "cox" && method == "mc_exact") {
    out <- cox_null_grad_norms_cpp(X, as.integer(n_simu), alpha)
  } else {
    out <- switch(
      method,
      mc_exact    = pdb_mc_exact_cpp(X, family$name, as.integer(n_simu), alpha),
      mc_gaussian = pdb_mc_gaussian_cpp(X, c_n, as.integer(n_simu), alpha),
      analytical  = .pdb_analytical(X, c_n, alpha)
    )
  }
  out$method <- method
  out$alpha  <- alpha
  out$n_simu <- if (method == "analytical") NA_integer_ else as.integer(n_simu)
  class(out) <- "pic.lambda_pdb"
  out
}

# Bonferroni closed-form quantile. O(1); not worth porting to C++.
.pdb_analytical <- function(X, c_n, alpha) {
  p <- ncol(X)
  n <- nrow(X)
  list(value      = stats::qnorm(1.0 - alpha / (2.0 * p)) * sqrt(c_n / n),
       statistics = NULL)
}

#' @export
print.pic.lambda_pdb <- function(x, ...) {
  cat("pic lambda_pdb (method = ", x$method,
      ", alpha = ", format(x$alpha, digits = 3), ")\n", sep = "")
  cat("  lambda = ", format(x$value, digits = 6), "\n", sep = "")
  invisible(x)
}
