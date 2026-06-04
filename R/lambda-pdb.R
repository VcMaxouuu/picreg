#' Pivotal Detection Boundary regularization selector
#'
#' Computes the data-driven regularization parameter \eqn{\hat\lambda^{\rm PDB}_\alpha} 
#' using the Pivotal Detection Boundary (PDB) principle. The selected value is 
#' defined as the empirical \eqn{(1 - \alpha)} quantile of a null-distribution gradient statistic
#' \deqn{\hat\lambda^{\rm PDB}_\alpha =q_{1-\alpha}\left(\left\|\nabla \ell_0\right\|_\infty
#' \right),} where \eqn{\ell_0} denotes the loss evaluated under the null model. 
#' 
#' Under the null \eqn{\beta = 0}, the gradient of the loss carries only sampling noise.
#' The smallest \eqn{\lambda} large enough to dominate this noise is the natural threshold
#' separating signal from noise, i.e. the value above which a coefficient should be
#' kept rather than shrunk to zero. Calibrating \eqn{\hat\lambda} this way has two
#' consequences. First, the quantile depends only on the design matrix \eqn{X},
#' the family, and the level \eqn{\alpha} (not on the response \eqn{y}), so
#' cross-validation is no longer needed. Second, and more importantly, it leads
#' to sharper support recovery than prediction-error-based selectors such as
#' cross-validated lasso: by targeting the noise level of the gradient directly,
#' PDB controls the inclusion of noise variables and more reliably identifies
#' the true non-zero coefficients. Computing the quantile requires only the
#' distribution of \eqn{\|\nabla \ell_0\|_\infty}, which `lambda_pdb()` estimates
#' through one of the three methods below.
#'
#' ## Details on `method` option
#' 
#' The empirical quantile is obtained using one of the three following methods:
#' \describe{
#'   \item{`"mc_exact"`}{
#'     Family-aware Monte Carlo. For each of the `n_simu` draws, a
#'     response vector is sampled under the null model
#'     (\eqn{\beta = 0}) of the chosen family; the family-specific
#'     gradient of the loss at \eqn{\beta = 0} is evaluated on the
#'     fixed design \eqn{X}, and its supremum norm is recorded. The
#'     empirical \eqn{(1 - \alpha)} quantile of the `n_simu` recorded
#'     norms is returned. Most accurate but slowest of the three.
#'   }
#'   \item{`"mc_gaussian"`}{
#'     Monte Carlo under the Gaussian approximation of the null
#'     gradient. A central-limit argument gives
#'     \eqn{\nabla \ell_0 \approx \mathcal{N}(0,\, c(n)\, \Sigma_X / n)}
#'     with \eqn{\Sigma_X = X^\top X / n}. Each of the `n_simu` draws
#'     samples directly from this Gaussian and records its supremum
#'     norm — no family-specific evaluation needed. Family-agnostic
#'     and noticeably faster than `"mc_exact"`; valid in the regime
#'     where the CLT kicks in (moderate to large \eqn{n}).
#'   }
#'   \item{`"analytical"`}{
#'     Closed-form Bonferroni bound on the Gaussian tail. Combining a
#'     union bound over the \eqn{p} coordinates of the gradient with
#'     the standard Gaussian tail bound gives
#'     \deqn{\hat\lambda_\alpha^{\rm analytical} =
#'           \Phi^{-1}\!\left(1 - \alpha / (2p)\right)\, \sqrt{c(n) / n}.}
#'     Deterministic and \eqn{O(1)} — no simulation. Conservative
#'     (slightly over-estimates the true quantile) and gets looser as
#'     \eqn{p} grows, but useful when speed matters or when Monte
#'     Carlo is overkill.
#'   }
#' }
#' The `"mc_gaussian"` and `"analytical"` methods use a variance scaling
#' factor \eqn{c(n)} depending on the family:
#'
#' \itemize{
#'   \item Gaussian / Binomial / Poisson / Exponential: \eqn{c(n) = 1}
#'   \item Gumbel: \eqn{c(n) = \exp(2(\gamma + 1))}, with \eqn{\gamma} the Euler-Mascheroni constant.
#'   \item Cox: \eqn{c(n) = 1 / (4 \log n)}
#' }
#'
#' ## Computational cost
#'
#' Both Monte Carlo methods are dominated by a single \eqn{p \times n_{\rm simu}}
#' matrix product \eqn{X^\top R}, where \eqn{R} stacks the simulated
#' residuals. This product is dispatched to BLAS as a single \eqn{\tt gemm},
#' which is essentially as fast as a Monte Carlo selector can be on dense
#' designs. For very large \eqn{n} or \eqn{p}, however, the constant
#' becomes large and `"mc_exact"` may be unnecessarily expensive:
#'
#' \itemize{
#'   \item As \eqn{n} grows, the central-limit approximation tightens
#'     and `"mc_gaussian"` gives essentially the same \eqn{\hat\lambda}
#'     as `"mc_exact"` at a fraction of the cost (no family-specific
#'     residual generation, fewer dependencies on \eqn{y}-draws).
#'   \item `"analytical"` is \eqn{O(1)} and useful as a quick upper
#'     bound — slightly conservative, but accurate enough for triage
#'     and for very high \eqn{p} where the Bonferroni tail is tight.
#' }
#'
#' Practical rule of thumb: prefer `"mc_exact"` for small to moderate
#' problems (default), `"mc_gaussian"` once \eqn{n} grows and the
#' design is dense, and `"analytical"` when even the Monte Carlo cost is
#' a concern.
#'
#' ## Details on `alpha` option
#'
#' The level \eqn{\alpha} is the nominal Type-I error of the test of
#' \eqn{H_0\!\colon \beta = 0}. By construction of the quantile,
#' \deqn{\Pr\!\left(\left\|\nabla \ell_0\right\|_\infty
#'   > \hat\lambda^{\rm PDB}_\alpha \,\big|\, H_0\right) = \alpha,}
#' so under the null model no variable enters the active set with
#' probability at most \eqn{\alpha}. With the default
#' \eqn{\alpha = 0.05}, this caps the false-discovery rate at \eqn{5\%}
#' under the null: when the data carry no signal, picreg returns the
#' empty support \eqn{95\%} of the time.
#'
#'
#' @param X Numeric design matrix. Columns should typically be standardized
#'   to zero mean and unit variance.
#' @param family A PIC family specification. Can be either a family object
#'   or a character string accepted by [get_family()].
#' @param n_simu Number of Monte Carlo simulations used by the stochastic
#'   methods. Ignored when `method = "analytical"`.
#' @param alpha Nominal tail probability used to define the quantile level.
#' @param method One of "mc_exact", "mc_gaussian", or "analytical".
#' @return An object of class `"pic.lambda_pdb"`.
#' \item{value}{Selected regularization parameter \eqn{\hat\lambda^{\rm PDB}_\alpha}.}
#' \item{statistics}{Simulated null statistics used to estimate the quantile. 
#' `NULL` for the analytical method.}
#' \item{method}{Estimation method used.}
#' \item{alpha}{Quantile level parameter.}
#' \item{n_simu}{Number of Monte Carlo simulations.}
#'
#' @examples
#' X <- scale(matrix(rnorm(100 * 20), 100, 20))
#' lam <- lambda_pdb(
#'   X,
#'   family = "gaussian",
#'   method = "mc_exact"
#' )
#' print(lam)
#'
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
    warning("lambda_pdb: design matrix X does not appear standardized; ",
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

# Closed-form KKT threshold lambda_max = ||grad_beta L(0, beta0*)||_infty.
# Smallest lambda that kills every coordinate for lasso (and SCAD/MCP, whose
# derivative at 0 also equals lambda). Used as the entry point of the
# warm-start regularization path in `pic()`.
.lambda_max <- function(X, y, family, fit_intercept) {
  if (family$name == "cox") {
    lambda_max_cox_cpp(X, y[, 1L], y[, 2L])
  } else {
    y_vec <- if (is.matrix(y)) y[, 1L] else as.numeric(y)
    lambda_max_cpp(X, y_vec, family$name, fit_intercept)
  }
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
