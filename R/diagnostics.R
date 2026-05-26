#' Phase-transition analysis of support recovery.
#'
#' For a fixed `(n, p)` configuration, varies the sparsity level `s`
#' from `0` to `s_max` and, for each `s`, estimates by Monte Carlo
#' (`m` replications) the probability that `pic()` recovers exactly
#' the support of the true coefficient vector. If several penalties are
#' supplied, one result is produced for each `(n, p, penalty)` combination.
#'
#' At each Monte Carlo replicate a design matrix is drawn from a
#' standard Gaussian, `s` features are sampled uniformly at random,
#' the response is generated under the chosen `type`, and one `pic`
#' fit is run for each requested penalty. The selected support is
#' compared to the truth and three metrics are stored:
#'
#' \describe{
#'   \item{`exact_recovery`}{1 if the selected set equals the true set, 0 otherwise.}
#'   \item{`tpr`}{\eqn{\left|\hat{S}\cap S\right| / |S|} - true positive rate. For `s = 0`, this is set to 1 when no true feature exists.}
#'   \item{`fdr`}{\eqn{\left|\hat{S}\setminus S\right| / \max{(|\hat{S}|, 1)}} - false discovery rate.}
#' }
#' 
#' ## Details on data sampling
#'
#' At each replicate, the design \eqn{X} is drawn iid
#' \eqn{\mathcal{N}(0, 1)}, the true support \eqn{S} is sampled
#' uniformly at random, and \eqn{\beta_j = } `beta_value` for
#' \eqn{j \in S}, \eqn{\beta_j = 0} otherwise. The linear predictor
#' is \eqn{\eta = X\beta}. The response \eqn{y} is then drawn
#' conditionally on \eqn{\eta} according to the requested family:
#'
#' \describe{
#'   \item{`"gaussian"`}{
#'     \eqn{y_i = \eta_i + \varepsilon_i} with
#'     \eqn{\varepsilon_i \sim \mathcal{N}(0, 1)}.
#'   }
#'   \item{`"binomial"`}{
#'     \eqn{y_i \sim \mathrm{Bernoulli}\!\left(\sigma(\eta_i)\right)},
#'     where \eqn{\sigma(z) = 1 / (1 + e^{-z})} is the logistic
#'     function.
#'   }
#'   \item{`"poisson"`}{
#'     \eqn{y_i \sim \mathrm{Poisson}\!\left(e^{\eta_i}\right)}.
#'   }
#'   \item{`"exponential"`}{
#'     \eqn{y_i \sim \mathrm{Exp}\!\left(\mathrm{rate} = e^{\eta_i}\right)}.
#'   }
#'   \item{`"gumbel"`}{
#'     \eqn{y_i = \eta_i + \varepsilon_i} with
#'     \eqn{\varepsilon_i \sim \mathrm{Gumbel}(0, 1)}, drawn as
#'     \eqn{-\log(-\log U_i)} for \eqn{U_i \sim \mathcal{U}(0, 1)}.
#'   }
#'   \item{`"cox"`}{
#'     Event times \eqn{T_i \sim \mathrm{Exp}\!\left(e^{\eta_i}\right)}
#'     and independent censoring times
#'     \eqn{C_i \sim \mathrm{Exp}(1)}. The response is the 2-column
#'     matrix \eqn{\bigl(\min(T_i, C_i),\, \mathbf{1}\{T_i \le C_i\}\bigr)}.
#'   }
#' }
#'
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
#' @param parallel Logical; if `TRUE`, distribute the `m` Monte Carlo
#'   replications of each `(n, p, penalty, s)` cell across multiple R
#'   processes via the `future` framework (`future::multisession` plan).
#' @param workers Integer; number of background R processes to use when
#'   `parallel = TRUE`. Ignored otherwise.
#'
#' @return An object of class `c("pic.phase_transition", "pic.diagnostic")`.
#' \item{s_grid}{`0:s_max`.}
#' \item{exact_recovery, tpr, fdr}{Matrices of shape `(length(n) * length(penalty), s_max + 1)` -
#'       one row per `(n, p, penalty)` curve, one column per sparsity level}
#' \item{curve_n, curve_p, curve_penalty}{Curve descriptors aligned with the rows of the metric matrices.}
#' \item{config}{A list of all configuration arguments for downstream plotting / reporting.}
#' \item{call}{The call.}
#'
#' @seealso [plot.pic.phase_transition()] for visualisation.
#'
#' @examples
#' \donttest{
#' pt <- phase_transition(n = 50, p = 100, type = "gaussian",
#'                        s_max = 8, m = 20,
#'                        penalty = c("lasso", "scad"), 
#'                        parallel = TRUE)
#' plot(pt)
#' }
#' @export
phase_transition <- function(
    n, p,
    type          = c("gaussian", "binomial", "poisson",
                      "exponential", "gumbel", "cox"),
    s_max,
    m             = 50,
    penalty       = "lasso",
    beta_value    = 3,
    lambda_method = "mc_exact",
    lambda_alpha  = 0.05,
    lambda_n_simu = 5000L,
    verbose       = FALSE,
    parallel      = FALSE,
    workers       = parallel::detectCores() - 1L
) {
  type <- match.arg(type)
  penalty <- match.arg(
    penalty,
    choices = c("lasso", "scad", "mcp"),
    several.ok = TRUE
  )
  
  n <- as.integer(n)
  p <- as.integer(p)
  
  if (length(n) != length(p))
    stop("`n` and `p` must have the same length.")
  
  s_max <- as.integer(s_max)
  m <- as.integer(m)
  
  s_grid <- 0:s_max
  K <- length(n)
  
  if (parallel) {
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)
  }
  
  results <- vector("list", K * length(penalty) * length(s_grid))
  row_id <- 0L
  
  for (k in seq_len(K)) {
    nk <- n[k]
    pk <- p[k]
    
    for (pen in penalty) {
      for (s in s_grid) {
        
        rep_results <- future.apply::future_lapply(
          seq_len(m),
          future.packages = "picreg",
          function(rep) {
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
              return(c(exact = 0, tpr = 0, fdr = 1))
            }
            
            sel <- sort(fit$selected)
            true <- sort(sim$true_features)
            
            c(
              exact = as.integer(identical(sel, true)),
              tpr = if (length(true) == 0L) 1 else length(intersect(sel, true)) / length(true),
              fdr = length(setdiff(sel, true)) / max(length(sel), 1L)

            )
          },
          future.seed = TRUE
        )
        
        rep_results <- do.call(rbind, rep_results)
        
        row_id <- row_id + 1L
        
        results[[row_id]] <- data.frame(
          n = nk,
          p = pk,
          penalty = pen,
          s = s,
          exact_recovery = mean(rep_results[, "exact"]),
          tpr = mean(rep_results[, "tpr"]),
          fdr = mean(rep_results[, "fdr"]),
          stringsAsFactors = FALSE
        )
        
        if (verbose) {
          message(sprintf(
            "(n=%d, p=%d, penalty=%s, s=%d) recovery=%.2f  tpr=%.2f  fdr=%.2f",
            nk, pk, pen, s,
            mean(rep_results[, "exact"]),
            mean(rep_results[, "tpr"]),
            mean(rep_results[, "fdr"])
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
#' @param ... Unused; present for S3 method consistency.
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
#' substitute for `mc_exact` in the asymptotic regime.
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
#' @return An object of class `c("pic.pdb_asymptotic", "pic.diagnostic")`.
#'  \item{n_grid, p, type, alpha, n_simu}{Configuration.}
#'  \item{stats_exact, stats_gaussian}{Lists of length `length(n_grid)` where
#'       each element is a numeric vector of length `n_simu` containing
#'       the simulated null statistics from the corresponding selector.}
#'  \item{lambda_exact, lambda_gaussian, lambda_analytical}{Numeric
#'       vectors of length `length(n_grid)` - the (1 - alpha) quantile
#'       under each selector at each `n`.}
#'  \item{call}{The call.}
#'
#' @seealso [plot.pic.pdb_asymptotic()] for visualisation.
#'
#' @examples
#' as_ <- pdb_asymptotic(n_grid = c(50, 200, 1000),
#'                       p = 200, type = "poisson")
#' plot(as_)
#' @export
pdb_asymptotic <- function(
    n_grid,
    p,
    type    = c("gaussian", "binomial", "poisson",
                "exponential", "gumbel", "cox"),
    alpha   = 0.05,
    n_simu  = 5000L,
    verbose = FALSE
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
#' @param ... Unused; present for S3 method consistency.
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