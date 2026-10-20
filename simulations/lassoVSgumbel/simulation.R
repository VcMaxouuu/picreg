# ============================================================
# Exp-Lasso (picreg family = "gumbel") vs square-root Lasso
# (picreg family = "gaussian") under Gumbel errors.
#
# Reported, for each sparsity level s:
#   - PESR (probability of exact support recovery) for both methods
#   - prediction norm ||X (beta_hat - beta_0)||_2 / sqrt(n) for both methods
#   - Gumbel location / scale estimates (Exp-Lasso only)
# ============================================================

library(picreg)
library(evd)
library(future)
library(future.apply)

# ============================================================
# GLOBAL PARAMETERS
# ============================================================

BETA_VALUE <- 1.0
N_SIM <- 1000L
LAMBDA_ALPHA <- 0.05
LAMBDA_C <- 1.0
LAMBDA_N_SIMU <- 5000L
BASE_SEED <- 2026L

sigma_grid <- 5.0
s_grid <- 0:20

gumbel_noise <- function(n) evd::rgumbel(n, loc = 0, scale = 1)

# picreg keeps the Gumbel profile scale MLE internal; fetch it once.
gumbel_profile <- utils::getFromNamespace("gumbel_profile_cpp", "picreg")


# Same convention as picreg's internal standardization (population sd, n divisor),
# so a lambda computed here is on the scale pic() actually optimizes on.
standardize_design <- function(X) {
  X <- as.matrix(X)
  X <- sweep(X, 2, colMeans(X), "-")
  column_norms <- sqrt(colMeans(X^2))
  column_norms[column_norms == 0] <- 1
  unclass(sweep(X, 2, column_norms, "/"))
}

generate_independent_uniform_design <- function(n, p, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(stats::runif(n * p, min = -sqrt(3), max = sqrt(3)),
              nrow = n, ncol = p)
  # Centered and scaled to unit sample variance at generation.
  standardize_design(X)
}

# ============================================================
# SQUARE-ROOT LASSO PIVOTAL LAMBDA
#
# picreg's "gaussian" family is the square-root Lasso:
#   loss(lp)     = sqrt(mean((y - lp)^2))
#   gradient     = X'(lp - y) / (sqrt(n) * ||y - lp||_2)
#
# so the pivotal penalty level is the (1 - alpha) quantile of
#   max_j |X_j' eps| / (sqrt(n) * ||eps||_2)
#
# under the null. The intercept absorbs the mean, hence the
# centering of eps. Lambda does not depend on sigma.
# ============================================================

lambda_sqrt_lasso_gumbel <- function(X, alpha = 0.05, c = 1.1,
                                     n_simu = 5000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  X <- as.matrix(X)
  n <- nrow(X)
  
  statistics <- replicate(n_simu, {
    epsilon <- gumbel_noise(n)
    epsilon <- epsilon - mean(epsilon)
    max(abs(crossprod(X, epsilon))) / (sqrt(n) * sqrt(sum(epsilon^2)))
  })
  
  c * as.numeric(stats::quantile(statistics,
                                 probs = 1 - alpha,
                                 names = FALSE, type = 8
  ))
}


compute_penalties <- function(X, alpha = 0.05, c = 1.1,
                              n_simu = 5000L, seed = NULL) {
  X_standardized <- standardize_design(X)
  
  lambda_sqrt <- lambda_sqrt_lasso_gumbel(
    X      = X_standardized,
    alpha  = alpha,
    c      = c,
    n_simu = n_simu,
    seed   = if (is.null(seed)) NULL else seed + 1000L
  )
  
  if (!is.null(seed)) set.seed(seed + 2000L)
  
  lambda_gumbel <- picreg::lambda_pdb(
    X_standardized,
    family = "gumbel",
    n_simu = n_simu,
    alpha  = alpha,
    method = "mc_exact"
  )$value
  
  c(gumbel = lambda_gumbel, gaussian = lambda_sqrt)
}


# ============================================================
# SIMULATION FOR ONE (n, p, sigma)
#
# DGP:  y = X beta_0 + sigma * eps,  eps ~ Gumbel(0, 1)
#       => true Gumbel location intercept = 0, true scale = sigma
#
# beta_0: BETA_VALUE on a random support of size s, 0 elsewhere.
#
# X and the penalties are held fixed across s, sigma and repetitions.
# ============================================================

simulate_one_sigma <- function(X, sigma, s_grid, lambda_by_family,
                               n_sim = 1000L, seed = NULL, verbose = TRUE) {
  X <- as.matrix(X)
  sigma <- as.numeric(sigma)
  s_grid <- as.integer(s_grid)
  n_sim <- as.integer(n_sim)
  
  lambda_by_family <- stats::setNames(
    as.numeric(lambda_by_family),
    names(lambda_by_family)
  )
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    stop("sigma must be one positive finite number.")
  }
  if (any(s_grid < 0L) || any(s_grid > p)) {
    stop("All s values must lie in [0, p].")
  }
  
  fitted_families <- c("gumbel", "gaussian")
  beta_value <- BETA_VALUE
  
  all_s_results <- vector("list", length(s_grid))
  
  for (s_id in seq_along(s_grid)) {
    s <- as.integer(s_grid[s_id])
    
    seed_current <- if (is.null(seed)) TRUE else seed + 100000L * s
    
    rep_results <- future.apply::future_lapply(
      seq_len(n_sim),
      function(rep_id) {
        # ------------------------------------------------
        # True beta and noise-free targets
        # ------------------------------------------------
        beta_0 <- numeric(p)
        true_support <- if (s > 0L) {
          sort(sample.int(p, size = s, replace = FALSE))
        } else {
          integer(0L)
        }
        beta_0[true_support] <- beta_value
        
        eta <- as.numeric(X %*% beta_0)
        
        # ------------------------------------------------
        # y = X beta_0 + sigma * Gumbel(0, 1)
        # ------------------------------------------------
        y <- eta + sigma * gumbel_noise(n)
        
        # ------------------------------------------------
        # Fit both methods to exactly the same y
        # ------------------------------------------------
        family_results <- lapply(fitted_families, function(fitted_family) {
          lambda_current <- unname(lambda_by_family[fitted_family])
          
          fit <- picreg::pic(
            X, y,
            family = fitted_family,
            lambda = lambda_current
          )
          
          # coef() returns the ORIGINAL-scale coefficients,
          # row 1 = intercept.
          coefficients <- as.numeric(stats::coef(fit))
          intercept_hat <- coefficients[1L]
          beta_hat <- coefficients[-1L]
          
          # ------------------------------------------------
          # Support recovery
          # ------------------------------------------------
          selected_support <- which(fit$beta != 0)
          
          exact_recovery <- as.integer(
            length(selected_support) == length(true_support) &&
              setequal(selected_support, true_support)
          )
          
          # ------------------------------------------------
          # Prediction norm ||X (beta_hat - beta_0)||_2 / sqrt(n),
          # i.e. the root mean square of the signal error. The
          # 1/sqrt(n) makes it comparable across sample sizes.
          #
          # beta_hat comes from coef(), i.e. the ORIGINAL scale of
          # X, which is the scale beta_0 lives on. The intercept is
          # deliberately absent: the two families estimate different
          # location functionals (the sqrt-Lasso intercept absorbs
          # sigma * Euler-gamma, the Gumbel one estimates the
          # location), so only the slope part is comparable.
          # ------------------------------------------------
          prediction_norm <- sqrt(mean((X %*% (beta_hat - beta_0))^2))
          
          # ------------------------------------------------
          # Gumbel location / scale (Exp-Lasso only)
          #   true location = 0, true scale = sigma
          # ------------------------------------------------
          if (fitted_family == "gumbel") {
            linear_predictor <- stats::predict(fit,
                                               newx = X,
                                               type = "link"
            )
            offset <- y[1L]
            scale_hat <- tryCatch(
              gumbel_profile(
                y - offset,
                linear_predictor - offset
              )$sigma,
              error = function(e) NA_real_
            )
            loc_hat <- intercept_hat
          } else {
            scale_hat <- NA_real_
            loc_hat <- NA_real_
          }
          
          data.frame(
            fitted_family    = fitted_family,
            lambda           = lambda_current,
            exact_recovery   = exact_recovery,
            prediction_norm  = prediction_norm,
            loc_hat          = loc_hat,
            scale_hat        = as.numeric(scale_hat),
            stringsAsFactors = FALSE
          )
        })
        
        do.call(rbind, family_results)
      },
      future.seed = seed_current,
      future.packages = c("picreg", "evd"),
      future.globals = list(
        X                = X,
        n                = n,
        p                = p,
        s                = s,
        sigma            = sigma,
        beta_value       = beta_value,
        fitted_families  = fitted_families,
        lambda_by_family = lambda_by_family,
        gumbel_noise     = gumbel_noise,
        gumbel_profile   = gumbel_profile
      )
    )
    
    rep_results <- do.call(rbind, rep_results)
    
    # ------------------------------------------------------
    # Aggregate over repetitions
    # ------------------------------------------------------
    finite_mean <- function(x) {
      x <- x[is.finite(x)]
      if (length(x)) mean(x) else NA_real_
    }
    
    summary_by_family <- do.call(rbind, lapply(
      split(rep_results, rep_results$fitted_family),
      function(z) {
        pesr <- mean(z$exact_recovery)
        
        data.frame(
          n                = n,
          p                = p,
          log_p_over_n     = log(p / n),
          sigma            = sigma,
          s                = s,
          n_sim            = n_sim,
          fitted_family    = z$fitted_family[1L],
          lambda           = z$lambda[1L],
          pesr             = pesr,
          pesr_se          = sqrt(pesr * (1 - pesr) / n_sim),
          prediction_norm  = finite_mean(z$prediction_norm),
          loc_hat          = finite_mean(z$loc_hat),
          scale_hat        = finite_mean(z$scale_hat),
          stringsAsFactors = FALSE
        )
      }
    ))
    
    all_s_results[[s_id]] <- summary_by_family
    
    if (verbose) {
      pick <- function(column, family) {
        summary_by_family[[column]][summary_by_family$fitted_family == family]
      }
      
      message(sprintf(
        paste0(
          "[sigma=%.1f | n=%d | p=%d | s=%2d] ",
          "PESR: Exp=%.3f sqrt=%.3f | ",
          "PredNorm: Exp=%.3f sqrt=%.3f | ",
          "Gumbel loc=%.3f (true 0) scale=%.3f (true %.3f)"
        ),
        sigma, n, p, s,
        pick("pesr", "gumbel"), pick("pesr", "gaussian"),
        pick("prediction_norm", "gumbel"), pick("prediction_norm", "gaussian"),
        pick("loc_hat", "gumbel"), pick("scale_hat", "gumbel"), sigma
      ))
    }
  }
  
  result <- do.call(rbind, all_s_results)
  rownames(result) <- NULL
  result
}


# ============================================================
# EXPERIMENTAL GRID
#
#   n = 450, log(p/n) = -2
#   n = 600, log(p/n) = -1.5
#   n = 750, log(p/n) = -1
# ============================================================

simulation_grid <- data.frame(
  n                   = c(600L, 800L, 1000L),
  log_p_over_n_target = c(-2, -1.5, -1)
)

simulation_grid$p <- as.integer(round(
  simulation_grid$n * exp(simulation_grid$log_p_over_n_target)
))

print(simulation_grid)


# ============================================================
# FIXED DESIGN AND PENALTIES
#
# The design is drawn iid Uniform[-sqrt(3), sqrt(3)] and reused for
# every sparsity level, repetition and sigma.
# ============================================================

designs <- vector("list", nrow(simulation_grid))
penalties <- vector("list", nrow(simulation_grid))

for (k in seq_len(nrow(simulation_grid))) {
  n_current <- simulation_grid$n[k]
  p_current <- simulation_grid$p[k]
  
  designs[[k]] <- generate_independent_uniform_design(
    n = n_current, p = p_current, seed = BASE_SEED + 10000L * k
  )
  
  penalties[[k]] <- compute_penalties(
    X      = designs[[k]],
    alpha  = LAMBDA_ALPHA,
    c      = LAMBDA_C,
    n_simu = LAMBDA_N_SIMU,
    seed   = BASE_SEED + 20000L * k
  )
  
  message(sprintf(
    "[n=%d, p=%d] lambda Exp-Lasso = %.6f | lambda sqrt-Lasso = %.6f",
    n_current, p_current,
    penalties[[k]]["gumbel"], penalties[[k]]["gaussian"]
  ))
}


# ============================================================
# RUN
# ============================================================

old_plan <- future::plan()

# Cap at the cores parallelly actually reports, otherwise it warns (soft
# limit) or refuses outright (hard limit) on smaller machines.
workers <- 12L#max(1L, min(12L, future::availableCores()))
future::plan(future::multisession, workers = workers)
message(sprintf("Parallel workers: %d", workers))

# Restoring the caller's plan can itself fail: the captured plan may ask for
# more workers than parallelly now allows (e.g. a 15-worker multisession left
# over from an earlier session). Never let teardown abort the run.
restore_plan <- function(previous_plan) {
  tryCatch(
    future::plan(previous_plan),
    error = function(e) {
      warning("Could not restore the previous future plan (",
              conditionMessage(e), "). Falling back to sequential.",
              call. = FALSE
      )
      future::plan(future::sequential)
    }
  )
}

get_output_directory <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    
    if (length(active_path) == 1L && nzchar(active_path)) {
      return(dirname(active_path))
    }
  }
  getwd()
}

output_dir <- get_output_directory()
message(sprintf("CSV output directory: %s", output_dir))

results <- list()

for (sigma_current in sigma_grid) {
  for (k in seq_len(nrow(simulation_grid))) {
    block <- simulate_one_sigma(
      X                = designs[[k]],
      sigma            = sigma_current,
      s_grid           = s_grid,
      lambda_by_family = penalties[[k]],
      n_sim            = N_SIM,
      # Same seed across sigma: identical Gumbel draws are reused.
      seed             = BASE_SEED + 30000L * k,
      verbose          = TRUE
    )
    block$log_p_over_n_target <- simulation_grid$log_p_over_n_target[k]
    results[[length(results) + 1L]] <- block
  }
}

results <- do.call(rbind, results)
rownames(results) <- NULL

results <- results[, c(
  "n", "p", "log_p_over_n_target", "log_p_over_n", "sigma", "s", "n_sim",
  "fitted_family", "lambda",
  "pesr", "pesr_se",
  "prediction_norm",
  "loc_hat", "scale_hat"
)]

output_file <- file.path(output_dir, "PESR_Gumbel.csv")
write.csv(results, file = output_file, row.names = FALSE)

message(sprintf("Saved: %s", output_file))

# Teardown only after the results are safely on disk.
restore_plan(old_plan)

print(results)
