library(picreg)
library(evd)
library(future)
library(future.apply)

P_FIXED <- 100L
S_FIXED <- 4L
SIGMA <- 5
BETA_VALUE <- 1
N_SIM <- 1000L
LAMBDA_ALPHA <- 0.05
LAMBDA_N_SIMU <- 2000L
BASE_SEED <- 2026L

n_grid <- seq(200L, 4000L, by = 200L)

gumbel_noise <- function(n) evd::rgumbel(n, loc = 0, scale = 1)
gumbel_profile <- utils::getFromNamespace("gumbel_profile_cpp", "picreg")

standardize_design <- function(X) {
  X <- sweep(X, 2, colMeans(X), "-")
  norms <- sqrt(colMeans(X^2))
  norms[norms == 0] <- 1
  unclass(sweep(X, 2, norms, "/"))
}

generate_design <- function(n, p, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(runif(n * p, -sqrt(3), sqrt(3)), nrow = n, ncol = p)
  standardize_design(X)
}

oracle_gumbel_mle <- function(y, X, beta_true) {
  z <- y - as.numeric(X %*% beta_true)
  fit <- tryCatch(evd::fgev(z, shape = 0), error = function(e) NULL)
  if (is.null(fit)) return(c(loc = NA_real_, scale = NA_real_))
  c(loc = unname(fit$estimate["loc"]),
    scale = unname(fit$estimate["scale"]))
}

simulate_one_n <- function(n, p, s, sigma, lambda, n_sim, X, seed = NULL) {
  
  reps <- future_lapply(seq_len(n_sim), function(rep_id) {
    beta_0 <- numeric(p)
    true_support <- sort(sample.int(p, s))
    beta_0[true_support] <- BETA_VALUE
    
    y <- as.numeric(X %*% beta_0) + sigma * gumbel_noise(n)
    
    # Exp-Lasso
    fit <- picreg::pic(X, y, family = "gumbel", lambda = lambda)
    loc_hat <- as.numeric(coef(fit))[1L]
    
    selected_support <- which(fit$beta != 0)
    exact_recovery <- as.integer(setequal(selected_support, true_support))
    
    eta_hat <- predict(fit, newx = X, type = "link")
    offset <- y[1L]
    scale_hat <- tryCatch(
      gumbel_profile(y - offset, eta_hat - offset)$sigma,
      error = function(e) NA_real_
    )
    
    # Oracle MLE using true beta
    oracle <- oracle_gumbel_mle(y, X, beta_0)
    
    data.frame(
      exact_recovery = exact_recovery,
      loc_hat = loc_hat,
      scale_hat = as.numeric(scale_hat),
      mle_loc_hat = oracle["loc"],
      mle_scale_hat = oracle["scale"]
    )
  },
  future.seed = if (is.null(seed)) TRUE else seed,
  future.packages = c("picreg", "evd"),
  future.globals = list(
    X = X, n = n, p = p, s = s, sigma = sigma, lambda = lambda,
    gumbel_noise = gumbel_noise, gumbel_profile = gumbel_profile,
    oracle_gumbel_mle = oracle_gumbel_mle, BETA_VALUE = BETA_VALUE
  ))
  
  reps <- do.call(rbind, reps)
  fmean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) mean(x) else NA_real_
  }
  # Monte Carlo variance of the estimator across repetitions. Expected to
  # decay like 1/n, so n * var should plateau at the asymptotic variance.
  fvar <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) >= 2L) stats::var(x) else NA_real_
  }
  
  data.frame(
    n = n,
    pesr = mean(reps$exact_recovery),
    mean_loc_hat = fmean(reps$loc_hat),
    mean_scale_hat = fmean(reps$scale_hat),
    mean_mle_loc_hat = fmean(reps$mle_loc_hat),
    mean_mle_scale_hat = fmean(reps$mle_scale_hat),
    var_loc_hat = fvar(reps$loc_hat),
    var_scale_hat = fvar(reps$scale_hat),
    var_mle_loc_hat = fvar(reps$mle_loc_hat),
    var_mle_scale_hat = fvar(reps$mle_scale_hat)
  )
}

old_plan <- future::plan()
future::plan(future::multisession, workers = 12L)

results <- vector("list", length(n_grid))

for (k in seq_along(n_grid)) {
  n <- n_grid[k]
  
  X <- generate_design(n, P_FIXED, BASE_SEED + 10000L * k)
  
  set.seed(BASE_SEED + 20000L * k)
  lambda <- picreg::lambda_pdb(
    standardize_design(X),
    family = "gumbel",
    n_simu = LAMBDA_N_SIMU,
    alpha = LAMBDA_ALPHA,
    method = "mc_exact"
  )$value
  
  results[[k]] <- simulate_one_n(
    n, P_FIXED, S_FIXED, SIGMA, lambda, N_SIM, X,
    BASE_SEED + 30000L * k
  )
  
  r <- results[[k]]
  message(sprintf(
    paste0("[n=%5d | lambda=%.6f] PESR=%.3f | ",
           "Lasso: loc=%+.4f scale=%.4f (var %.2e / %.2e) | ",
           "MLE: loc=%+.4f scale=%.4f (var %.2e / %.2e)"),
    n, lambda, r$pesr,
    r$mean_loc_hat, r$mean_scale_hat, r$var_loc_hat, r$var_scale_hat,
    r$mean_mle_loc_hat, r$mean_mle_scale_hat,
    r$var_mle_loc_hat, r$var_mle_scale_hat
  ))
}

results <- do.call(rbind, results)
rownames(results) <- NULL

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
output_file <- file.path(output_dir, "asymptotic_convergence_gumbel.csv")
write.csv(results, output_file, row.names = FALSE)

tryCatch(
  future::plan(old_plan),
  error = function(e) future::plan(future::sequential)
)

print(results)