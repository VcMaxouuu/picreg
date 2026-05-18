# Smoke test: Cox recovery rate with fixed lambda.
# Run from repository root:
#   Rscript picreg/dev/smoke-gumbel.R

library(picreg)

simulate_cox_data <- function(n, p, s, X = NULL, h0 = 1,
                              target_censoring = 0.3,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  if (is.null(X)) {
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  }
  
  beta <- numeric(p)
  true_features <- sample.int(p, size = s, replace = FALSE)
  beta[true_features] <- 3 * sample(c(-1, 1), size = s, replace = TRUE)
  
  eta <- as.numeric(X %*% beta)
  event_rate <- h0 * exp(eta)
  
  U <- runif(n)
  T <- -log(U) / event_rate
  
  censoring_error <- function(lambda_c) {
    mean(1 - exp(-lambda_c * T)) - target_censoring
  }
  
  lambda_c <- uniroot(
    censoring_error,
    interval = c(1e-12, 1e12)
  )$root
  
  C <- rexp(n, rate = lambda_c)
  
  time <- pmin(T, C)
  event <- as.integer(T <= C)
  
  list(
    X = X,
    y = time,
    delta = event,
    true_features = true_features,
    beta = beta
  )
}

cat("---- Cox smoke test with fixed lambda ----\n")

n <- 500
p <- 75
s <- 24
n_sim <- 100
target_censoring <- 0.3

# First run: estimate lambda once with SCAD
data0 <- simulate_cox_data(
  n = n,
  p = p,
  s = s,
  target_censoring = target_censoring,
  seed = 1
)

X_fixed <- data0$X
Y0 <- cbind(time = data0$y, event = data0$delta)

model0 <- pic(
  X_fixed,
  Y0,
  family = "cox",
  penalty = "scad",
  lambda_method = "mc_exact"
)

lambda_fixed <- model0$lambda

cat(sprintf("Fixed lambda = %.6f\n", lambda_fixed))

recoveries <- logical(n_sim)
event_rates <- numeric(n_sim)
censoring_rates <- numeric(n_sim)
n_selected <- integer(n_sim)
n_true_in_selected <- integer(n_sim)

for (i in seq_len(n_sim)) {
  dat <- simulate_cox_data(
    n = n,
    p = p,
    s = s,
    X = X_fixed,
    target_censoring = target_censoring,
    seed = i + 1000
  )
  
  Y <- cbind(time = dat$y, event = dat$delta)
  
  fit <- pic(
    X_fixed,
    Y,
    family = "cox",
    penalty = "lasso",
    lambda = lambda_fixed
  )
  
  selected <- which(fit$beta != 0)
  
  recoveries[i] <- setequal(selected, dat$true_features)
  event_rates[i] <- mean(dat$delta)
  censoring_rates[i] <- mean(dat$delta == 0)
  n_selected[i] <- length(selected)
  n_true_in_selected[i] <- sum(dat$true_features %in% selected)
}

cat(sprintf("n = %d, p = %d, s = %d\n", n, p, s))
cat(sprintf("n_sim = %d\n", n_sim))

cat(sprintf(
  "Exact recovery rate = %.2f%% (%d / %d)\n",
  100 * mean(recoveries),
  sum(recoveries),
  n_sim
))

cat(sprintf(
  "Average event rate = %.3f\n",
  mean(event_rates)
))