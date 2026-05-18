devtools::load_all()   # run from the project root (picreg/)
library(evd)

# -------------------------------------------------------------------
# Simulation settings
# -------------------------------------------------------------------

m <- 50
n <- 500
p <- 200
s <- 5

simulate_cox_data <- function(n, p, s, h0 = 1, target_censoring = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  
  beta <- numeric(p)
  true_features <- sample.int(p, s, replace = FALSE)
  beta[true_features] <- 3 * sample(c(-1, 1), s, replace = TRUE)
  
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
    time = time,
    event = event,
    y = cbind(time = time, event = event),
    true_features = sort(true_features),
    beta = beta
  )
}

correct_selection <- logical(m)
t0 <- Sys.time()

# -------------------------------------------------------------------
# Monte Carlo loop
# -------------------------------------------------------------------

for (k in seq_len(m)) {
  
  # Generate design
  data <- simulate_cox_data(n, p, s)
  X <- data$X
  y <- data$y
  true_features <- data$true_features
  
  # Fit model
  fit <- pic(
    X,
    y,
    family = "cox",
    penalty = "lasso",
    lambda = 0.0355
  )
  
  # Compare supports
  selected <- sort(fit$selected)
  
  correct_selection[k] <- identical(selected, true_features)
}

t1 <- Sys.time()

# -------------------------------------------------------------------
# Results
# -------------------------------------------------------------------

mean_recovery <- mean(correct_selection)
elapsed <- as.numeric(t1 - t0, units = "secs")

cat("\n")
cat(sprintf(
  "Exact support recovery rate: %.3f\n",
  mean_recovery
))

cat(sprintf(
  "Elapsed time: %.2f seconds\n",
  elapsed
))