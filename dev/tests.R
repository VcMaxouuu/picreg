devtools::load_all()   # run from the project root (picreg/)
library(evd)

# -------------------------------------------------------------------
# Simulation settings
# -------------------------------------------------------------------

m <- 50
n <- 500
p <- 200
s <- 5

correct_selection <- logical(m)
t0 <- Sys.time()

# -------------------------------------------------------------------
# Monte Carlo loop
# -------------------------------------------------------------------

for (k in seq_len(m)) {
  
  # Generate design
  X <- matrix(rnorm(n * p), n, p)
  
  # True support
  true_features <- sort(sample.int(p, s))
  
  beta <- numeric(p)
  beta[true_features] <- 3
  
  # Generate response
  noise <- rnorm(n)
  y <- as.numeric(X %*% beta + noise)
  
  # Fit model
  fit <- pic(
    X,
    y,
    family = "gaussian",
    penalty = "scad"
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