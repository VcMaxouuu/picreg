# Smoke test: SCAD recovery rate on a simple Gaussian sparse-regression problem.
# Run from the repository root after `R CMD INSTALL picreg`:
#   Rscript picreg/dev/smoke-scad.R
library(picreg)

run_one <- function(s, p = 100, n = 100, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  X <- matrix(rnorm(n * p), n, p)
  true_features <- sample.int(p, s)
  
  beta <- numeric(p)
  beta[true_features] <- 3
  
  y <- as.numeric(X %*% beta + rnorm(n))
  
  fit <- pic(X, y, family = "gaussian", penalty = "scad")
  
  selected <- which(fit$beta != 0)
  recovery <- setequal(true_features, selected)
  
  list(
    s = s,
    recovery = recovery,
    n_selected = length(selected),
    n_true_in_selected = sum(true_features %in% selected),
    lambda = fit$lambda
  )
}

cat("---- SCAD smoke test ----\n")

s <- 4
n_sim <- 1

results <- lapply(seq_len(n_sim), function(i) {
  run_one(s = s, seed = i)
})

recoveries <- vapply(results, function(res) res$recovery, logical(1))
n_selected <- vapply(results, function(res) res$n_selected, integer(1))
n_true_in_selected <- vapply(results, function(res) res$n_true_in_selected, integer(1))

cat(sprintf("s = %d\n", s))
cat(sprintf("n_sim = %d\n", n_sim))
cat(sprintf(
  "Exact recovery rate = %.2f%% (%d / %d)\n",
  100 * mean(recoveries),
  sum(recoveries),
  n_sim
))
cat(sprintf(
  "Average number of selected variables = %.2f\n",
  mean(n_selected)
))
cat(sprintf(
  "Average number of true variables selected = %.2f / %d\n",
  mean(n_true_in_selected),
  s
))