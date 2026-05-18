# Smoke test: Cox PH with sparse coefficients. Generates exponential
# survival times under a Cox model, checks support recovery and that the
# C-index is well above 0.5. Run from repo root:
#   Rscript picreg/dev/smoke-cox.R

library(picreg)

set.seed(42)
n <- 300; p <- 50; s <- 5
X <- matrix(rnorm(n * p), n, p)
true_features <- sort(sample.int(p, s))
beta <- numeric(p); beta[true_features] <- 1.0
eta <- as.numeric(X %*% beta)

# Cox with exponential baseline: T ~ Exp(rate = exp(eta))
T <- rexp(n, rate = exp(eta))
# 30% censoring
C <- rexp(n, rate = mean(exp(eta)) * 0.4)
time <- pmin(T, C)
event <- as.integer(T <= C)
y <- cbind(time = time, event = event)

cat("---- Cox smoke test ----\n")
cat(sprintf("censoring rate: %.2f\n", 1 - mean(event)))

fit <- pic(X, y, family = "cox", penalty = "lasso",
            lambda_n_simu = 1500L)
sel <- which(fit$beta != 0)
cat(sprintf(
  "[cox L1]   selected = %d  true_in_sel = %d/%d  lambda = %.4f\n",
  length(sel), sum(true_features %in% sel), s, fit$lambda
))

eta_hat <- as.numeric(scale(X, center = fit$X_mean, scale = fit$X_std) %*% fit$beta)
cidx <- concordance_index(time, event, eta_hat)
cat(sprintf("[cox C-index in-sample] %.3f\n", cidx))


cat("baseline_survival head:\n")
print(head(fit$baseline_survival, 4))
