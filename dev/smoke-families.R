# Quick smoke check across families. Run from the repo root:
#   Rscript picreg/dev/smoke-families.R

library(picreg)

set.seed(7)
n <- 200; p <- 50; s <- 5
X <- matrix(rnorm(n * p), n, p)
true_features <- sort(sample.int(p, s))
beta <- numeric(p); beta[true_features] <- 1.5
eta <- as.numeric(X %*% beta)

cat("---- Family smoke checks ----\n")

# Gaussian
y_g <- eta + rnorm(n)
fg <- pic(X, y_g, family = "gaussian", penalty = "lasso",
           lambda_n_simu = 1500L)
cat(sprintf(
  "[gaussian L1]   selected = %d  true_in_sel = %d/%d  lambda = %.4f\n",
  sum(fg$beta != 0), sum(true_features %in% which(fg$beta != 0)), s, fg$lambda
))

# Binomial
p_y <- 1 / (1 + exp(-eta))
y_b <- rbinom(n, 1, p_y)
fb <- pic(X, y_b, family = "binomial", penalty = "lasso",
           lambda_n_simu = 1500L)
cat(sprintf(
  "[binomial L1]   selected = %d  true_in_sel = %d/%d  lambda = %.4f\n",
  sum(fb$beta != 0), sum(true_features %in% which(fb$beta != 0)), s, fb$lambda
))

# Poisson
mu_p <- exp(0.5 * eta + 0.5)
y_p <- rpois(n, mu_p)
fp <- pic(X, y_p, family = "poisson", penalty = "lasso",
           lambda_n_simu = 1500L)
cat(sprintf(
  "[poisson  L1]   selected = %d  true_in_sel = %d/%d  lambda = %.4f\n",
  sum(fp$beta != 0), sum(true_features %in% which(fp$beta != 0)), s, fp$lambda
))

# Predict sanity check
yhat_g <- predict(fg, newx = X, type = "response")
r2 <- 1 - sum((y_g - yhat_g) ^ 2) / sum((y_g - mean(y_g)) ^ 2)
mse <- mean((y_g - yhat_g) ^ 2)
cat(sprintf("[gaussian]  R2 = %.3f  MSE = %.3f\n", r2, mse))

yhat_b <- predict(fb, newx = X, type = "class")
acc <- mean(yhat_b == y_b)
cat(sprintf("[binomial]  accuracy = %.3f\n", acc))
