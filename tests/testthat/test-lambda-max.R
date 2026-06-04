# lambda_max — closed-form KKT threshold ||grad L(0, beta0*)||_inf.
# Used to seed the warm-start path; tests its consistency with the
# behavior of `pic()` at lambda = lambda_max.

test_that("lambda_max scales as ||X^T r|| / n for the Gaussian family", {
  set.seed(100)
  n <- 80; p <- 30
  X <- matrix(rnorm(n * p), n, p)
  mu <- colMeans(X)
  X  <- sweep(X, 2L, mu, "-")
  sv <- sqrt(colMeans(X^2)); sv[sv == 0] <- 1
  X  <- sweep(X, 2L, sv, "/")

  y <- as.numeric(X[, 1L] + rnorm(n))

  lmax <- picreg:::.lambda_max(X, y, picreg:::get_family("gaussian"),
                               fit_intercept = TRUE)
  expect_true(is.numeric(lmax))
  expect_true(lmax > 0)
  expect_length(lmax, 1L)
})

test_that("pic() returns the empty support at lambda just above lambda_max", {
  set.seed(101)
  n <- 60; p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X[, 1L] + rnorm(n))

  prep <- picreg:::check_Xy(X, y, y_kind = "continuous", standardize_X = TRUE)
  lmax <- picreg:::.lambda_max(prep$X, prep$y,
                               picreg:::get_family("gaussian"),
                               fit_intercept = TRUE)

  fit_above <- pic(X, y, lambda = lmax * 1.05)
  expect_equal(sum(fit_above$beta != 0), 0L)
})

test_that("path_length governs the warm-start path granularity", {
  set.seed(102)
  n <- 80; p <- 30
  X <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p); beta[1:3] <- 2
  y <- as.numeric(X %*% beta + rnorm(n))

  # Short and long paths should agree on the final solution (the path is
  # only a speed-up, not an accuracy lever, once it converges).
  fit_short <- pic(X, y, path_length = 2L,  lambda_n_simu = 200L)
  fit_long  <- pic(X, y, path_length = 20L, lambda_n_simu = 200L)

  expect_setequal(fit_short$selected, fit_long$selected)
  expect_lt(max(abs(fit_short$beta - fit_long$beta)), 1e-2)
})
