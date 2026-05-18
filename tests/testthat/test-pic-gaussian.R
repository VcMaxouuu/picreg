test_that("pic() fits a sparse Gaussian model and recovers the support", {
  set.seed(1)
  n <- 80; p <- 60; s <- 3
  X <- matrix(rnorm(n * p), n, p)
  true_features <- sort(sample.int(p, s))
  beta <- numeric(p); beta[true_features] <- 3
  y <- as.numeric(X %*% beta + rnorm(n))

  fit <- pic(X, y, family = "gaussian", penalty = "lasso",
             lambda_n_simu = 500L)

  expect_s3_class(fit, "pic")
  expect_s3_class(fit, "pic.gaussian")
  expect_length(fit$beta, p)
  expect_true(fit$lambda > 0)
  expect_true(all(true_features %in% fit$selected))
})

test_that("pic() returns the right object structure", {
  set.seed(2)
  n <- 50; p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X[, 1] + rnorm(n))

  fit <- pic(X, y, family = "gaussian", lambda_n_simu = 200L)

  expect_true(all(c("beta", "intercept", "df", "selected", "lambda",
                    "family", "penalty", "lambda_pdb",
                    "n_samples", "n_features") %in% names(fit)))
  expect_equal(fit$n_samples, n)
  expect_equal(fit$n_features, p)
  expect_equal(fit$df, length(fit$selected))
})

test_that("predict.pic returns sensible outputs for type = link / response", {
  set.seed(3)
  n <- 60; p <- 15
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X[, 1] + rnorm(n))

  fit <- pic(X, y, family = "gaussian", lambda_n_simu = 200L)

  link <- predict(fit, X, type = "link")
  resp <- predict(fit, X, type = "response")

  expect_length(link, n)
  expect_length(resp, n)
  expect_equal(link, resp)   # Gaussian: identity link
})

test_that("coef.pic returns a named vector with intercept first", {
  set.seed(4)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X[, 1] + rnorm(n))

  fit <- pic(X, y, family = "gaussian", lambda_n_simu = 200L)
  co  <- coef(fit)

  expect_length(co, p + 1L)
  expect_equal(names(co)[1L], "(Intercept)")
})
