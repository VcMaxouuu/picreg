test_that("all non-Cox families fit without error on small data", {
  set.seed(11)
  n <- 60; p <- 20
  X <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p); beta[1L:2L] <- 2
  eta  <- as.numeric(X %*% beta)

  responses <- list(
    gaussian    = eta + rnorm(n),
    binomial    = rbinom(n, 1, 1 / (1 + exp(-eta))),
    poisson     = rpois(n, lambda = exp(pmin(eta, 5))),
    exponential = rexp(n, rate = exp(pmin(eta, 5))),
    gumbel      = eta - log(-log(runif(n)))
  )

  for (fam in names(responses)) {
    fit <- pic(X, responses[[fam]], family = fam, penalty = "lasso",
               lambda_n_simu = 200L)
    expect_s3_class(fit, paste0("pic.", fam))
    expect_true(fit$lambda > 0,
                label = sprintf("lambda > 0 for family = %s", fam))
  }
})

test_that("get_family() rejects unknown family names", {
  expect_error(get_family("normal"))
  expect_error(get_family("weibull"))
  expect_error(get_family(123))
})

test_that("binomial predict(type='class') returns 0/1 integers", {
  set.seed(12)
  n <- 60; p <- 15
  X <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p); beta[1L] <- 2
  prob <- 1 / (1 + exp(-as.numeric(X %*% beta)))
  y    <- rbinom(n, 1, prob)

  fit <- pic(X, y, family = "binomial", lambda_n_simu = 200L)
  cls <- predict(fit, X, type = "class")

  expect_true(all(cls %in% c(0L, 1L)))
  expect_length(cls, n)
})

test_that("predict(type='class') errors on non-binomial fits", {
  set.seed(13)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X[, 1] + rnorm(n))

  fit <- pic(X, y, family = "gaussian", lambda_n_simu = 200L)
  expect_error(predict(fit, X, type = "class"), "binomial")
})
