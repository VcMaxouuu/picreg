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

test_that("Gumbel conditional scale is the profile MLE", {
  set.seed(14)
  residual_cases <- list(
    typical  = rnorm(100) - log(-log(runif(100))),
    skewed   = c(-100, rep(1, 49)),
    constant = rep(3, 30),
    tiny     = c(-1e-10, 0, 1e-10)
  )

  for (r in residual_cases) {
    out <- picreg:::gumbel_profile_cpp(r, numeric(length(r)))
    score <- -length(r) * out$sigma +
      sum(r * (1 - exp(-r / out$sigma)))
    expect_lt(abs(score) / (length(r) * out$sigma), 1e-10)

    objective <- function(log_sigma) {
      sigma <- exp(log_sigma)
      z <- r / sigma
      log(sigma) + mean(z + exp(-z))
    }
    scale_ref <- max(abs(r))
    reference <- optimize(objective,
                          log(scale_ref) + c(-10, 10),
                          tol = 1e-10)
    expect_equal(out$sigma, exp(reference$minimum), tolerance = 1e-7)
  }
})

test_that("Gumbel profiled residual is the numerical gradient", {
  set.seed(15)
  n <- 30
  y <- 1.5 + 2 * (-log(-log(runif(n))))
  lp <- rnorm(n)
  out <- picreg:::gumbel_profile_cpp(y, lp)
  eps <- 1e-6

  finite_difference <- vapply(seq_len(n), function(i) {
    lp_plus <- lp_minus <- lp
    lp_plus[i] <- lp_plus[i] + eps
    lp_minus[i] <- lp_minus[i] - eps
    (picreg:::gumbel_profile_cpp(y, lp_plus)$value -
       picreg:::gumbel_profile_cpp(y, lp_minus)$value) / (2 * eps)
  }, numeric(1))

  expect_equal(as.numeric(out$residual), finite_difference, tolerance = 1e-6)
  expect_equal(out$value, exp(out$raw_loss), tolerance = 1e-13)
})

test_that("Gumbel null location and scale are joint MLEs", {
  set.seed(16)
  y <- 4 + 2.5 * (-log(-log(runif(80))))
  out <- picreg:::gumbel_null_fit_cpp(y)
  z <- (y - out$intercept) / out$sigma

  location_score <- mean(exp(-z) - 1) / out$sigma
  scale_score <- -length(y) * out$sigma +
    sum((y - out$intercept) * (1 - exp(-z)))
  expect_lt(abs(location_score), 1e-11)
  expect_lt(abs(scale_score) / (length(y) * out$sigma), 1e-11)
  expect_lt(abs(sum(out$residual)), 1e-11)

  # Location/scale equivariance is essential to the pivotal construction.
  shifted <- picreg:::gumbel_null_fit_cpp(-7 + 3 * y)
  expect_equal(shifted$intercept, -7 + 3 * out$intercept, tolerance = 1e-10)
  expect_equal(shifted$sigma, 3 * out$sigma, tolerance = 1e-10)
  expect_equal(as.numeric(shifted$residual), as.numeric(out$residual),
               tolerance = 1e-10)
})

test_that("Gumbel lambda_max uses the joint null MLE", {
  set.seed(17)
  n <- 70; p <- 12
  X <- scale(matrix(rnorm(n * p), n, p), center = TRUE, scale = FALSE)
  y <- 2 - 1.3 * log(-log(runif(n)))
  null <- picreg:::gumbel_null_fit_cpp(y)

  expected <- max(abs(crossprod(X, null$residual)))
  observed <- picreg:::lambda_max_cpp(X, y, "gumbel", TRUE)
  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("Gumbel reports an undefined degenerate MLE", {
  expect_error(
    picreg:::gumbel_profile_cpp(1:5, 1:5),
    "MLE is undefined when all residuals are zero"
  )
})

test_that("Gumbel fitting is stable under a large response shift", {
  set.seed(18)
  n <- 80; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(0.8 * X[, 1] - 0.5 * X[, 2] - log(-log(runif(n))))

  fit <- pic(X, y, family = "gumbel", lambda = 0.08,
             path_length = 3L, tol = 1e-7)
  shifted <- pic(X, y + 1e10, family = "gumbel", lambda = 0.08,
                 path_length = 3L, tol = 1e-7)

  expect_equal(shifted$beta, fit$beta, tolerance = 2e-5)
  expect_equal(shifted$intercept - 1e10, fit$intercept,
               tolerance = 2e-5)
})
