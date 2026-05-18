test_that("pic() fits a Cox model and produces baseline functions", {
  set.seed(41)
  n <- 100; p <- 20
  X <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p); beta[1L:2L] <- 1
  eta  <- as.numeric(X %*% beta)
  t_event <- rexp(n, rate = exp(eta))
  t_cens  <- rexp(n, rate = 1)
  y <- cbind(time  = pmin(t_event, t_cens),
             event = as.integer(t_event <= t_cens))

  fit <- pic(X, y, family = "cox", penalty = "lasso",
             lambda_n_simu = 200L)

  expect_s3_class(fit, "pic.cox")
  expect_null(fit$intercept)
  expect_true(!is.null(fit$baseline_cumulative_hazard))
  expect_true(!is.null(fit$baseline_survival))
})

test_that("concordance_index returns a scalar in [0, 1]", {
  set.seed(42)
  n <- 60
  times       <- rexp(n, 1)
  events      <- rbinom(n, 1, 0.7)
  predictions <- rnorm(n)
  ci <- concordance_index(times, events, predictions)
  expect_length(ci, 1L)
  expect_true(ci >= 0 && ci <= 1)
})

test_that("predict_survival_function returns time and survival of the right shape", {
  set.seed(43)
  n <- 80; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  eta <- as.numeric(X[, 1L])
  t_event <- rexp(n, rate = exp(eta))
  t_cens  <- rexp(n, rate = 1)
  y <- cbind(pmin(t_event, t_cens),
             as.integer(t_event <= t_cens))

  fit <- pic(X, y, family = "cox", lambda_n_simu = 200L)

  newx <- X[1L:5L, , drop = FALSE]
  sf <- predict_survival_function(fit, newx)
  expect_true(is.list(sf))
  expect_true(all(c("time", "survival") %in% names(sf)))
  expect_equal(ncol(sf$survival), 5L)
  expect_equal(nrow(sf$survival), length(sf$time))
})
