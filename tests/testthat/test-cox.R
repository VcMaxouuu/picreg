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

test_that("predict(type='survival') returns time and survival of the right shape", {
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
  sf <- predict(fit, newx, type = "survival")
  expect_true(is.list(sf))
  expect_true(all(c("time", "survival") %in% names(sf)))
  expect_equal(ncol(sf$survival), 5L)
  expect_equal(nrow(sf$survival), length(sf$time))
})

test_that("concordance_index matches a brute-force count, ties included", {
  set.seed(44)
  n <- 150
  times <- round(rexp(n), 1)
  events <- rbinom(n, 1, 0.6)
  predictions <- round(rnorm(n), 1)

  num <- 0; den <- 0
  for (i in which(events == 1)) {
    for (j in seq_len(n)) {
      if (i != j && (times[i] < times[j] ||
                     (times[i] == times[j] && events[j] == 0))) {
        den <- den + 1
        num <- num + (predictions[i] > predictions[j]) +
          0.5 * (predictions[i] == predictions[j])
      }
    }
  }
  expect_equal(concordance_index(times, events, predictions), num / den)
  expect_equal(concordance_index(times, rep(0, n), predictions), 0)
  expect_error(concordance_index(times, events, predictions[-1]),
               "same length")
})

test_that("cox_partial_log_likelihood matches the Breslow definition", {
  set.seed(45)
  n <- 120
  times <- round(rexp(n), 1)
  events <- rbinom(n, 1, 0.6)
  predictions <- rnorm(n)

  brute <- -sum(vapply(which(events == 1), function(i) {
    predictions[i] - log(sum(exp(predictions[times >= times[i]])))
  }, numeric(1L))) / n
  expect_equal(cox_partial_log_likelihood(times, events, predictions), brute)
})

test_that("Cox predict(type='response') is the relative risk; coef() has no intercept", {
  data(CoxExample)
  fit <- pic(CoxExample$X, CoxExample$y, family = "cox", lambda = 0.05)
  newx <- CoxExample$X[1:4, ]
  expect_equal(predict(fit, newx, type = "response"),
               exp(predict(fit, newx, type = "link")))

  co <- coef(fit)
  expect_equal(dim(co), c(ncol(CoxExample$X), 1L))
  expect_false("(Intercept)" %in% rownames(co))
  expect_equal(rownames(co), colnames(CoxExample$X))
})
