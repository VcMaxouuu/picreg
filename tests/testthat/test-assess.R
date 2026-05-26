# assess.pic: family-appropriate predictive metrics + optional support
# recovery diagnostics.

test_that("assess returns a 2-column pic.assess data.frame", {
  set.seed(301)
  n <- 60; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X[, 1L] + rnorm(n))

  fit <- pic(X, y, lambda_n_simu = 200L)
  out <- assess(fit, X, y)

  expect_s3_class(out, "pic.assess")
  expect_s3_class(out, "data.frame")
  expect_equal(names(out), c("metric", "value"))
  expect_type(out$metric, "character")
  expect_type(out$value, "double")
})

test_that("Gaussian assess reports MSE, MAE, R2", {
  set.seed(302)
  n <- 80; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X[, 1L] + rnorm(n))

  fit <- pic(X, y, lambda_n_simu = 200L)
  out <- assess(fit, X, y)
  expect_equal(out$metric, c("MSE", "MAE", "R2"))
})

test_that("Binomial assess reports accuracy, AUC, deviance", {
  set.seed(303)
  n <- 200; p <- 8
  X <- matrix(rnorm(n * p), n, p)
  b <- numeric(p); b[1] <- 2
  y <- rbinom(n, 1L, plogis(X %*% b))

  fit <- pic(X, y, family = "binomial", lambda_n_simu = 200L)
  out <- assess(fit, X, y)
  expect_equal(out$metric, c("accuracy", "AUC", "deviance"))

  acc <- out$value[out$metric == "accuracy"]
  expect_gt(acc, 0.5)
  expect_lte(acc, 1)
})

test_that("Cox assess reports c_index and partial_log_likelihood", {
  set.seed(304)
  data(CoxExample)
  fit <- pic(CoxExample$X, CoxExample$y,
             family = "cox", lambda_n_simu = 200L)
  out <- assess(fit, CoxExample$X, CoxExample$y)
  expect_equal(out$metric, c("c_index", "partial_log_likelihood"))
  c_idx <- out$value[out$metric == "c_index"]
  expect_gte(c_idx, 0)
  expect_lte(c_idx, 1)
})

test_that("true_features appends support-recovery metrics", {
  set.seed(305)
  data(QuickStartExample)
  fit <- pic(QuickStartExample$X, QuickStartExample$y, lambda_n_simu = 200L)

  out <- assess(fit, QuickStartExample$X, QuickStartExample$y,
                true_features = paste0("gene_", 1:5))

  for (m in c("exact_recovery", "tpr", "fdr", "f1")) {
    expect_true(m %in% out$metric)
  }
  # Indices should also work.
  preproc <- attr(fit, "preproc")
  true_idx <- match(paste0("gene_", 1:5), preproc$feature_names)
  out2 <- assess(fit, QuickStartExample$X, QuickStartExample$y,
                 true_features = true_idx)
  for (m in c("exact_recovery", "tpr", "fdr", "f1"))
    expect_equal(out$value[out$metric == m],
                 out2$value[out2$metric == m])
})

test_that("assess errors on bad Cox newy", {
  set.seed(306)
  data(CoxExample)
  fit <- pic(CoxExample$X, CoxExample$y,
             family = "cox", lambda_n_simu = 200L)
  expect_error(
    assess(fit, CoxExample$X, CoxExample$y[, 1L]),
    "2-column matrix"
  )
})
