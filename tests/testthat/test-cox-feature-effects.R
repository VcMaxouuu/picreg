# feature_effects_on_survival: only needs the fit (preproc cache).
# Also covers the predict(type='survival') output shape and the
# composition with plot_survival_curves().

test_that("feature_effects_on_survival composes with predict(type='survival')", {
  set.seed(401)
  data(CoxExample)
  fit <- pic(CoxExample$X, CoxExample$y,
             family = "cox", lambda_n_simu = 200L)

  # Default values: cached grid (quantiles or unique vals).
  fx <- feature_effects_on_survival(fit, idx = "gene_1")
  expect_named(fx, c("time", "survival"))
  expect_true(is.numeric(fx$time))
  expect_true(is.matrix(fx$survival))
  expect_equal(nrow(fx$survival), length(fx$time))
  expect_true(all(grepl("^gene_1 = ", colnames(fx$survival))))
  expect_true(all(fx$survival >= 0 & fx$survival <= 1))
})

test_that("feature_effects_on_survival accepts user values and integer idx", {
  set.seed(402)
  data(CoxExample)
  fit <- pic(CoxExample$X, CoxExample$y,
             family = "cox", lambda_n_simu = 200L)

  vals <- c(-2, 0, 2)
  fx <- feature_effects_on_survival(fit, idx = "gene_2", values = vals)
  expect_equal(ncol(fx$survival), length(vals))

  # Same call via integer index — should give identical result.
  idx_int <- match("gene_2", colnames(CoxExample$X))
  fx2 <- feature_effects_on_survival(fit, idx = idx_int, values = vals)
  expect_equal(fx$survival, fx2$survival)
})

test_that("feature_effects_on_survival rejects non-selected features", {
  set.seed(403)
  data(CoxExample)
  fit <- pic(CoxExample$X, CoxExample$y,
             family = "cox", lambda_n_simu = 200L)
  expect_error(
    feature_effects_on_survival(fit, idx = "noise_1"),
    "selected support"
  )
})

test_that("predict(type='survival') returns a list of (time, survival)", {
  set.seed(404)
  data(CoxExample)
  fit <- pic(CoxExample$X, CoxExample$y,
             family = "cox", lambda_n_simu = 200L)

  sf <- predict(fit, newx = CoxExample$X[1:3, ], type = "survival")
  expect_named(sf, c("time", "survival"))
  expect_equal(ncol(sf$survival), 3L)
  expect_equal(nrow(sf$survival), length(sf$time))
})

test_that("plot_survival_curves runs without error", {
  set.seed(405)
  data(CoxExample)
  fit <- pic(CoxExample$X, CoxExample$y,
             family = "cox", lambda_n_simu = 200L)
  sf <- predict(fit, newx = CoxExample$X[1:4, ], type = "survival")

  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp); on.exit(unlink(tmp), add = TRUE)
  expect_silent(plot_survival_curves(sf))
  dev.off()
})

test_that("feature_effects_on_survival works without standardization", {
  data(CoxExample)
  fit_std <- pic(CoxExample$X, CoxExample$y, family = "cox", lambda = 0.05)
  fit_raw <- suppressWarnings(
    pic(CoxExample$X, CoxExample$y, family = "cox", lambda = 0.05,
        standardize = FALSE)
  )
  feat <- fit_raw$selected[1L]
  fx <- feature_effects_on_survival(fit_raw, idx = feat, values = c(-1, 1))
  expect_equal(ncol(fx$survival), 2L)

  # The single batched prediction matches one prediction per value.
  X_mean <- colMeans(CoxExample$X)
  row <- matrix(X_mean, 1L, dimnames = list(NULL, colnames(CoxExample$X)))
  row[1L, feat] <- 1
  expect_equal(unname(fx$survival[, 2L]),
               as.numeric(predict(fit_raw, row, type = "survival")$survival))
})
