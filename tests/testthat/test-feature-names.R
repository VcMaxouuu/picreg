# Feature names captured at preprocessing time are propagated to
# fit$selected, coef(fit), and the preproc attribute used by
# feature_effects_on_survival().

test_that("feature names from a data.frame propagate through the fit", {
  set.seed(201)
  n <- 60; p <- 8
  df <- as.data.frame(matrix(rnorm(n * p), n, p))
  names(df) <- paste0("v_", letters[seq_len(p)])
  y <- as.numeric(df$v_a + rnorm(n))

  fit <- pic(df, y, lambda_n_simu = 200L)

  # fit$selected returns names (not integers).
  expect_type(fit$selected, "character")
  expect_true(all(fit$selected %in% names(df)))

  # coef(fit) uses the same names.
  co <- coef(fit)
  expect_equal(co$variable[-1L], names(df))

  # preproc stores them.
  preproc <- attr(fit, "preproc")
  expect_equal(preproc$feature_names, names(df))
})

test_that("matrices without colnames get V1..Vp labelling", {
  set.seed(202)
  n <- 50; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X[, 1L] + rnorm(n))

  fit <- pic(X, y, lambda_n_simu = 200L)

  expect_type(fit$selected, "integer")
  expect_equal(coef(fit)$variable[-1L], paste0("V", seq_len(p)))
  expect_null(attr(fit, "preproc")$feature_names)
})

test_that("feature_values are cached on Cox fits only", {
  set.seed(203)
  n <- 80; p <- 10

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("g", seq_len(p))
  y_gauss <- as.numeric(X[, 1L] + rnorm(n))
  fit_g <- pic(X, y_gauss, lambda_n_simu = 200L)

  # Gaussian fit must NOT carry the per-column value grids.
  expect_null(attr(fit_g, "preproc")$feature_values)

  # Cox fit must carry them (one entry per training column).
  times  <- rexp(n, rate = exp(X[, 1L]))
  censor <- rexp(n, rate = 1)
  y_cox  <- cbind(time  = pmin(times, censor),
                  event = as.integer(times <= censor))
  fit_c  <- pic(X, y_cox, family = "cox", lambda_n_simu = 200L)

  fv <- attr(fit_c, "preproc")$feature_values
  expect_type(fv, "list")
  expect_length(fv, p)
  expect_true(all(vapply(fv, is.numeric, logical(1L))))
})
