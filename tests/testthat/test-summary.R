test_that("summary.pic returns a structured summary with a non-zero coef table", {
  set.seed(7)
  n <- 60; p <- 12
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("g", seq_len(p))
  y <- as.numeric(X[, 1] - X[, 2] + rnorm(n))

  fit <- pic(X, y, family = "gaussian", lambda_n_simu = 200L)
  s   <- summary(fit)

  expect_s3_class(s, "summary.pic")
  expect_equal(s$family, "gaussian")
  expect_equal(s$n_features, p)
  expect_true(is.data.frame(s$coefficients))
  expect_named(s$coefficients, c("variable", "coefficient"))

  # The non-zero table holds exactly df features (intercept excluded).
  expect_equal(nrow(s$coefficients), s$df)
  expect_false("(Intercept)" %in% s$coefficients$variable)

  # Printing works and is labeled.
  expect_output(print(s), "pic fit summary")
})

test_that("summary.pic handles the empty-support case", {
  set.seed(8)
  n <- 40; p <- 8
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)                       # pure noise: likely empty support

  fit <- pic(X, y, family = "gaussian", lambda_n_simu = 200L)
  s   <- summary(fit)

  expect_s3_class(s, "summary.pic")
  expect_true(nrow(s$coefficients) == s$df)
  expect_output(print(s), "pic fit summary")
})
