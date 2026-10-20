test_that("check_X rejects non-finite, non-numeric, and too-narrow matrices", {
  expect_error(check_X(matrix("a", 5, 5)), "numeric")
  expect_error(check_X(matrix(NA_real_, 5, 5)), "NA")
  expect_error(check_X(matrix(0, 0, 5)), "row")
  expect_error(check_X(matrix(0, 5, 1)), "columns")
})

test_that("check_y enforces the response kind", {
  expect_error(check_y(c(0.5, 1, 0), y_kind = "binary"), "0 and 1")
  expect_error(check_y(c(-1, 1, 2), y_kind = "positive"), "non-negative")
  expect_error(check_y(c(1, NA, 2), y_kind = "continuous"), "NA")
  expect_error(check_y(matrix(1:4, 2, 2), y_kind = "continuous"), "1-D")
})

test_that("Gumbel rejects a constant response with no location-scale MLE", {
  set.seed(901)
  X <- matrix(rnorm(40), 10, 4)
  expect_error(
    pic(X, rep(2, 10), family = "gumbel", lambda = 0.1),
    "at least two distinct response values"
  )
})

test_that("survival y must be a 2-column matrix with valid event indicators", {
  expect_error(check_y(c(1, 2, 3), y_kind = "survival"), "2-column")
  expect_error(check_y(cbind(c(-1, 2), c(0, 1)), y_kind = "survival"),
               "non-negative")
  expect_error(check_y(cbind(c(1, 2), c(0, 2)), y_kind = "survival"),
               "0 or 1")
})

test_that("pic() validates its scalar arguments", {
  set.seed(902)
  X <- matrix(rnorm(200), 40, 5)
  y <- as.numeric(X[, 1] + rnorm(40))

  expect_error(pic(X, y, lambda = c(0.1, 0.2)), "`lambda` must be NULL")
  expect_error(pic(X, y, lambda = -1), "`lambda` must be NULL")
  expect_error(pic(X, y, lambda = "a"), "`lambda` must be NULL")
  expect_error(pic(X, y, lambda_alpha = 0), "`lambda_alpha`")
  expect_error(pic(X, y, lambda_alpha = 1), "`lambda_alpha`")
  expect_error(pic(X, y, lambda_n_simu = 0), "`lambda_n_simu`")
  expect_error(pic(X, y, lambda_n_simu = 10.5), "`lambda_n_simu`")
  expect_error(pic(X, y, lambda_seed = c(1, 2)), "`lambda_seed`")
  expect_error(pic(X, y, tol = 0), "`tol`")
  expect_error(pic(X, y, path_length = 0), "`path_length`")
  expect_error(pic(X, y, maxit = -1), "`maxit`")
  expect_error(pic(X, y, standardize = NA), "`standardize`")
  expect_error(pic(X, y, intercept = "yes"), "`intercept`")
  expect_error(pic(X, y, penalty = "scad", scad_a = c(3, 4)), "`scad_a`")
  expect_error(pic(X, y, penalty = "mcp", mcp_gamma = NA), "`mcp_gamma`")
})

test_that("lambda_pdb() validates its arguments", {
  X <- scale(matrix(rnorm(200), 40, 5))
  expect_error(lambda_pdb(as.data.frame(X), "gaussian"), "numeric matrix")
  expect_error(lambda_pdb(X, "gaussian", alpha = 1.5), "`alpha`")
  expect_error(lambda_pdb(X, "gaussian", n_simu = 0), "`n_simu`")
  expect_error(lambda_pdb(X, "gaussian", seed = "a"), "`seed`")
})

test_that("predict() reports a column mismatch clearly", {
  set.seed(903)
  X <- matrix(rnorm(300), 30, 10)
  fit <- pic(X, X[, 1] + rnorm(30), lambda_n_simu = 100L)
  expect_error(predict(fit, X[, 1:5]), "must have 10 columns")
})

test_that("a constant column does not trigger the standardization warning", {
  set.seed(904)
  X <- matrix(rnorm(500), 50, 10)
  X[, 3] <- 0.1
  y <- as.numeric(X[, 1] + rnorm(50))
  expect_no_warning(fit <- pic(X, y, lambda_n_simu = 200L))
  expect_equal(fit$beta[3], 0)
  # A raw, non-standardized design still warns.
  expect_warning(lambda_pdb(X, "gaussian", n_simu = 50L), "standardized")
})

test_that("binomial accepts logical and two-level factor responses", {
  data(BinomialExample)
  X <- BinomialExample$X
  y <- BinomialExample$y

  fit_num <- pic(X, y, family = "binomial", lambda = 0.1)
  fit_lgl <- pic(X, y == 1, family = "binomial", lambda = 0.1)
  fit_fac <- pic(X, factor(y, levels = 0:1, labels = c("no", "yes")),
                 family = "binomial", lambda = 0.1)
  expect_equal(fit_lgl$beta, fit_num$beta)
  expect_equal(fit_fac$beta, fit_num$beta)
  expect_equal(attr(fit_fac, "preproc")$y_levels, c("no", "yes"))

  # assess() maps a factor response with the training levels.
  y_fac <- factor(ifelse(y == 1, "yes", "no"), levels = c("no", "yes"))
  expect_equal(assess(fit_fac, X, y_fac), assess(fit_num, X, y))

  expect_error(pic(X, factor(sample(letters[1:3], 300, TRUE)),
                   family = "binomial", lambda = 0.1), "two levels")
  expect_error(pic(X, c(NA, y[-1] == 1), family = "binomial", lambda = 0.1),
               "NA")
})
