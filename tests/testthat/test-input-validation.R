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

test_that("survival y must be a 2-column matrix with valid event indicators", {
  expect_error(check_y(c(1, 2, 3), y_kind = "survival"), "2-column")
  expect_error(check_y(cbind(c(-1, 2), c(0, 1)), y_kind = "survival"),
               "non-negative")
  expect_error(check_y(cbind(c(1, 2), c(0, 2)), y_kind = "survival"),
               "0 or 1")
})
