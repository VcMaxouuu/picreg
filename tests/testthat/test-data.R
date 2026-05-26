# Sanity checks on the bundled datasets used by the vignette.

test_that("QuickStartExample loads with expected shape and labels", {
  data(QuickStartExample)
  expect_named(QuickStartExample, c("X", "y"))
  expect_equal(dim(QuickStartExample$X), c(100L, 30L))
  expect_length(QuickStartExample$y, 100L)
  nm <- colnames(QuickStartExample$X)
  expect_equal(sum(grepl("^gene_", nm)), 5L)
  expect_equal(sum(grepl("^noise_", nm)), 25L)
})

test_that("CoxExample loads with expected shape and (time, event) format", {
  data(CoxExample)
  expect_named(CoxExample, c("X", "y"))
  expect_equal(ncol(CoxExample$X), 50L)
  expect_equal(nrow(CoxExample$X), 250L)
  expect_equal(ncol(CoxExample$y), 2L)
  expect_equal(colnames(CoxExample$y), c("time", "event"))
  expect_true(all(CoxExample$y[, "event"] %in% c(0L, 1L)))
})
