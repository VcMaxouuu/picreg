test_that("print.pic.phase_transition lists every curve", {
  set.seed(501)
  pt <- phase_transition(n = c(30, 40), p = c(15, 15), type = "gaussian",
                         s_max = 1L, m = 2L, penalty = c("lasso", "scad"),
                         lambda_n_simu = 100L)
  out <- capture.output(print(pt))
  expect_true(any(grepl("curves    : 4", out)))
  expect_true(any(grepl("n = 40, p = 15, penalty = scad", out)))
  expect_equal(sum(pt$results$n_failed), 0L)
})

test_that("failed fits are excluded from the averages and reported", {
  local_mocked_bindings(pic = function(...) stop("boom"))
  expect_warning(
    pt <- phase_transition(n = 20, p = 10, type = "gaussian", s_max = 1L,
                           m = 2L, lambda_n_simu = 50L),
    "4 of 4 fits failed.*boom"
  )
  expect_equal(pt$results$n_failed, c(2L, 2L))
  expect_true(all(is.nan(pt$results$exact_recovery)))
})

test_that("parallel phase_transition restores the caller's future plan", {
  skip_on_cran()
  skip_if_not_installed("future.apply")
  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)
  set.seed(502)
  pt <- phase_transition(n = 20, p = 10, type = "gaussian", s_max = 1L,
                         m = 2L, lambda_n_simu = 50L,
                         parallel = TRUE, workers = 1L)
  expect_s3_class(pt, "pic.phase_transition")
  expect_true(inherits(future::plan(), "sequential"))
})
