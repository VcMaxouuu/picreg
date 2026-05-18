test_that("lambda_pdb() returns a scalar lambda and metadata for the three methods", {
  set.seed(31)
  n <- 80; p <- 40
  X <- matrix(rnorm(n * p), n, p)
  X <- sweep(X, 2L, colMeans(X), "-")
  sv <- sqrt(colMeans(X^2)); sv[sv == 0] <- 1
  X <- sweep(X, 2L, sv, "/")

  fam <- get_family("gaussian")

  for (m in c("mc_exact", "mc_gaussian", "analytical")) {
    out <- lambda_pdb(X, fam, n_simu = 200L, alpha = 0.05, method = m)
    expect_s3_class(out, "pic.lambda_pdb")
    expect_true(out$value > 0, info = sprintf("method = %s", m))
    expect_equal(out$method, m)
    expect_equal(out$alpha, 0.05)
  }
})

test_that("lambda_pdb warns when X is not standardised", {
  set.seed(32)
  X <- matrix(rnorm(50 * 10, mean = 5, sd = 2), 50, 10)
  expect_warning(
    lambda_pdb(X, "gaussian", n_simu = 100L, method = "analytical"),
    "standardised"
  )
})

test_that("pdb_summary() returns the lambda and Monte Carlo statistics", {
  set.seed(33)
  n <- 60; p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X[, 1] + rnorm(n))

  fit <- pic(X, y, family = "gaussian", lambda_n_simu = 200L,
             lambda_method = "mc_exact")
  s <- pdb_summary(fit)

  expect_equal(s$method, "mc_exact")
  expect_equal(s$lambda, fit$lambda)
  expect_true("q95" %in% names(s))
})
