test_that("the three penalties all fit and produce positive lambdas", {
  set.seed(21)
  n <- 80; p <- 30
  X <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p); beta[1L:3L] <- 3
  y <- as.numeric(X %*% beta + rnorm(n))

  fits <- lapply(c("lasso", "scad", "mcp"), function(pen) {
    set.seed(21)   # identical PDB seed across penalties
    pic(X, y, family = "gaussian", penalty = pen)
  })

  lambdas <- vapply(fits, `[[`, numeric(1L), "lambda")
  expect_true(all(lambdas > 0))
  # PDB only depends on (X, family, alpha, n_simu) — reseeding makes
  # the three calls draw the same Monte Carlo statistics.
  expect_equal(lambdas[1L], lambdas[2L])
  expect_equal(lambdas[1L], lambdas[3L])
})

test_that("get_penalty() rejects unknown penalty names and invalid params", {
  expect_error(get_penalty("ridge"))
  expect_error(get_penalty(NULL))
  expect_error(get_penalty("scad", scad_a = 1.5))
  expect_error(get_penalty("mcp",  mcp_gamma = 0.5))
})
