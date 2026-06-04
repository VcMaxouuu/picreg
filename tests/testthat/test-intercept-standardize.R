test_that("pic() allows intercept=FALSE with standardize=TRUE (assumes centered data)", {
  set.seed(1)
  n <- 50; p <- 10
  X <- scale(matrix(rnorm(n * p), n, p), center = TRUE, scale = FALSE)
  y <- as.numeric(X[, 1] + rnorm(n)); y <- y - mean(y)

  # No longer an error: the user is responsible for centring when intercept = FALSE.
  expect_error(
    pic(X, y, family = "gaussian", intercept = FALSE, standardize = TRUE,
        lambda_n_simu = 200L),
    NA
  )
})

test_that("pic() with intercept=TRUE, standardize=TRUE, lambda=0 matches lm()", {
  set.seed(11)
  n <- 80; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X %*% rnorm(p) + rnorm(n))

  fit <- pic(X, y, family = "gaussian", intercept = TRUE, standardize = TRUE,
             lambda = 0, tol = 1e-10, maxit = 50000L)
  b_pic <- as.numeric(coef(fit))

  df <- data.frame(y = y, X)
  b_lm <- unname(coef(lm(y ~ ., data = df)))

  expect_equal(b_pic, b_lm, tolerance = 1e-4)
})

test_that("pic() with intercept=FALSE, standardize=FALSE, lambda=0 matches lm(y ~ . - 1)", {
  set.seed(12)
  n <- 80; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X %*% rnorm(p) + rnorm(n))

  fit <- suppressWarnings(
    pic(X, y, family = "gaussian", intercept = FALSE, standardize = FALSE,
        lambda = 0, tol = 1e-10, maxit = 50000L)
  )
  coefs <- coef(fit)
  b_pic <- setNames(as.numeric(coefs), rownames(coefs))
  b_pic <- b_pic[names(b_pic) != "(Intercept)"]

  df <- data.frame(y = y, X)
  b_lm <- coef(lm(y ~ . - 1, data = df))

  expect_equal(unname(b_pic), unname(b_lm), tolerance = 1e-4)
})

test_that("Cox ignores intercept silently regardless of standardize", {
  set.seed(13)
  n <- 60; p <- 8
  X <- matrix(rnorm(n * p), n, p)
  time  <- rexp(n, rate = 0.5)
  event <- rbinom(n, 1, 0.7)
  y <- cbind(time = time, event = event)

  expect_error(
    pic(X, y, family = "cox", intercept = TRUE, standardize = TRUE,
        lambda_n_simu = 200L),
    NA
  )
  expect_error(
    pic(X, y, family = "cox", intercept = FALSE, standardize = TRUE,
        lambda_n_simu = 200L),
    NA
  )
})
