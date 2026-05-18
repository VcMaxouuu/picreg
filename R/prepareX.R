#' Validate and optionally standardise the design matrix.
#'
#' Standardises columns to zero mean and unit variance (population sd; `n` divisor)
#' when `standardize_X` is `TRUE`. When `X_mean` and `X_std` are supplied they
#' are reapplied without recomputing — used at prediction time on new data.
#'
#' @param X A numeric matrix or coercible (data.frame).
#' @param standardize_X Logical; standardise columns of `X`.
#' @param X_mean Optional pre-computed column means.
#' @param X_std Optional pre-computed column standard deviations.
#' @return A list with components `X`, `X_mean`, `X_std`.
#' @keywords internal
check_X <- function(X,
                    standardize_X = TRUE,
                    X_mean = NULL,
                    X_std = NULL) {
  X <- as.matrix(X)
  
  if (!is.numeric(X)) {
    stop("X must contain only numeric values.")
  }
  if (nrow(X) == 0L) {
    stop("X must have at least one row.")
  }
  if (ncol(X) < 2L) {
    stop("X must have at least 2 columns.")
  }
  
  storage.mode(X) <- "double"
  
  if (any(!is.finite(X))) {
    stop("X contains NA, NaN, or Inf values.")
  }
  
  if (standardize_X) {
    if (is.null(X_mean)) {
      X_mean <- colMeans(X)
    }
    
    if (is.null(X_std)) {
      X_centered <- sweep(X, 2L, X_mean, FUN = "-", check.margin = FALSE)
      X_std <- sqrt(colMeans(X_centered^2))
      X_std[X_std == 0] <- 1
    }
    
    X <- sweep(X, 2L, X_mean, FUN = "-", check.margin = FALSE)
    X <- sweep(X, 2L, X_std, FUN = "/", check.margin = FALSE)
  }
  
  list(X = X,
       X_mean = X_mean,
       X_std = X_std)
}
#' Validate response according to its distributional kind.
#'
#' @param y Response vector or 2-column matrix (survival).
#' @param n_samples Expected number of observations.
#' @param y_kind One of `"continuous"`, `"positive"`, `"binary"`, `"survival"`.
#' @return Validated `y` (possibly coerced to matrix for survival).
#' @keywords internal
check_y <- function(y,
                    n_samples = NULL,
                    y_kind = "continuous") {
  if (y_kind == "survival") {
    if (!is.matrix(y) || ncol(y) != 2L)
      stop("Survival y must be a 2-column matrix (time, event).")
    if (any(!is.finite(y)))
      stop("y contains NA, NaN, or Inf values.")
    if (!is.null(n_samples) && nrow(y) != n_samples)
      stop("Length of y does not match the number of rows in X.")
    if (any(y[, 1L] < 0))
      stop("Survival times must be non-negative.")
    if (!all(unique(y[, 2L]) %in% c(0, 1)))
      stop("Survival event indicators must be 0 or 1.")
    storage.mode(y) <- "double"
    return(y)
  }
  
  if (is.matrix(y) && (ncol(y) == 1L || nrow(y) == 1L))
    y <- as.numeric(y)
  if (!is.numeric(y) || !is.null(dim(y)))
    stop("y must be a 1-D numeric vector.")
  if (any(!is.finite(y)))
    stop("y contains NA, NaN, or Inf values.")
  if (!is.null(n_samples) && length(y) != n_samples)
    stop("Length of y does not match the number of rows in X.")
  
  if (y_kind == "positive") {
    if (any(y < 0))
      stop("y must be non-negative.")
    return(as.numeric(y))
  }
  if (y_kind == "binary") {
    if (!all(unique(y) %in% c(0, 1)))
      stop("Binary y must contain only 0 and 1.")
    return(as.integer(y))
  }
  if (y_kind == "continuous")
    return(as.numeric(y))
  stop(sprintf("Unknown y_kind '%s'.", y_kind))
}

#' Joint validation/preprocessing of `(X, y)`.
#'
#' For survival data, rows of `X` and `y` are sorted by ascending time.
#'
#' @keywords internal
check_Xy <- function(X,
                     y,
                     y_kind,
                     standardize_X = TRUE,
                     X_mean = NULL,
                     X_std = NULL) {
  px <- check_X(X,
                standardize_X = standardize_X,
                X_mean = X_mean,
                X_std = X_std)
  y <- check_y(y, n_samples = nrow(px$X), y_kind = y_kind)
  
  if (y_kind == "survival") {
    ord <- order(y[, 1L])
    px$X <- px$X[ord, , drop = FALSE]
    y <- y[ord, , drop = FALSE]
  }
  
  list(
    X = px$X,
    X_mean = px$X_mean,
    X_std = px$X_std,
    y = y
  )
}
