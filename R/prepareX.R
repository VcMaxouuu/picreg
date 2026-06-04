#' Validate and optionally standardize the design matrix.
#'
#' Standardizes columns to zero mean and unit variance (population sd; `n` divisor)
#' when `standardize_X` is `TRUE`. When `X_mean` and `X_std` are supplied they
#' are reapplied without recomputing — used at prediction time on new data.
#' Column names of `X` (or column names of the source data frame) are captured
#' and returned as `feature_names`; they are propagated through `pic()` to
#' annotate `fit$selected`, `coef()`, and the coefficient plot.
#'
#' @param X A numeric matrix or coercible (data.frame).
#' @param standardize_X Logical; standardize columns of `X`.
#' @param X_mean Optional pre-computed column means.
#' @param X_std Optional pre-computed column standard deviations.
#' @return A list with components `X`, `X_mean`, `X_std`, `feature_names`
#'   (the column names of `X` if any, otherwise `NULL`).
#' @keywords internal
check_X <- function(X,
                    standardize_X = TRUE,
                    X_mean = NULL,
                    X_std = NULL) {
  # A bare numeric vector is interpreted as one observation, not one
  # feature.
  if (is.numeric(X) && is.null(dim(X))) {
    X <- matrix(X, nrow = 1L, dimnames = list(NULL, names(X)))
  }
  X <- as.matrix(X)
  feature_names <- colnames(X)

  if (!is.numeric(X)) {
    stop("`X` must contain only numeric values.", call. = FALSE)
  }
  if (nrow(X) == 0L) {
    stop("`X` must have at least one row.", call. = FALSE)
  }
  if (ncol(X) < 2L) {
    stop("`X` must have at least 2 columns.", call. = FALSE)
  }

  storage.mode(X) <- "double"

  if (any(!is.finite(X))) {
    stop("`X` contains NA, NaN, or Inf values.", call. = FALSE)
  }

  # Fit-time only: pre-compute a small set of representative values per
  # column on the original scale, used downstream by
  # `feature_effects_on_survival()` so the caller does not have to keep
  # the training X around. Heuristic mirrors the runtime default of
  # that function: unique values when there are at most five (handy
  # for ordinal / categorical covariates), otherwise four equispaced
  # empirical quantiles. Suppressed at predict time, detected via
  # supplied X_mean / X_std.
  feature_values <- NULL
  if (is.null(X_mean) && is.null(X_std)) {
    feature_values <- lapply(seq_len(ncol(X)), function(j) {
      col <- X[, j]
      uv  <- unique(col)
      if (length(uv) <= 5L)
        sort(uv)
      else
        as.numeric(stats::quantile(col,
                                   seq(0, 1, length.out = 4L),
                                   names = FALSE))
    })
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
       X_std = X_std,
       feature_names = feature_names,
       feature_values = feature_values)
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
      stop("Survival `y` must be a 2-column matrix (time, event).", call. = FALSE)
    if (any(!is.finite(y)))
      stop("`y` contains NA, NaN, or Inf values.", call. = FALSE)
    if (!is.null(n_samples) && nrow(y) != n_samples)
      stop("Length of `y` does not match the number of rows in `X`.", call. = FALSE)
    if (any(y[, 1L] < 0))
      stop("Survival times must be non-negative.", call. = FALSE)
    if (!all(unique(y[, 2L]) %in% c(0, 1)))
      stop("Survival event indicators must be 0 or 1.", call. = FALSE)
    storage.mode(y) <- "double"
    return(y)
  }
  
  if (is.matrix(y) && (ncol(y) == 1L || nrow(y) == 1L))
    y <- as.numeric(y)
  if (!is.numeric(y) || !is.null(dim(y)))
    stop("`y` must be a 1-D numeric vector.", call. = FALSE)
  if (any(!is.finite(y)))
    stop("`y` contains NA, NaN, or Inf values.", call. = FALSE)
  if (!is.null(n_samples) && length(y) != n_samples)
    stop("Length of `y` does not match the number of rows in `X`.", call. = FALSE)
  
  if (y_kind == "positive") {
    if (any(y < 0))
      stop("`y` must be non-negative.", call. = FALSE)
    return(as.numeric(y))
  }
  if (y_kind == "binary") {
    if (!all(unique(y) %in% c(0, 1)))
      stop("Binary `y` must contain only 0 and 1.", call. = FALSE)
    return(as.integer(y))
  }
  if (y_kind == "continuous")
    return(as.numeric(y))
  stop(sprintf("Unknown `y_kind` '%s'.", y_kind), call. = FALSE)
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
    y = y,
    feature_names = px$feature_names,
    feature_values = px$feature_values
  )
}
