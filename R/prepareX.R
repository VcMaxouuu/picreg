#' Validate and optionally standardize the design matrix.
#'
#' Standardizes columns to zero mean and unit variance (population sd; `n` divisor)
#' when `standardize_X` is `TRUE`. When `X_mean` and `X_std` are supplied they
#' are reapplied without recomputing — used at prediction time on new data.
#' Column means are computed whenever `X_mean` is not supplied, even without
#' standardization. Column names of `X` (or column names of the source data
#' frame) are captured and returned as `feature_names`; they are propagated
#' through `pic()` to annotate `fit$selected`, `coef()`, and the coefficient
#' plot.
#'
#' @param X A numeric matrix or coercible (data.frame).
#' @param standardize_X Logical; standardize columns of `X`.
#' @param X_mean Optional pre-computed column means.
#' @param X_std Optional pre-computed column standard deviations.
#' @param keep_feature_values Logical; also return a small grid of
#'   representative values per column (used by Cox fits).
#' @param expected_cols Optional required number of columns.
#' @return A list with components `X`, `X_mean`, `X_std`, `feature_names`
#'   (the column names of `X` if any, otherwise `NULL`) and `feature_values`
#'   (`NULL` unless `keep_feature_values = TRUE`).
#' @keywords internal
check_X <- function(X,
                    standardize_X = TRUE,
                    X_mean = NULL,
                    X_std = NULL,
                    keep_feature_values = FALSE,
                    expected_cols = NULL) {
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
  if (!is.null(expected_cols) && ncol(X) != expected_cols) {
    stop(sprintf("`newx` must have %d columns, as the training `X`; got %d.",
                 expected_cols, ncol(X)), call. = FALSE)
  }

  storage.mode(X) <- "double"

  if (any(!is.finite(X))) {
    stop("`X` contains NA, NaN, or Inf values.", call. = FALSE)
  }

  # Cox fits only: pre-compute a small set of representative values per
  # column on the original scale, used downstream by
  # `feature_effects_on_survival()` so the caller does not have to keep
  # the training X around. Heuristic mirrors the runtime default of
  # that function: unique values when there are at most five (handy
  # for ordinal / categorical covariates), otherwise four equispaced
  # empirical quantiles.
  feature_values <- NULL
  if (keep_feature_values) {
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

  # Column means are kept even without standardization: the Cox profile
  # curves of `feature_effects_on_survival()` hold covariates at them.
  if (is.null(X_mean)) {
    X_mean <- colMeans(X)
  }

  if (standardize_X) {
    if (is.null(X_std)) {
      X_centered <- sweep(X, 2L, X_mean, FUN = "-", check.margin = FALSE)
      X_std <- sqrt(colMeans(X_centered^2))
      # Constant columns: rounding in the mean can leave a spurious sd of
      # order eps * |mean|. They are centered to zero and left unscaled.
      X_std[X_std <= 1e-12 * abs(X_mean)] <- 1
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
#' A binary response may also be logical (`FALSE`/`TRUE` mapped to `0`/`1`)
#' or a two-level factor (first level mapped to `0`, second to `1`).
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

  if (y_kind == "binary")
    y <- .as_binary01(y)
  
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

# Map a logical or two-level factor response to 0/1, as glmnet does:
# FALSE/TRUE -> 0/1, first/second factor level -> 0/1. `levels` imposes the
# training levels when mapping a new response. Other inputs are returned
# unchanged for the numeric checks of `check_y()`.
.as_binary01 <- function(y, levels = NULL) {
  if (is.logical(y)) {
    if (anyNA(y))
      stop("`y` contains NA values.", call. = FALSE)
    return(as.integer(y))
  }
  if (is.factor(y)) {
    lv <- if (is.null(levels)) levels(y) else levels
    if (length(lv) != 2L)
      stop("A factor `y` must have exactly two levels for the binomial ",
           "family.", call. = FALSE)
    idx <- match(as.character(y), lv)
    if (anyNA(idx))
      stop("`y` contains NA values or values outside the levels ",
           paste0("'", lv, "'", collapse = ", "), ".", call. = FALSE)
    return(idx - 1L)
  }
  y
}

#' Joint validation/preprocessing of `(X, y)`.
#'
#' For survival data, rows of `X` and `y` are sorted by ascending time.
#' For a factor binary response, its two levels are returned as `y_levels`.
#'
#' @keywords internal
check_Xy <- function(X,
                     y,
                     y_kind,
                     standardize_X = TRUE,
                     X_mean = NULL,
                     X_std = NULL,
                     keep_feature_values = FALSE) {
  px <- check_X(X,
                standardize_X = standardize_X,
                X_mean = X_mean,
                X_std = X_std,
                keep_feature_values = keep_feature_values)
  y_levels <- if (y_kind == "binary" && is.factor(y)) levels(y) else NULL
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
    feature_values = px$feature_values,
    y_levels = y_levels
  )
}
