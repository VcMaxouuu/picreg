#' S3 methods for fitted pic objects.
#'
#' The fitted object has class `c("pic.<family>", "pic")` —
#'  methods dispatch on the second class for shared behavior
#'  and on the first class for family-specific overrides.
#'
#' @name pic_methods
NULL

#' Linear predictor / response prediction for a pic fit.
#'
#' @param object A fitted `pic` object.
#' @param newx Matrix of new values at which predictions are to be made.
#' @param type `"link"` (linear predictor) or `"response"` (default; family `g`-link applied).
#' @param ... Unused; present for S3 method consistency.
#' @return A numeric vector.
#' @export
predict.pic <- function(object, newx,
                        type = c("response", "link", "class"),
                        ...) {
  type <- match.arg(type)
  preproc <- attr(object, "preproc")
  px <- check_X(newx,
                standardize_X = isTRUE(preproc$standardize),
                X_mean = preproc$X_mean, X_std = preproc$X_std)
  X <- px$X

  eta <- as.numeric(X %*% object$beta)
  if (!is.null(object$intercept))
    eta <- eta + object$intercept

  if (type == "link") return(eta)
  if (type == "class") {
    if (!inherits(object, "pic.binomial"))
      stop("type = 'class' is only defined for binomial fits.")
    p <- object$family$g$fn(eta)
    return(as.integer(p >= 0.5))
  }
  # type == "response"
  object$family$g$fn(eta)
}

#' Coefficients of a fitted pic model.
#'
#' Returns the fitted coefficients as a one-column **sparse matrix** of class
#' `"dgCMatrix"` (from the \pkg{Matrix} package), mirroring the output of
#' `coef()` for \pkg{glmnet}. The first row is the intercept, labeled
#' `"(Intercept)"` (value `0` when no intercept was fitted); the remaining
#' rows are the predictors. Row names are taken from the column names of `X`
#' when available (matrix `colnames(X)` or data-frame column names) and
#' otherwise default to `V1, ..., Vp`. Zero coefficients are stored implicitly
#' and printed as `"."`, which keeps the display compact in high dimensions.
#'
#' Internally the model is fitted on a standardized design matrix, so the
#' raw coefficients live on the standardized scale. By default this method
#' rescales them back to the **original scale** of `X` — the values to
#' plug into the un-standardized design for prediction — via
#' \deqn{beta\_orig = beta / s \quad intercept\_orig = intercept - sum(m * beta\_orig)}
#' where `m` and `s` are the column mean and standard deviation. Pass
#' `standardized = TRUE` to skip the rescaling and return the raw fit
#' values (identical to `fit$beta`).
#'
#' @param object A fitted `pic` object.
#' @param standardized Logical; if `TRUE`, return the coefficients on the
#'   standardized scale used internally during fitting. Default `FALSE`
#'   (return coefficients on the original scale of `X`).
#' @param ... Unused; present for S3 method consistency.
#' @return A sparse `(p + 1)` by `1` matrix of class `"dgCMatrix"`, with row
#'   names `c("(Intercept)", <variables>)` and column name `"coefficient"`.
#'   Use `as.numeric()` to obtain a plain numeric vector, or
#'   `which(coef(fit) != 0)` to list the selected entries.
#' @importFrom Matrix Matrix
#' @export
coef.pic <- function(object, standardized = FALSE, ...) {
  preproc <- attr(object, "preproc")
  nm <- if (!is.null(preproc$feature_names))
          preproc$feature_names
        else paste0("V", seq_along(object$beta))

  beta <- object$beta
  has_intercept <- !is.null(object$intercept)
  intercept <- if (has_intercept) object$intercept else 0.0

  unstd <- !standardized && isTRUE(preproc$standardize) &&
    !is.null(preproc$X_std) && !is.null(preproc$X_mean)
  if (unstd) {
    beta <- beta / preproc$X_std
    if (has_intercept) {
      intercept <- intercept - sum(preproc$X_mean * beta)
    }
  }

  out <- Matrix::Matrix(c(intercept, as.numeric(beta)),
                        ncol = 1L, sparse = TRUE)
  dimnames(out) <- list(c("(Intercept)", nm), "coefficient")
  out
}

#' @export
print.pic <- function(x, ...) {
  cat("pic fit (", class(x)[1L], ")\n", sep = "")
  cat("  family   : ", x$family$name, "\n", sep = "")
  cat("  penalty  : ", x$penalty$name, "\n", sep = "")
  cat("  lambda   : ", format(x$lambda, digits = 6), "\n", sep = "")
  n_sel <- sum(x$beta != 0)
  cat("  selected : ", n_sel, " / ", length(x$beta), "\n", sep = "")
  if (!is.null(x$intercept))
    cat("  intercept: ", format(x$intercept, digits = 6), "\n", sep = "")
  invisible(x)
}

#' Summarize a fitted pic model.
#'
#' Returns a structured summary of the fit: the family, penalty, selected
#' \eqn{\hat\lambda}, the problem dimensions, the number of selected
#' variables, the intercept, and a table of the **non-zero** coefficients
#' ordered by decreasing absolute value. Coefficients are returned on the
#' original scale of `X` by default (pass `standardized = TRUE` for the
#' internal standardized scale, matching `fit$beta`).
#'
#' @param object A fitted `pic` object.
#' @param standardized Logical; if `TRUE`, summarize the standardized
#'   coefficients used internally. Default `FALSE` (original scale of `X`).
#' @param ... Unused; present for S3 method consistency.
#' @return An object of class `"summary.pic"`: a list with elements
#'   `family`, `penalty`, `lambda`, `n_samples`, `n_features`, `df`,
#'   `intercept`, `standardized`, and `coefficients` (a two-column
#'   `data.frame` of the non-zero coefficients).
#' @examples
#' data(QuickStartExample)
#' fit <- pic(QuickStartExample$X, QuickStartExample$y)
#' summary(fit)
#' @export
summary.pic <- function(object, standardized = FALSE, ...) {
  co <- coef.pic(object, standardized = standardized)
  cv <- as.numeric(co)
  names(cv) <- rownames(co)

  intercept_val <- if (!is.null(object$intercept)) cv[["(Intercept)"]] else NA_real_
  feat <- cv[names(cv) != "(Intercept)"]
  nz   <- feat[feat != 0]
  nz   <- nz[order(abs(nz), decreasing = TRUE)]

  tab <- data.frame(
    variable    = names(nz),
    coefficient = unname(nz),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  out <- list(
    family       = object$family$name,
    penalty      = object$penalty$name,
    lambda       = object$lambda,
    n_samples    = object$n_samples,
    n_features   = object$n_features,
    df           = object$df,
    intercept    = intercept_val,
    standardized = standardized,
    coefficients = tab
  )
  class(out) <- "summary.pic"
  out
}

#' Print a `summary.pic` object.
#'
#' @param x A `summary.pic` object.
#' @param digits Number of significant digits for the printed values.
#' @param ... Unused; present for S3 method consistency.
#' @return Invisibly returns `x`.
#' @export
print.summary.pic <- function(x, digits = 4L, ...) {
  cat("pic fit summary\n")
  cat("  family    : ", x$family, "\n", sep = "")
  cat("  penalty   : ", x$penalty, "\n", sep = "")
  cat("  lambda    : ", format(x$lambda, digits = digits), "\n", sep = "")
  cat("  dimensions: n = ", x$n_samples, ", p = ", x$n_features, "\n", sep = "")
  cat("  selected  : ", x$df, " / ", x$n_features, "\n", sep = "")
  if (!is.na(x$intercept))
    cat("  intercept : ", format(x$intercept, digits = digits), "\n", sep = "")

  scale_lbl <- if (isTRUE(x$standardized)) "standardized scale" else "original scale"
  if (nrow(x$coefficients) > 0L) {
    cat("\n  Non-zero coefficients (", scale_lbl, "):\n", sep = "")
    tab <- x$coefficients
    tab$coefficient <- round(tab$coefficient, digits)
    print.data.frame(tab, row.names = FALSE)
  } else {
    cat("\n  No variables selected.\n")
  }
  invisible(x)
}

#' Summary of the PDB lambda selector.
#'
#' Prints a formatted summary of the Pivotal Detection Boundary selector
#' used by [pic()] to choose \eqn{\hat\lambda}: method, nominal level,
#' Monte Carlo size, selected \eqn{\hat\lambda}, and a compact view of
#' the null distribution when Monte Carlo was run. For models fitted
#' with `method = "analytical"` or with a user-supplied `lambda`, only
#' the selector metadata is shown.
#'
#' @param model A fitted `pic` object.
#' @param digits Number of significant digits used in the distribution table.
#'
#' @return Invisibly returns `NULL`. Called for its side effect (printing).
#'
#' @examples
#' data(QuickStartExample)
#' fit <- pic(QuickStartExample$X, QuickStartExample$y,
#'            family = "gaussian", penalty = "lasso")
#' pdb_summary(fit)
#' @export
pdb_summary <- function(model, digits = 4L) {
  if (!inherits(model, "pic"))
    stop("`model` must be a fitted pic object.")
  
  pdb <- model$lambda_pdb
  
  cat("PDB lambda selector\n")
  cat(strrep("-", 19L), "\n", sep = "")
  
  if (is.null(pdb)) {
    cat("  method       : user_supplied\n")
    cat("  lambda_hat   : ", format(model$lambda, digits = digits),
        "\n", sep = "")
    return(invisible(NULL))
  }
  
  cat("  method       : ", pdb$method, "\n", sep = "")
  cat("  alpha        : ", format(pdb$alpha, digits = 3L), "\n", sep = "")
  if (!is.na(pdb$n_simu))
    cat("  n_simu       : ", format(pdb$n_simu, big.mark = ","),
        "\n", sep = "")
  cat("  lambda_hat   : ", format(pdb$value, digits = digits),
      "\n", sep = "")
  
  s <- pdb$statistics
  if (!is.null(s)) {
    qs <- stats::quantile(s, c(0.05, 0.25, 0.50, 0.75, 0.95), names = FALSE)
    q  <- c(min = min(s), q05 = qs[1L], q25 = qs[2L], median = qs[3L],
            q75 = qs[4L], q95 = qs[5L], max = max(s))
    
    cat("\n  Null distribution:\n")
    print(round(q, digits), print.gap = 2L)
    cat(sprintf("\n  mean = %.*f      sd = %.*f\n",
                digits, mean(s), digits, stats::sd(s)))
  }
  
  invisible(NULL)
}