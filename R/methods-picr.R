#' S3 methods for fitted pic objects.
#'
#' The fitted object has class `c("pic.<family>", "pic")` —
#'  methods dispatch on the second class for shared behaviour
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
#' Returns a two-column data frame with the variable name in the first
#' column (`variable`) and the fitted coefficient in the second
#' (`coefficient`). The first row is always the intercept, labelled
#' `"(Intercept)"`; when no intercept was fitted its value is `0`.
#' Variable names are taken from the column names of `X` if any were
#' supplied (matrix `colnames(X)` or data-frame column names), and
#' otherwise default to `V1, ..., Vp`.
#'
#' Internally the model is fitted on a standardised design matrix, so the
#' raw coefficients live on the standardised scale. By default this method
#' rescales them back to the **original scale** of `X` — the values to
#' plug into the un-standardised design for prediction — via
#' \deqn{beta\_orig = beta / s \quad intercept\_orig = intercept - sum(m * beta\_orig)}
#' where `m` and `s` are the column mean and standard deviation. Pass
#' `standardized = TRUE` to skip the rescaling and return the raw fit
#' values (identical to `fit$beta`).
#'
#' @param object A fitted `pic` object.
#' @param standardized Logical; if `TRUE`, return the coefficients on the
#'   standardised scale used internally during fitting. Default `FALSE`
#'   (return coefficients on the original scale of `X`).
#' @param ... Unused; present for S3 method consistency.
#' @return A `data.frame` with columns `variable` (character) and
#'   `coefficient` (numeric), of length `p + 1` (intercept + features).
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

  out <- data.frame(
    variable    = c("(Intercept)", nm),
    coefficient = c(intercept, as.numeric(beta)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  class(out) <- c("pic.coef", "data.frame")
  out
}

#' Print a `pic.coef` table.
#'
#' Pretty-prints the two-column coefficient table returned by
#' [coef.pic()] without the integer row numbers that `print.data.frame`
#' would otherwise add. The underlying object is still a `data.frame`,
#' so any subsetting / column access works as usual.
#'
#' @param x A `pic.coef` object.
#' @param ... Forwarded to [print.data.frame()].
#' @return Invisibly returns `x`.
#' @export
print.pic.coef <- function(x, ...) {
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
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