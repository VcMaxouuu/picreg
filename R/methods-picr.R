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
#' @param type One of `"link"` (linear predictor),
#'   `"response"` (default; family `g`-link applied),
#'   `"coefficients"` (return fitted coefficients incl. intercept)
#' @param ... Ignored.
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
#' Returns a named numeric vector `(intercept, b1, ..., bp)`. When no
#' intercept was fitted the first entry is `0`.
#'
#' Internally the model is fitted on a standardised design matrix, so the
#' raw coefficients live on the standardised scale. By default this method
#' returns the **original-scale** coefficients — the values to plug into
#' the un-standardised `X` for prediction. Pass `standardized = TRUE` to
#' get the raw fit values instead.
#'
#' Conversion: with mean `m` and standard deviation `s` per column,
#' `beta_orig = beta / s` and
#' `intercept_orig = intercept - sum(m * beta_orig)`.
#'
#' @param object A fitted `pic` object.
#' @param standardized Logical; if `TRUE`, return the coefficients on the
#'   standardised scale used internally during fitting. Default `FALSE`
#'   (return coefficients on the original scale of `X`).
#' @param ... Ignored.
#' @export
coef.pic <- function(object, standardized = FALSE, ...) {
  nm <- paste0("V", seq_along(object$beta))

  beta <- object$beta
  has_intercept <- !is.null(object$intercept)
  intercept <- if (has_intercept) object$intercept else 0.0

  preproc <- attr(object, "preproc")
  unstd <- !standardized && isTRUE(preproc$standardize) &&
    !is.null(preproc$X_std) && !is.null(preproc$X_mean)
  if (unstd) {
    beta <- beta / preproc$X_std
    if (has_intercept) {
      intercept <- intercept - sum(preproc$X_mean * beta)
    }
  }

  out <- c(intercept, as.numeric(beta))
  names(out) <- c("(Intercept)", nm)
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


#' Summary statistics of the PDB null-distribution selector.
#'
#' Returns a named list of descriptive statistics characterising the
#' Monte Carlo null distribution of the gradient-norm pivot used by the
#' Pivotal Detection Boundary selector. When the model was fitted with
#' `"analytical"` or with a user-supplied `lambda`, only the scalar
#' selector value (and metadata when available) is returned.
#'
#' Components, in order:
#' \describe{
#'   \item{method}{The PDB method (`"mc_exact"`, `"mc_gaussian"`,
#'     `"analytical"`, or `"user_supplied"`).}
#'   \item{alpha}{Nominal level used.}
#'   \item{n_simu}{Number of Monte Carlo draws.}
#'   \item{lambda}{Selected (or fixed) \eqn{\hat\lambda}.}
#'   \item{mean,sd,min,q05,q25,median,q75,q95,max}{Summary statistics of
#'     the Monte Carlo null distribution (only present when MC was run).}
#' }
#'
#' @param model A fitted `pic` object.
#' @return A named list of statistics.
#' @examples
#' \dontrun{
#' fit <- pic(X, y, family = "gaussian", penalty = "lasso")
#' s <- pdb_summary(fit)
#' s$lambda
#' s$mean
#' s$q95   # upper 5% of the null
#' }
#' @export
pdb_summary <- function(model) {
  if (!inherits(model, "pic"))
    stop("`model` must be a fitted pic object.")

  pdb <- model$lambda_pdb
  if (is.null(pdb)) {
    return(list(method = "user_supplied", lambda = model$lambda))
  }

  out <- list(
    method = pdb$method,
    alpha  = pdb$alpha,
    n_simu = pdb$n_simu,
    lambda = pdb$value
  )

  s <- pdb$statistics
  if (!is.null(s)) {
    qs <- stats::quantile(s, c(0.05, 0.25, 0.50, 0.75, 0.95), names = FALSE)
    out$mean   <- mean(s)
    out$sd     <- stats::sd(s)
    out$min    <- min(s)
    out$q05    <- qs[1L]
    out$q25    <- qs[2L]
    out$median <- qs[3L]
    out$q75    <- qs[4L]
    out$q95    <- qs[5L]
    out$max    <- max(s)
  }
  out
}
