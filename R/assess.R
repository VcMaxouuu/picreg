#' Assess performance of a `pic` fit.
#'
#' Given a test set, reports a small set of family-appropriate
#' predictive metrics. Optionally appends support-recovery diagnostics
#' when the true active set is known.
#'
#' The metrics depend on the family:
#' \describe{
#'   \item{`"gaussian"`}{MSE, MAE, R-squared.}
#'   \item{`"binomial"`}{accuracy, AUC, binomial deviance.}
#'   \item{`"poisson"`}{MSE, MAE, Poisson deviance.}
#'   \item{`"exponential"`}{MSE, MAE, Exponential deviance.}
#'   \item{`"gumbel"`}{MSE, MAE, deviance computed from the per-sample
#'         negative log-likelihood with a moment estimate of the scale
#'         parameter \eqn{\hat\sigma = \mathrm{sd}(y - \hat\eta)\sqrt{6}/\pi}.}
#'   \item{`"cox"`}{Harrell's C-index and the Breslow partial
#'         log-likelihood (negative, normalized by `n`).}
#' }
#'
#' When `true_features` is non-`NULL`, four support-recovery metrics
#' are appended (independent of `newx` / `newy`):
#' `exact_recovery`, `tpr` (true-positive rate / sensitivity), `fdr`
#' (false-discovery rate) and `f1` (harmonic mean of precision and
#' recall). Names are accepted in addition to integer positions when
#' the fit carries column names.
#'
#' @param object A fitted `pic` object.
#' @param newx Numeric design matrix at which predictions are
#'   evaluated, with the same columns as the training data.
#' @param newy Response on the new observations. Numeric vector for
#'   all families except Cox; for Cox, a two-column matrix
#'   `(time, event)` matching the training response format.
#' @param true_features Optional integer or character vector listing
#'   the indices (or names) of the true active variables. When
#'   supplied, support-recovery metrics are appended.
#'
#' @return A two-column `data.frame` with columns `metric` (character)
#'   and `value` (numeric).
#'
#' @param ... Unused; present for S3 method consistency.
#'
#' @examples
#' data(QuickStartExample)
#' X <- QuickStartExample$X
#' y <- QuickStartExample$y
#' fit <- pic(X, y, family = "gaussian", penalty = "scad")
#' assess(fit, X, y, true_features = paste0("gene_", 1:5))
#' @export
assess <- function(object, ...) UseMethod("assess")

#' @rdname assess
#' @export
assess.pic <- function(object, newx, newy, true_features = NULL, ...) {
  family <- object$family$name

  if (family == "cox") {
    if (!is.matrix(newy) || ncol(newy) != 2L)
      stop("For Cox fits, `newy` must be a 2-column matrix (time, event).")
    metrics <- .assess_cox(object, newx, newy)
  } else {
    newy <- as.numeric(newy)
    eta  <- predict(object, newx = newx, type = "link")
    mu   <- object$family$g$fn(eta)
    metrics <- switch(family,
      gaussian    = .assess_gaussian(newy, mu),
      binomial    = .assess_binomial(newy, mu),
      poisson     = .assess_poisson(newy, mu),
      exponential = .assess_exponential(newy, mu),
      gumbel      = .assess_gumbel(newy, eta)
    )
  }

  if (!is.null(true_features)) {
    metrics <- rbind(metrics, .support_recovery(object, true_features))
  }

  rownames(metrics) <- NULL
  class(metrics) <- c("pic.assess", "data.frame")
  metrics
}

#' @export
print.pic.assess <- function(x, ...) {
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}


# ---- family-specific helpers -----------------------------------------------
# Each returns a 2-column data.frame (metric, value).

.assess_gaussian <- function(y, mu) {
  resid <- y - mu
  ss_tot <- sum((y - mean(y))^2)
  ss_res <- sum(resid^2)
  r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
  data.frame(
    metric = c("MSE", "MAE", "R2"),
    value  = c(mean(resid^2), mean(abs(resid)), r2)
  )
}

.assess_binomial <- function(y, p) {
  eps <- 1e-15
  p_clipped <- pmin(pmax(p, eps), 1 - eps)
  y_hat <- as.integer(p >= 0.5)
  data.frame(
    metric = c("accuracy", "AUC", "deviance"),
    value  = c(
      mean(y_hat == y),
      .auc(y, p),
      -2 * mean(y * log(p_clipped) + (1 - y) * log(1 - p_clipped))
    )
  )
}

.assess_poisson <- function(y, mu) {
  mu <- pmax(mu, 1e-15)
  ind <- y > 0
  dev_terms <- ifelse(ind, y * log(y / mu), 0) - (y - mu)
  data.frame(
    metric = c("MSE", "MAE", "deviance"),
    value  = c(mean((y - mu)^2), mean(abs(y - mu)), 2 * mean(dev_terms))
  )
}

.assess_exponential <- function(y, mu) {
  mu <- pmax(mu, 1e-15)
  y_safe <- pmax(y, 1e-15)
  dev_terms <- y_safe / mu - log(y_safe / mu) - 1
  data.frame(
    metric = c("MSE", "MAE", "deviance"),
    value  = c(mean((y - mu)^2), mean(abs(y - mu)), 2 * mean(dev_terms))
  )
}

.assess_gumbel <- function(y, eta) {
  resid <- y - eta
  # Moment estimator: Var(Gumbel(0, sigma)) = (pi^2 / 6) sigma^2.
  sigma_hat <- max(stats::sd(resid) * sqrt(6) / pi, 1e-15)
  z <- resid / sigma_hat
  nll <- log(sigma_hat) + z + exp(-z)
  data.frame(
    metric = c("MSE", "MAE", "deviance"),
    value  = c(mean(resid^2), mean(abs(resid)), 2 * mean(nll))
  )
}

.assess_cox <- function(object, newx, newy) {
  risk <- predict(object, newx = newx, type = "link")
  times  <- newy[, 1L]
  events <- newy[, 2L]
  data.frame(
    metric = c("c_index", "partial_log_likelihood"),
    value  = c(
      concordance_index(times, events, risk),
      cox_partial_log_likelihood(times, events, risk)
    )
  )
}


# ---- support recovery ------------------------------------------------------

.support_recovery <- function(object, true_features) {
  preproc <- attr(object, "preproc")
  nm_fit  <- preproc$feature_names

  # Resolve true_features to integer column positions of the training X.
  if (is.character(true_features)) {
    if (is.null(nm_fit))
      stop("Character `true_features` requires the fit to have feature names.")
    pos <- match(true_features, nm_fit)
    if (any(is.na(pos)))
      stop("Some entries of `true_features` were not found in the fitted model: ",
           paste(true_features[is.na(pos)], collapse = ", "))
    true_idx <- sort(unique(as.integer(pos)))
  } else {
    true_idx <- sort(unique(as.integer(true_features)))
  }

  selected_idx <- sort(which(object$beta != 0))

  inter <- length(intersect(selected_idx, true_idx))
  precision <- if (length(selected_idx) == 0L) NA_real_
               else inter / length(selected_idx)
  recall    <- if (length(true_idx) == 0L) NA_real_
               else inter / length(true_idx)
  f1 <- if (is.na(precision) || is.na(recall) ||
            (precision + recall) == 0) 0.0
        else 2 * precision * recall / (precision + recall)
  fdr <- if (length(selected_idx) == 0L) 0.0
         else length(setdiff(selected_idx, true_idx)) / length(selected_idx)
  exact <- as.integer(identical(selected_idx, true_idx))

  data.frame(
    metric = c("exact_recovery", "tpr", "fdr", "f1"),
    value  = c(exact, recall, fdr, f1)
  )
}


# ---- internal: AUC --------------------------------------------------------

.auc <- function(y, prob) {
  y <- as.integer(y == 1L)
  if (length(unique(y)) < 2L) return(NA_real_)
  ranks <- rank(prob, ties.method = "average")
  n1 <- sum(y == 1L)
  n0 <- length(y) - n1
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
