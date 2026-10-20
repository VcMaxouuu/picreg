#' Survival utilities for the Cox family.
#'
#' Provides Harrell's C-index, the Breslow partial log-likelihood, baseline
#' cumulative hazard / survival, and feature-effect curves.
#'
#' @name pic_survival
NULL

#' Harrell's concordance index.
#'
#' A pair \eqn{(i, j)} is comparable when \eqn{i} has an event and either
#' \eqn{t_i < t_j}, or \eqn{t_i = t_j} and \eqn{j} is censored. It is
#' concordant when \eqn{i} has the higher risk score; tied scores count
#' \eqn{1/2}. Computed in \eqn{O(n \log n)}.
#'
#' @param times Survival times (length `n`).
#' @param events Event indicators (0 or 1; length `n`).
#' @param predictions Risk scores; higher = higher risk = shorter survival.
#' @return Numeric scalar in `[0, 1]` (`0` when no pair is comparable).
#' @export
concordance_index <- function(times, events, predictions) {
  times <- as.numeric(times)
  events <- as.numeric(events)
  predictions <- as.numeric(predictions)
  if (length(events) != length(times) || length(predictions) != length(times))
    stop("`times`, `events`, and `predictions` must have the same length.",
         call. = FALSE)
  concordance_index_cpp(times, events, predictions)
}

#' Breslow partial log-likelihood (negative, normalized by `n`).
#'
#' @inheritParams concordance_index
#' @return A numeric scalar: the negative Breslow partial log-likelihood
#'   for the Cox proportional-hazards model, normalized by the sample
#'   size `n`. Lower values indicate a better fit.
#' @export
cox_partial_log_likelihood <- function(times, events, predictions) {
  ord <- order(times, decreasing = TRUE)
  t_s <- times[ord]; p_s <- predictions[ord]; e_s <- as.numeric(events[ord])
  # Tie groups by decreasing time; exponentials shifted by max(p) for
  # numerical stability.
  m <- max(p_s)
  grp <- rowsum(cbind(p_s * e_s, exp(p_s - m), e_s), t_s, reorder = FALSE)
  log_risk <- log(cumsum(grp[, 2L])) + m
  d <- grp[, 3L]
  ev <- d > 0
  -(sum(grp[, 1L]) - sum(d[ev] * log_risk[ev])) / length(times)
}

#' Breslow baseline cumulative hazard and survival.
#'
#' @inheritParams concordance_index
#' @return A data.frame with columns `time`, `cumulative_hazard`, `survival`.
#' @export
baseline_functions <- function(times, events, predictions) {
  exp_pred <- exp(as.numeric(predictions))
  df <- data.frame(durations = as.numeric(times),
                   exp_pred = exp_pred,
                   events   = as.numeric(events))
  agg <- aggregate(cbind(exp_pred, events) ~ durations, data = df, FUN = sum)
  ord <- order(agg$durations, decreasing = TRUE)
  agg <- agg[ord, , drop = FALSE]
  risk_set_sum <- cumsum(agg$exp_pred)
  d <- agg$events
  hazard_inc <- ifelse(d > 0, d / risk_set_sum, 0)
  H0_desc <- rev(cumsum(rev(hazard_inc)))
  out <- data.frame(time = rev(agg$durations),
                    cumulative_hazard = rev(H0_desc),
                    survival = exp(-rev(H0_desc)))
  out
}

#' Effect of one feature on the Cox survival curve.
#'
#' Builds the family of survival curves obtained by varying a single
#' covariate while holding the others at their column means. The
#' returned object has the same shape as [predict.pic()] with
#' `type = "survival"`, so it composes directly with [plot_survival_curves()].
#'
#' Both the per-column mean row and a default grid of representative
#' values (used when `values = NULL`) are cached on the fit by `pic()`
#' under `attr(fit, "preproc")`, so the training design matrix does
#' not need to be passed back in. The cached default grid uses the
#' unique values of the column when there are at most five of them
#' (handy for ordinal / categorical covariates), and the four
#' equispaced empirical quantiles (`0`, `1/3`, `2/3`, `1`) otherwise.
#'
#' @param object A fitted `pic.cox` object.
#' @param idx Index of the feature to vary. Either an integer column
#'   position in the training `X` or, if `X` carried column names, the
#'   variable name. The feature must lie in the model's selected
#'   support; otherwise the curve would be flat in `v` and the call
#'   is rejected.
#' @param values Optional numeric vector of values to evaluate. When
#'   `NULL` (default), the cached grid described in *Details* is
#'   used.
#' @return A list with components `time` (length `K`) and `survival`
#'   (matrix `K x length(values)`, one column per evaluated value).
#'   Column names of `survival` are formatted as
#'   `"<feature_name> = <value>"` and are picked up automatically by
#'   [plot_survival_curves()] for the legend.
#'
#' @seealso [predict.pic()], [plot_survival_curves()].
#' @export
feature_effects_on_survival <- function(object, idx, values = NULL) {
  if (!inherits(object, "pic.cox"))
    stop("`object` must be a fitted pic.cox model.")

  preproc <- attr(object, "preproc")
  if (is.null(preproc) || is.null(preproc$X_mean) ||
      is.null(preproc$feature_values)) {
    stop("Cached preprocessing summaries not found on this fit. ",
         "Refit the model with the current version of `pic()`.")
  }
  X_mean <- preproc$X_mean
  nm_fit <- preproc$feature_names
  p      <- length(X_mean)

  # Resolve `idx` into an integer column position and a display name.
  if (is.character(idx)) {
    if (length(idx) != 1L)
      stop("`idx` must be a single integer or a single feature name.")
    if (is.null(nm_fit))
      stop("Character `idx` requires the fit to have feature names.")
    pos <- match(idx, nm_fit)
    if (is.na(pos))
      stop(sprintf("Feature '%s' not found in the fitted model.", idx))
    idx_int      <- pos
    feature_name <- idx
  } else {
    idx_int <- as.integer(idx)
    if (length(idx_int) != 1L || is.na(idx_int) ||
        idx_int < 1L || idx_int > p)
      stop("`idx` must be a valid column position in the training X.")
    feature_name <- if (!is.null(nm_fit)) nm_fit[idx_int]
                    else paste0("feature_", idx_int)
  }

  # Verified against the integer support (robust to whether selected is
  # stored as indices or as names).
  selected_int <- which(object$beta != 0)
  if (!(idx_int %in% selected_int)) {
    stop(sprintf(
      "Feature '%s' is not in the model's selected support; its effect on survival is null.",
      feature_name))
  }

  if (is.null(values)) {
    values <- preproc$feature_values[[idx_int]]
  }
  values <- as.numeric(values)

  # One profile subject per value: covariates at their means, the chosen
  # one set to the value.
  Xv <- matrix(X_mean, nrow = length(values), ncol = p, byrow = TRUE,
               dimnames = list(NULL, nm_fit))
  Xv[, idx_int] <- values
  sf <- predict(object, Xv, type = "survival")

  surv <- sf$survival
  colnames(surv) <- sprintf("%s = %s", feature_name,
                            formatC(values, digits = 6, format = "g"))

  list(time = sf$time, survival = surv)
}
