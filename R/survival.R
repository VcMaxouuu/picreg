#' Survival utilities for the Cox family.
#'
#' Provides Harrell's C-index, the Breslow partial log-likelihood, baseline
#' cumulative hazard / survival, and feature-effect curves.
#'
#' @name pic_survival
NULL

#' Harrell's concordance index.
#'
#' @param times Survival times (length `n`).
#' @param events Event indicators (0 or 1; length `n`).
#' @param predictions Risk scores; higher = higher risk = shorter survival.
#' @return Numeric scalar in `[0, 1]`.
#' @export
concordance_index <- function(times, events, predictions) {
  events <- as.numeric(events)
  ev_idx <- which(events == 1)
  if (length(ev_idx) == 0L) return(0.0)

  ti <- times[ev_idx]; pi_ <- predictions[ev_idx]
  tj <- times;          pj <- predictions; ej <- events

  comp <- outer(ti, tj, "<") |
    (outer(ti, tj, "==") & matrix(ej == 0, nrow = length(ti), ncol = length(tj), byrow = TRUE))
  # remove self-comparisons
  for (k in seq_along(ev_idx)) comp[k, ev_idx[k]] <- FALSE

  n_comp <- sum(comp)
  if (n_comp == 0L) return(0.0)

  conc <- outer(pi_, pj, ">")
  ties <- outer(pi_, pj, "==")
  score <- sum(ifelse(comp, ifelse(conc, 1.0, 0.0), 0.0)) +
    0.5 * sum(ifelse(comp, ifelse(ties, 1.0, 0.0), 0.0))
  score / n_comp
}

#' Breslow partial log-likelihood (negative, normalised by `n`).
#'
#' @inheritParams concordance_index
#' @export
cox_partial_log_likelihood <- function(times, events, predictions) {
  ord <- order(times, decreasing = TRUE)
  t_s <- times[ord]; p_s <- predictions[ord]; e_s <- as.numeric(events[ord])
  exp_p <- exp(p_s)
  utimes <- unique(t_s)
  G <- length(utimes)

  group_uncens <- numeric(G); group_exp <- numeric(G); group_events <- integer(G)
  for (i in seq_len(G)) {
    mask <- (t_s == utimes[i])
    group_uncens[i] <- sum(p_s[mask] * e_s[mask])
    group_exp[i]    <- sum(exp_p[mask])
    group_events[i] <- sum(e_s[mask])
  }
  cum_exp <- cumsum(group_exp)
  log_term <- ifelse(group_events > 0, group_events * log(cum_exp), 0)
  loss_contrib <- group_uncens - log_term
  -sum(loss_contrib) / length(times)
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

#' Survival curves for new data from a fitted Cox pic model.
#'
#' @param object A fitted `pic.cox` model.
#' @param newdata New design matrix (rows = subjects).
#' @return A list with `time` (length `K`) and `survival` (matrix `K x m`,
#'   one column per subject).
#' @export
predict_survival_function <- function(object, newdata) {
  if (!inherits(object, "pic.cox"))
    stop("predict_survival_function() is only defined for pic.cox models.")
  if (is.null(object$baseline_survival))
    stop("Baseline survival not stored on the fit; was the model fitted?")

  preproc <- attr(object, "preproc")
  px <- check_X(newdata, standardize_X = isTRUE(preproc$standardize),
                X_mean = preproc$X_mean, X_std = preproc$X_std)
  eta <- as.numeric(px$X %*% object$beta)
  S0  <- object$baseline_survival$survival
  t   <- object$baseline_survival$time
  surv <- outer(S0, exp(eta), `^`)
  list(time = t, survival = surv)
}
