#' Plotting utilities for pic fits and diagnostics.
#'
#' All pic plots share a common minimalist aesthetic:
#'
#' * grayscale by default, with a single dark-navy accent reserved for the
#'   value the user is meant to focus on (e.g. \eqn{\hat\lambda_{PDB}});
#' * an L-shaped axis frame (no full box), short outward tick marks,
#'   sans-serif typography and tight margins;
#' * curves distinguished by line type and point shape — not colour —
#'   whenever multiple series are shown on the same panel.
#'
#' Entry points:
#'
#' * `plot(fit)` — lollipop plot of the non-zero coefficients.
#' * `plot(fit$lambda_pdb)` — histogram of the PDB null distribution
#'   with a vertical line at the selected \eqn{\hat\lambda}.
#' * [plot_baseline()] — Cox-only: baseline cumulative hazard and
#'   baseline survival, two panels.
#' * `plot(phase_transition(...))` — recovery curves.
#' * `plot(pdb_asymptotic(...))` — null-distribution convergence panels.
#'
#' @name pic_plots
NULL

# ---- shared style ---------------------------------------------------------
#
# Internal helpers used by every plot method to keep the visual contract
# uniform. The accent colour is reserved for *one* element per plot: the
# value the user should look at first.

.pic_accent <- "#08306B"   # dark navy — sparing accent
.pic_grey   <- "#5A5A5A"   # neutral grey for secondary elements
.pic_light  <- "#BFBFBF"   # subtle grey for fills / minor grid

# Set common plotting parameters. Returns the previous par() so callers
# can restore with on.exit().
.pic_par <- function(mar = c(4, 4, 2.4, 1), mfrow = NULL) {
  old <- graphics::par(no.readonly = TRUE)
  args <- list(
    mar       = mar,
    mgp       = c(2.3, 0.55, 0),
    las       = 1,
    family    = "sans",
    tcl       = -0.3,
    cex.main  = 1.0,
    font.main = 1,
    bty       = "n"
  )
  if (!is.null(mfrow)) args$mfrow <- mfrow
  do.call(graphics::par, args)
  old
}

# Draw the L-shaped frame and axes. Call after graphics::plot(..., axes = FALSE).
.pic_axes <- function(side_y = 2L) {
  graphics::axis(1, lwd = 0, lwd.ticks = 0.8)
  graphics::axis(side_y, lwd = 0, lwd.ticks = 0.8)
  graphics::box(bty = "l", lwd = 0.8)
}

# Curve style for K series on the same panel. ≤ 5: vary line type and
# point shape with a single (black) ink. > 5: switch to a grayscale ramp.
.pic_curve_style <- function(K) {
  if (K <= 5L) {
    list(
      col = rep("black", K),
      lty = c(1, 2, 4, 5, 3)[seq_len(K)],
      pch = c(19, 17, 15, 18, 1)[seq_len(K)]
    )
  } else {
    list(
      col = grDevices::grey.colors(K, start = 0.0, end = 0.7),
      lty = rep(1L, K),
      pch = rep(19, K)
    )
  }
}


# ---- coefficients lollipop ------------------------------------------------

#' Horizontal lollipop plot of the non-zero coefficients of a pic fit.
#'
#' One row per selected variable, sorted by descending absolute coefficient
#' value (largest at the top). Each variable is drawn as a horizontal
#' segment from zero to its fitted value. Sign is encoded by point shape
#' (filled circle = positive, hollow circle = negative).
#'
#' @param x A fitted `pic` object.
#' @param standardized Logical; if `TRUE`, plot the standardised
#'   coefficients used internally during fitting. Default `FALSE`
#'   (coefficients on the original scale of `X`).
#' @param max_features Optional cap on the number of features displayed
#'   (the strongest are kept).
#' @param ... Additional graphical parameters forwarded to
#'   [graphics::plot()] for the empty frame.
#' @return Invisibly returns the plotted (named) coefficient vector.
#' @export
plot.pic <- function(x, standardized = FALSE, max_features = NULL, ...) {
  co <- coef.pic(x, standardized = standardized)
  beta <- co[-1L]
  p_total <- length(beta)
  nz <- which(beta != 0)
  if (length(nz) == 0L) {
    message("No selected variables to plot.")
    return(invisible(NULL))
  }
  vals <- beta[nz]
  vals <- vals[order(abs(vals), decreasing = TRUE)]
  if (!is.null(max_features) && length(vals) > max_features) {
    vals <- vals[seq_len(max_features)]
  }

  K <- length(vals)
  pch_vec <- ifelse(vals >= 0, 19L, 1L)   # filled / hollow distinguishes sign
  y_pos   <- rev(seq_len(K))

  left_mar <- max(4, 0.55 * max(nchar(names(vals))))
  old <- .pic_par(mar = c(4, left_mar, 2.6, 1))
  on.exit(graphics::par(old))

  rng <- range(c(0, vals))
  pad <- max(diff(rng), 1e-12) * 0.05
  xlim <- c(rng[1L] - pad, rng[2L] + pad)

  graphics::plot(
    NULL,
    xlim = xlim,
    ylim = c(0.5, K + 0.5),
    xlab = if (standardized) "coefficient (standardised)" else "coefficient",
    ylab = "",
    yaxt = "n",
    axes = FALSE,
    main = sprintf("%d / %d variables selected (%s)",
                   length(nz), p_total, x$family$name),
    ...
  )
  graphics::axis(1, lwd = 0, lwd.ticks = 0.8)
  graphics::axis(2, at = y_pos, labels = names(vals),
                 las = 1L, tick = FALSE, line = -0.5)
  graphics::box(bty = "l", lwd = 0.8)

  graphics::abline(v = 0, col = .pic_grey, lwd = 0.8)
  graphics::segments(0, y_pos, vals, y_pos, col = "black", lwd = 1.4)
  graphics::points(vals, y_pos, pch = pch_vec, col = "black",
                   bg = "white", cex = 1.1)

  invisible(vals)
}


# ---- PDB null distribution ------------------------------------------------

#' Plot the PDB null distribution.
#'
#' Histogram of the simulated null gradient-norm statistics, with a
#' vertical line at the selected \eqn{\hat\lambda}.
#'
#' @param x A `pic.lambda_pdb` object (typically `fit$lambda_pdb`).
#' @param breaks Number of histogram bins (default 40).
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
plot.pic.lambda_pdb <- function(x, breaks = 40L, ...) {
  if (is.null(x$statistics)) {
    message("No simulated statistics to plot (method = ", x$method, ").")
    return(invisible(NULL))
  }

  old <- .pic_par()
  on.exit(graphics::par(old))

  h <- graphics::hist(x$statistics, breaks = breaks, plot = FALSE)
  xlim <- range(c(h$breaks, x$value))
  ylim <- c(0, max(h$density) * 1.08)

  graphics::plot(
    h, freq = FALSE,
    col = .pic_light, border = "white",
    xlim = xlim, ylim = ylim, axes = FALSE,
    main = sprintf("PDB null distribution (method = %s, alpha = %.3f)",
                   x$method, x$alpha),
    xlab = "null gradient-norm statistic",
    ylab = "density"
  )
  .pic_axes()
  graphics::abline(v = x$value, col = .pic_accent, lwd = 1.8)
  graphics::legend(
    "topright",
    legend = bquote(hat(lambda) == .(formatC(x$value, digits = 4, format = "g"))),
    col = .pic_accent, lwd = 1.8, bty = "n", cex = 0.9
  )
  invisible(x)
}


# ---- Cox baseline ---------------------------------------------------------

#' Plot Cox baseline cumulative hazard and baseline survival.
#'
#' Two-panel step plot for a fitted `pic.cox` model: the Breslow
#' baseline cumulative hazard \eqn{H_0(t)} on top and the baseline
#' survival \eqn{S_0(t) = \exp(-H_0(t))} below.
#'
#' @param model A fitted `pic.cox` object.
#' @return Invisibly returns `NULL`.
#' @export
plot_baseline <- function(model) {
  if (!inherits(model, "pic.cox"))
    stop("`model` must be a pic.cox fit.")
  if (is.null(model$baseline_survival) ||
      is.null(model$baseline_cumulative_hazard))
    stop("Baseline functions not stored on the fit; was the model fitted?")

  bh <- model$baseline_cumulative_hazard
  bs <- model$baseline_survival

  old <- .pic_par(mar = c(3.4, 4, 2.2, 1), mfrow = c(2L, 1L))
  on.exit(graphics::par(old))

  graphics::plot(
    bh$time, bh$cumulative_hazard,
    type = "s", lwd = 1.6, col = "black",
    xlab = "time", ylab = expression(H[0](t)),
    main = "Baseline cumulative hazard",
    axes = FALSE
  )
  .pic_axes()

  graphics::plot(
    bs$time, bs$survival,
    type = "s", lwd = 1.6, col = "black",
    ylim = c(0, 1),
    xlab = "time", ylab = expression(S[0](t)),
    main = "Baseline survival",
    axes = FALSE
  )
  .pic_axes()

  invisible(NULL)
}


# ---- phase transition -----------------------------------------------------

#' Phase-transition plot for a `pic.phase_transition` object.
#'
#' Plots the chosen recovery metric as a function of the sparsity
#' level `s`. Curves are distinguished by line type and point shape
#' (grayscale ramp beyond five curves).
#'
#' @param x An object returned by [phase_transition()].
#' @param metric One of `"exact_recovery"` (default), `"tpr"`, `"fdr"`.
#' @param ... Additional graphical parameters forwarded to [graphics::plot()].
#'
#' @return Invisibly returns `x`.
#' @export
plot.pic.phase_transition <- function(
  x,
  metric = c("exact_recovery", "tpr", "fdr"),
  ...
) {
  metric <- match.arg(metric)
  dat <- x$results
  cfg <- x$config

  ylab <- switch(metric,
                 exact_recovery = "PESR",
                 tpr            = "TPR",
                 fdr            = "FDR")

  penalties <- unique(dat$penalty)
  configs   <- unique(dat[c("n", "p")])

  old <- if (length(penalties) > 1L && nrow(configs) > 1L) {
    nr <- ceiling(sqrt(nrow(configs)))
    nc <- ceiling(nrow(configs) / nr)
    .pic_par(mfrow = c(nr, nc))
  } else {
    .pic_par()
  }
  on.exit(graphics::par(old))

  plot_one <- function(subdat, title, curve_var) {
    curves <- unique(subdat[[curve_var]])
    K <- length(curves)
    sty <- .pic_curve_style(K)

    graphics::plot(
      NA, NA,
      xlim = range(subdat$s),
      ylim = c(0, 1),
      xlab = "s",
      ylab = ylab,
      main = title,
      axes = FALSE,
      ...
    )
    .pic_axes()
    graphics::abline(h = c(0, 1), col = .pic_light, lwd = 0.6)

    for (i in seq_along(curves)) {
      d <- subdat[subdat[[curve_var]] == curves[i], ]
      d <- d[order(d$s), ]
      graphics::lines(
        d$s, d[[metric]],
        type = "b",
        pch  = sty$pch[i],
        col  = sty$col[i],
        lty  = sty$lty[i],
        lwd  = 1.4,
        cex  = 0.9,
        bg   = "white"
      )
    }

    graphics::legend(
      "topright",
      legend = curves,
      col    = sty$col,
      lty    = sty$lty,
      pch    = sty$pch,
      lwd    = 1.4,
      bty    = "n",
      cex    = 0.82,
      seg.len = 2.4
    )
  }

  if (length(penalties) == 1L) {
    dat$curve <- sprintf("n=%d, p=%d", dat$n, dat$p)
    plot_one(dat,
             sprintf("Phase transition — %s, %s", cfg$type, penalties[1]),
             "curve")
  } else if (nrow(configs) == 1L) {
    plot_one(dat,
             sprintf("Phase transition — %s, n=%d, p=%d",
                     cfg$type, configs$n[1], configs$p[1]),
             "penalty")
  } else {
    for (i in seq_len(nrow(configs))) {
      subdat <- dat[dat$n == configs$n[i] & dat$p == configs$p[i], ]
      plot_one(subdat,
               sprintf("%s — n=%d, p=%d", cfg$type, configs$n[i], configs$p[i]),
               "penalty")
    }
  }

  invisible(x)
}


# ---- PDB asymptotic -------------------------------------------------------

#' Diagnostic plot of the PDB asymptotic behaviour.
#'
#' Multi-panel histogram comparison of the simulated null gradient-norm
#' statistic under the `"mc_exact"` (light grey fill) and `"mc_gaussian"`
#' (dashed outline) selectors, one panel per `n` in `n_grid`. Two
#' vertical lines are added per panel:
#'
#' * \eqn{\hat\lambda_{PDB}} — solid navy, the empirical (1 - alpha)
#'   quantile of `"mc_exact"` (what pic would actually use).
#' * \eqn{\hat\lambda_{analytical}} — dashed black, the Bonferroni
#'   closed-form bound.
#'
#' @param x An object returned by [pdb_asymptotic()].
#' @param breaks Number of histogram bins (default 40).
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
plot.pic.pdb_asymptotic <- function(x, breaks = 40L, ...) {
  K <- length(x$n_grid)

  nr <- if (K == 1L) 1L else if (K <= 2L) 1L else if (K <= 4L) 2L
        else if (K <= 6L) 2L else ceiling(sqrt(K))
  nc <- ceiling(K / nr)

  old <- .pic_par(mar = c(3.7, 3.8, 2.4, 1), mfrow = c(nr, nc))
  on.exit(graphics::par(old))

  for (k in seq_len(K)) {
    se <- x$stats_exact[[k]]
    sg <- x$stats_gaussian[[k]]

    rng <- range(c(se, sg))
    pad <- diff(rng) * 0.05
    brks <- seq(rng[1L] - pad, rng[2L] + pad, length.out = breaks + 1L)

    h_e <- graphics::hist(se, breaks = brks, plot = FALSE)
    h_g <- graphics::hist(sg, breaks = brks, plot = FALSE)

    xlim <- range(c(brks, x$lambda_exact[k], x$lambda_analytical[k]))
    ylim <- c(0, max(h_e$density, h_g$density) * 1.12)

    graphics::plot(
      h_e, freq = FALSE,
      col = .pic_light, border = "white",
      xlim = xlim, ylim = ylim, axes = FALSE,
      main = sprintf("n = %d  (p = %d)", x$n_grid[k], x$p),
      xlab = "null gradient-norm statistic",
      ylab = "density"
    )
    # mc_gaussian overlay: outlined steps, no fill — visually subordinate.
    graphics::lines(
      stats::stepfun(h_g$breaks, c(0, h_g$density, 0)),
      do.points = FALSE, lty = 2, lwd = 1.1, col = "black"
    )
    .pic_axes()

    graphics::abline(v = x$lambda_exact[k],      col = .pic_accent, lwd = 1.8)
    graphics::abline(v = x$lambda_analytical[k], col = "black",      lwd = 1.1, lty = 3)

    if (k == 1L) {
      graphics::legend(
        "topright",
        legend = c("mc_exact", "mc_gaussian",
                   expression(hat(lambda)[PDB]),
                   expression(hat(lambda)[analytical])),
        fill   = c(.pic_light, NA, NA, NA),
        border = c("white",     NA, NA, NA),
        col    = c(NA, "black", .pic_accent, "black"),
        lwd    = c(NA, 1.1, 1.8, 1.1),
        lty    = c(NA, 2, 1, 3),
        bty    = "n",
        cex    = 0.78,
        seg.len = 2.4
      )
    }
  }

  invisible(x)
}
