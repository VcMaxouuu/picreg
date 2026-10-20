library(survival)
library(glmnet)
library(picreg)

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

## ----- helper: design from all columns except `drop` -----------------------
make_clean <- function(df, drop, time, event) {
  covs <- setdiff(names(df), drop)
  sub  <- df[, covs, drop = FALSE]
  ok   <- complete.cases(sub)                      # drop incomplete rows
  sub  <- sub[ok, , drop = FALSE]

  X <- model.matrix(~ ., data = sub)               # includes "(Intercept)"
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]

  list(X = X, time = time[ok], event = event[ok])
}

load_pbc <- function() {
  d <- survival::pbc
  d$event <- as.integer(d$status == 2)             # 0=cens,1=transplant,2=dead
  make_clean(d, drop = c("id", "time", "status", "event"),
             time = d$time, event = d$event)
}

## ----- load and standardize PBC data ---------------------------------------
data <- load_pbc()

X <- data$X
X <- sweep(X, 2, colMeans(X), "-")
sX <- sqrt(colMeans(X^2)); sX[sX == 0] <- 1
X <- sweep(X, 2, sX, "/")

y     <- data$time
delta <- data$event
surv  <- Surv(y, delta)
p     <- ncol(X)
varnames <- colnames(X)

## ----- (1) PIC-lasso -------------------------------------------------------
fit_pic_lasso <- tryCatch(
  suppressWarnings(pic(X, cbind(time = y, status = delta),
                       family = "cox", penalty = "lasso")),
  error = function(e) NULL
)
beta_pic_lasso <- if (!is.null(fit_pic_lasso)) as.numeric(fit_pic_lasso$beta) else rep(NA, p)

## ----- (2) PIC-scad --------------------------------------------------------
fit_pic_scad <- tryCatch(
  suppressWarnings(pic(X, cbind(time = y, status = delta),
                       family = "cox", penalty = "scad")),
  error = function(e) NULL
)
beta_pic_scad <- if (!is.null(fit_pic_scad)) as.numeric(fit_pic_scad$beta) else rep(NA, p)

## ----- (3) BIC forward selection -------------------------------------------
fit_bic_forward <- function(X, y, delta, maxsteps = 50) {
  p <- ncol(X); d <- sum(delta)
  if (d < 2) return(integer(0))
  surv <- Surv(y, delta)

  L0 <- tryCatch(suppressWarnings(coxph(surv ~ 1)$loglik[1]),
                 error = function(e) NA_real_)
  if (is.na(L0)) return(integer(0))

  S <- integer(0); BIC_cur <- -2 * L0
  repeat {
    cand <- setdiff(seq_len(p), S)
    if (length(cand) == 0 || length(S) >= maxsteps) break
    bestj <- NA_integer_; bestBIC <- BIC_cur
    for (j in cand) {
      fit <- tryCatch(
        suppressWarnings(coxph(surv ~ X[, c(S, j), drop = FALSE])),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      BICj <- -2 * fit$loglik[2] + log(d) * (length(S) + 1)
      if (is.finite(BICj) && BICj < bestBIC) { bestBIC <- BICj; bestj <- j }
    }
    if (is.na(bestj)) break
    S <- c(S, bestj); BIC_cur <- bestBIC
  }
  S
}

S_bic <- fit_bic_forward(X, y, delta)
beta_bic <- rep(0, p)
if (length(S_bic) > 0) {
  refit <- tryCatch(suppressWarnings(coxph(surv ~ X[, S_bic, drop = FALSE])),
                    error = function(e) NULL)
  if (!is.null(refit)) beta_bic[S_bic] <- coef(refit)
}

## ----- (4) cv.glmnet -------------------------------------------------------
cvfit <- tryCatch(
  suppressWarnings(cv.glmnet(X, surv, family = "cox")),
  error = function(e) NULL
)
beta_glmnet <- if (!is.null(cvfit)) as.numeric(coef(cvfit, s = "lambda.min")) else rep(NA, p)

## ----- (5) Cox + Benjamini-Hochberg ----------------------------------------
alpha <- 0.05
beta_bh <- rep(0, p)
if (sum(delta) > p) {
  fit_full <- tryCatch(suppressWarnings(coxph(surv ~ X)), error = function(e) NULL)
  if (!is.null(fit_full)) {
    cf <- coef(fit_full)
    if (all(is.finite(cf)) && max(abs(cf)) <= 50) {
      sm   <- summary(fit_full)$coefficients
      padj <- p.adjust(sm[, "Pr(>|z|)"], method = "BH")
      S_bh <- which(padj < alpha)
      if (length(S_bh) > 0) beta_bh[S_bh] <- cf[S_bh]
    }
  }
}

## ----- (6) MLE (full Cox, all variables) -----------------------------------
beta_mle <- rep(NA, p)
if (sum(delta) > p) {
  fit_mle <- tryCatch(suppressWarnings(coxph(surv ~ X)), error = function(e) NULL)
  if (!is.null(fit_mle)) beta_mle <- as.numeric(coef(fit_mle))
}

## ----- assemble results table ----------------------------------------------
results <- data.frame(
  Variable   = varnames,
  PIC_lasso  = round(beta_pic_lasso, 3),
  PIC_scad   = round(beta_pic_scad, 3),
  BIC        = round(beta_bic, 3),
  cv.glmnet  = round(beta_glmnet, 3),
  Cox_BH     = round(beta_bh, 3),
  MLE        = round(beta_mle, 3),
  stringsAsFactors = FALSE
)

print(results, row.names = FALSE)

## ----- save to CSV ---------------------------------------------------------
out_path <- file.path(script_dir, "results/pbc_coefficients.csv")
write.csv(results, out_path, row.names = FALSE)
cat("\nSaved ->", normalizePath(out_path), "\n")

## ----- composite plot: baselines + feature effects (PIC-scad) -------------
if (!is.null(fit_pic_scad) && inherits(fit_pic_scad, "pic.cox")) {

  selected_idx <- which(fit_pic_scad$beta != 0)
  selected_names <- varnames[selected_idx]
  n_selected <- length(selected_idx)

  # MLE baseline functions (full Cox model)
  has_mle <- !is.null(fit_mle) && all(is.finite(beta_mle))
  if (has_mle) {
    bl_mle <- baseline_functions(y, delta, as.numeric(X %*% beta_mle))
  }

  if (n_selected > 0) {
    # Layout: 3 columns, baselines + feature effects
    n_panels <- 2L + n_selected
    nc <- min(3L, n_panels)
    nr <- ceiling(n_panels / nc)

    pdf(file.path(script_dir, "results/pbc_pic_scad_plots.pdf"),
        width = 3 * nc, height = 2.8 * nr)

    par(mfrow = c(nr, nc),
        mar = c(4, 4, 2.8, 1),
        oma = c(0, 0, 2.5, 0),
        mgp = c(2.3, 0.55, 0),
        las = 1, bty = "n", tcl = -0.3,
        cex.main = 1.1, cex.lab = 1.05, cex.axis = 0.95,
        font.main = 1)

    pic_col <- "black"
    mle_col <- "firebrick3"

    # --- Panel 1: baseline cumulative hazard ---
    bh <- fit_pic_scad$baseline_cumulative_hazard
    ylim_h <- range(bh$cumulative_hazard)
    if (has_mle) ylim_h <- range(c(ylim_h, bl_mle$cumulative_hazard))
    plot(bh$time, bh$cumulative_hazard,
         type = "s", lwd = 1.8, col = pic_col,
         ylim = ylim_h,
         xlab = "time", ylab = expression(H[0](t)),
         main = "Baseline cumulative hazard", axes = FALSE)
    if (has_mle) {
      lines(bl_mle$time, bl_mle$cumulative_hazard,
            type = "s", lwd = 1.4, col = mle_col, lty = 2)
    }
    axis(1, lwd = 0, lwd.ticks = 0.8)
    axis(2, lwd = 0, lwd.ticks = 0.8)
    box(bty = "l", lwd = 0.8)
    legend("topleft",
           legend = if (has_mle) c("PIC-scad", "MLE") else "PIC-scad",
           col    = if (has_mle) c(pic_col, mle_col) else pic_col,
           lty    = if (has_mle) c(1, 2) else 1,
           lwd    = if (has_mle) c(1.8, 1.4) else 1.8,
           bty = "n", cex = 0.95)

    # --- Panel 2: baseline survival ---
    bs <- fit_pic_scad$baseline_survival
    plot(bs$time, bs$survival,
         type = "s", lwd = 1.8, col = pic_col,
         ylim = c(0, 1),
         xlab = "time", ylab = expression(S[0](t)),
         main = "Baseline survival", axes = FALSE)
    if (has_mle) {
      lines(bl_mle$time, bl_mle$survival,
            type = "s", lwd = 1.4, col = mle_col, lty = 2)
    }
    axis(1, lwd = 0, lwd.ticks = 0.8)
    axis(2, lwd = 0, lwd.ticks = 0.8)
    box(bty = "l", lwd = 0.8)
    legend("topright",
           legend = if (has_mle) c("PIC-scad", "MLE") else "PIC-scad",
           col    = if (has_mle) c(pic_col, mle_col) else pic_col,
           lty    = if (has_mle) c(1, 2) else 1,
           lwd    = if (has_mle) c(1.8, 1.4) else 1.8,
           bty = "n", cex = 0.95)

    # --- Feature effect panels ---
    lty_pool <- c(1L, 2L, 3L, 4L, 5L, 6L)
    pch_pool <- c(16L, 17L, 15L, 18L, 8L, 4L)

    for (i in seq_along(selected_idx)) {
      sf <- feature_effects_on_survival(fit_pic_scad, selected_idx[i])
      S  <- sf$survival
      k  <- ncol(S)
      cols <- if (k == 1L) "black"
              else grDevices::gray(seq(0, 0.55, length.out = k))
      lty_vec <- lty_pool[((seq_len(k) - 1L) %% length(lty_pool)) + 1L]
      pch_vec <- pch_pool[((seq_len(k) - 1L) %% length(pch_pool)) + 1L]

      matplot(sf$time, S,
              type = "s", lty = lty_vec, lwd = 1.6, col = cols,
              ylim = c(0, 1),
              xlab = "time", ylab = "survival",
              main = selected_names[i], axes = FALSE)
      axis(1, lwd = 0, lwd.ticks = 0.8)
      axis(2, lwd = 0, lwd.ticks = 0.8)
      box(bty = "l", lwd = 0.8)

      # sparse markers
      n_marks <- min(8L, nrow(S))
      mark_idx <- unique(round(seq(1, nrow(S), length.out = n_marks)))
      for (j in seq_len(k)) {
        points(sf$time[mark_idx], S[mark_idx, j],
               pch = pch_vec[j], col = cols[j], bg = "white", cex = 0.95)
      }

      # round legend values to 3 digits
      leg_labels <- sub("= .+", "", colnames(S))
      leg_values <- as.numeric(sub(".+= ", "", colnames(S)))
      leg_labels <- paste0(leg_labels, "= ", formatC(leg_values, digits = 2, format = "g"))

      legend("topright", legend = leg_labels,
             col = cols, lty = lty_vec, lwd = 1.6,
             pch = pch_vec, bty = "n", cex = 0.85, seg.len = 2.4)
    }

    dev.off()
    cat("Plot saved ->",
        normalizePath(file.path(script_dir, "results/pbc_pic_scad_plots.pdf")),
        "\n")
  }
}
