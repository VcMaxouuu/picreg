# =============================================================================
# One loader per survival dataset. Each returns list(X, time, event).
# Covariates = all columns EXCEPT an explicit exclusion list (ids, dates,
# the survival columns themselves, and post-baseline / leakage variables).
# =============================================================================

library(survival)

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

## ----- 1. NCCTG Lung Cancer ------------------------------------------------
load_lung <- function() {
  d <- survival::lung
  d$event <- as.integer(d$status == 2)             # 1=censored, 2=dead
  make_clean(d, drop = c("time", "status", "event", "inst"),
             time = d$time, event = d$event)
}

## ----- 2. Mayo PBC (randomized subjects) -----------------------------------
load_pbc <- function() {
  d <- survival::pbc
  d$event <- as.integer(d$status == 2)             # 0=cens,1=transplant,2=dead
  make_clean(d, drop = c("id", "time", "status", "event"),
             time = d$time, event = d$event)
}

## ----- 4. German Breast Cancer Study Group ---------------------------------
load_gbsg <- function() {
  d <- survival::gbsg
  d$event <- as.integer(d$status == 1)
  make_clean(d, drop = c("pid", "rfstime", "status", "event"),
             time = d$rfstime, event = d$event)
}

load_chop <- function(path = file.path(script_dir, "ChopData.csv")) {
  d <- read.csv(path)
  make_clean(d, drop = c("time", "event", "X"),
             time = d$time, event = d$event)
}

## ----- NKI 70-gene breast cancer signature (package 'penalized') -----------
load_nki70 <- function() {
  data(nki70, package = "penalized")
  d <- nki70
  
  d$time  <- d$time
  d$event <- as.integer(d$event)
  d$Grade <- as.integer(d$Grade)
  
  make_clean(d, drop = c("time", "event"),
             time = d$time, event = d$event)
}

load_mcl <- function(path = file.path(script_dir, "mclData.csv")) {
  d <- read.csv(path)
  make_clean(d, drop = c("time", "event", "X"),
             time = d$time, event = d$event)
}

## ----- load all -------------------------------------------------
datasets <- list(lung = load_lung(), pbc = load_pbc(),
                 gbsg = load_gbsg(), chop = load_chop(),
                 nki70 = load_nki70(), mcl = load_mcl())


library(glmnet)
library(picreg)

fit_glmnet <- function(X, y, delta) {
  ysurv <- Surv(y, delta)
  cvfit <- tryCatch(
    suppressWarnings(cv.glmnet(X, ysurv, family = "cox")),
    error = function(e) NULL
  )
  if (is.null(cvfit)) return(integer(0))
  b <- as.numeric(coef(cvfit, s = "lambda.min"))
  which(b != 0)
}

# (2) BIC forward selection; BIC = -2 loglik + log(d) * (#covariates), d = #events
fit_bic_forward <- function(X, y, delta, maxsteps = 50) {
  p <- ncol(X); d <- sum(delta)
  if (d < 2) return(integer(0))
  surv <- Surv(y, delta)
  
  # null model partial loglik: BIC_0 = -2 * loglik_partial(0)
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

# (3) full Cox p-values + Benjamini-Hochberg  (only well-posed when d > p)
fit_bh <- function(X, y, delta, alpha) {
  surv <- Surv(y, delta)
  if (sum(delta) <= ncol(X)) return(integer(0))       # ill-posed: d <= p
  fit <- tryCatch(suppressWarnings(coxph(surv ~ X)), error = function(e) NULL)
  if (is.null(fit)) return(integer(0))
  cf <- coef(fit)
  if (any(!is.finite(cf)) || max(abs(cf)) > 50) return(integer(0))   # non-convergence
  sm   <- summary(fit)$coefficients
  padj <- p.adjust(sm[, "Pr(>|z|)"], method = "BH")
  which(padj < alpha)
}

fit_pic <- function(X, y, delta, penalty) {
  ymat <- cbind(time = y, status = delta)
  fit <- tryCatch(
    suppressWarnings(pic(X, ymat, family = "cox", penalty = penalty)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(integer(0))
  which(fit$beta != 0)
}


M     <- 100        # Monte Carlo splits
alpha <- 0.05       # level for BH and for the PDB
prop_train <- 0.70  # train fraction

## ----- stratified 70/30 split on the censoring indicator -------------------
stratified_split <- function(event, prop) {
  idx0 <- which(event == 0); idx1 <- which(event == 1)
  tr0  <- sample(idx0, round(prop * length(idx0)))
  tr1  <- sample(idx1, round(prop * length(idx1)))
  train <- sort(c(tr0, tr1))
  list(train = train, test = setdiff(seq_along(event), train))
}

## ----- Harrell C-index of a linear predictor on test data ------------------
# higher linear predictor = higher risk = shorter survival
c_index <- function(time, event, lp) {
  if (length(unique(lp)) <= 1) return(0.5)            # no discrimination
  cc <- survival::concordance(Surv(time, event) ~ lp, reverse = TRUE)
  as.numeric(cc$concordance)
}

## ----- fit each selector on train, return support + test linear predictor --
# Each returns list(support = <indices>, lp = <test linear predictor>).
# The lp is built by refitting an UNPENALIZED Cox on the selected support
# if the support is empty or the refit is ill-posed, lp is constant (C-index = 0.5).

refit_lp <- function(Xtr, ytr, dtr, Xte, S) {
  if (length(S) == 0) return(rep(0, nrow(Xte)))
  if (length(S) >= sum(dtr)) S <- S[seq_len(max(1, sum(dtr) - 1))]  # keep d>p
  fit <- tryCatch(
    suppressWarnings(coxph(Surv(ytr, dtr) ~ Xtr[, S, drop = FALSE])),
    error = function(e) NULL)
  if (is.null(fit) || any(!is.finite(coef(fit)))) return(rep(0, nrow(Xte)))
  as.numeric(Xte[, S, drop = FALSE] %*% coef(fit))
}

run_selectors <- function(Xtr, ytr, dtr, Xte, yte, dte) {
  out <- list()
  
  S <- fit_glmnet(Xtr, ytr, dtr)
  out$glmnet <- list(support = S, lp = refit_lp(Xtr, ytr, dtr, Xte, S))
  
  S <- fit_bic_forward(Xtr, ytr, dtr)
  out$bic <- list(support = S, lp = refit_lp(Xtr, ytr, dtr, Xte, S))
  
  S <- fit_bh(Xtr, ytr, dtr, alpha)
  out$bh <- list(support = S, lp = refit_lp(Xtr, ytr, dtr, Xte, S))
  
  S <- fit_pic(Xtr, ytr, dtr, "lasso")
  out$pic_lasso <- list(support = S, lp = refit_lp(Xtr, ytr, dtr, Xte, S))
  
  S <- fit_pic(Xtr, ytr, dtr, "scad")
  out$pic_scad <- list(support = S, lp = refit_lp(Xtr, ytr, dtr, Xte, S))
  
  out
}

## ----- one Monte Carlo split: returns a one-row data.frame -----------------
run_split <- function(ds, sim_id) {
  sp <- stratified_split(ds$event, prop_train)
  Xtr <- ds$X[sp$train, , drop = FALSE]; ytr <- ds$time[sp$train]; dtr <- ds$event[sp$train]
  Xte <- ds$X[sp$test , , drop = FALSE]; yte <- ds$time[sp$test ]; dte <- ds$event[sp$test ]
  
  res <- run_selectors(Xtr, ytr, dtr, Xte, yte, dte)
  
  row <- data.frame(sim = sim_id)
  for (m in names(res)) {
    row[[paste0(m, "_nselected")]] <- length(res[[m]]$support)
    row[[paste0(m, "_cindex")]]    <- round(c_index(yte, dte, res[[m]]$lp), 4)
  }
  row
}

## ----- run M splits for one dataset, save one CSV --------------------------
library(parallel)
ncores <- max(1, detectCores() - 1)
print(ncores)

run_dataset <- function(ds, name) {
  cat(sprintf("[%s] n=%d p=%d\n", name, nrow(ds$X), ncol(ds$X)))
  
  rows <- mclapply(seq_len(M), function(m) run_split(ds, m), mc.cores = ncores)
  
  # guard: drop splits that errored (mclapply returns a try-error for those)
  ok   <- vapply(rows, function(r) is.data.frame(r), logical(1))
  if (any(!ok)) cat(sprintf("[%s] %d split(s) failed\n", name, sum(!ok)))
  results <- do.call(rbind, rows[ok])
  
  out_csv <- file.path(script_dir, sprintf("results/results_%s.csv", name))
  write.csv(results, out_csv, row.names = FALSE)
  cat("Saved ->", normalizePath(out_csv), "\n\n")
  invisible(results)
}

## ----- run on every dataset ------------------------------------------------
for (nm in names(datasets)) {
  run_dataset(datasets[[nm]], nm)
}