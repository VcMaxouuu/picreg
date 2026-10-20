# =============================================================================
# Phase transition of support recovery for different Cox models.
#
# For a fixed standardized Gaussian design X (n x p), and for each sparsity
# s = 0,...,smax, run M Monte Carlo replicates. Each replicate:
#   - draws beta_s with s nonzero entries ~ U(1,3) with random sign,
#   - generates survival times  T = -log(U) / (h0 * exp(X beta_s)),  h0 = 0.1,
#   - censors to ~30% censoring,
#   - fits the selectors and records PESR, FDR, TPR.
# Metrics are averaged over the M replicates and saved, one file per case.
#
# Models:
#   (1) glmnet    : cv.glmnet(family="cox") at lambda.min
#   (2) bic       : Cox BIC via forward selection (survival::coxph), BIC on #events
#   (3) bh        : full Cox + p-values, Benjamini-Hochberg adjusted (p.adjust)
#   (4) pic_lasso : PIC-Cox, lasso penalty, at the pivotal detection boundary
#   (5) pic_scad  : PIC-Cox, scad penalty,  at the pivotal detection boundary
# =============================================================================

library(glmnet)
library(survival)
library(parallel)
library(picreg)

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

## ---------------------- parameters common to all cases ---------------------
M            <- 200             # Monte Carlo replicates per s
h0           <- 0.1             # baseline hazard
target_event <- 0.70            # 1 - censoring rate (=> ~30% censoring)
alpha        <- 0.05            # level for the BH selection
maxsteps     <- 30              # cap on forward-selection steps
n_simu_pdb   <- 5000L           # Monte Carlo draws for the PIC detection boundary
ncores       <- max(1, detectCores() - 1)   # set to 1 on Windows

## ------------------------- data generation ---------------------------------
# Returns y, delta (censoring at ~30%), and the true support.
generate_data <- function(X, s, target_event) {
  n <- nrow(X); p <- ncol(X)
  beta <- numeric(p)
  if (s > 0) {
    idx <- sample.int(p, s)
    beta[idx] <- runif(s, 1, 3) * sample(c(-1, 1), s, replace = TRUE)
  }
  support <- which(beta != 0)
  
  eta <- as.numeric(X %*% beta)               # linear predictor
  U   <- runif(n)
  T   <- -log(U) / (h0 * exp(eta))            # T ~ Exp(rate = h0 exp(eta))
  
  # calibrate an Exp censoring rate so that E[#events]/n = target_event:
  #   E[delta_i | T_i] = exp(-rate * T_i)  ->  solve mean(exp(-rate T)) = target
  g <- function(rate) mean(exp(-rate * T)) - target_event
  rate_c <- tryCatch(
    uniroot(g, interval = c(1e-8, 1e4), extendInt = "yes")$root,
    error = function(e) 1.0
  )
  C     <- rexp(n, rate = rate_c)
  y     <- pmin(T, C)
  delta <- as.numeric(T <= C)
  
  list(y = y, delta = delta, support = support)
}

## ----------------------------- selectors -----------------------------------
# (1) cross-validated Cox lasso at lambda.min
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
fit_bic_forward <- function(X, y, delta, maxsteps) {
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

# (4)/(5) PIC-Cox with lasso / scad penalty, at the pivotal detection boundary lam.
#   y must be passed to pic() as an n x 2 matrix: column 1 = time, column 2 = status.
fit_pic <- function(X, y, delta, penalty, lam) {
  ymat <- cbind(time = y, status = delta)
  fit <- tryCatch(
    suppressWarnings(pic(X, ymat, family = "cox", penalty = penalty, lambda = lam)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(integer(0))
  which(fit$beta != 0)
}

## ----------------------------- metrics -------------------------------------
metrics <- function(Shat, Strue, p) {
  Shat <- unique(Shat); Strue <- unique(Strue)
  tp    <- length(intersect(Shat, Strue))
  fp    <- length(setdiff(Shat, Strue))
  nsel  <- length(Shat); ntrue <- length(Strue)
  c(PESR = as.numeric(setequal(Shat, Strue)),
    FDR  = if (nsel  == 0) 0 else fp / nsel,
    TPR  = if (ntrue == 0) 1 else tp / ntrue)
}

## ------------------------ one Monte Carlo replicate ------------------------
# Fits all selectors on the SAME dataset for a fair comparison.
run_one <- function(s, X, lam) {
  p   <- ncol(X)
  dat <- generate_data(X, s, target_event)
  rbind(
    glmnet    = metrics(fit_glmnet(X, dat$y, dat$delta),                dat$support, p),
    bic       = metrics(fit_bic_forward(X, dat$y, dat$delta, maxsteps), dat$support, p),
    bh        = metrics(fit_bh(X, dat$y, dat$delta, alpha),             dat$support, p),
    pic_lasso = metrics(fit_pic(X, dat$y, dat$delta, "lasso", lam),     dat$support, p),
    pic_scad  = metrics(fit_pic(X, dat$y, dat$delta, "scad",  lam),     dat$support, p)
  )
}

## --------------------- run one (n, p, smax) experiment ---------------------
run_case <- function(n, p, smax, seed, tag) {
  set.seed(seed)                              # reproducible per case
  
  # fixed standardized Gaussian design for this case
  X <- matrix(rnorm(n * p), n, p)
  X <- sweep(X, 2, colMeans(X), "-")          # center
  sX <- sqrt(colMeans(X^2)); sX[sX == 0] <- 1
  X <- sweep(X, 2, sX, "/")                    # scale to unit (population) variance
  
  # PIC detection boundary, computed once on this design
  lam <- lambda_pdb(X, family = "cox", method = "mc_exact", n_simu = n_simu_pdb)$value
  cat(sprintf("[%s]  n=%d p=%d  lambda = %.5f\n", tag, n, p, lam))
  
  results <- data.frame()
  for (s in 0:smax) {
    t0   <- Sys.time()
    sims <- mclapply(seq_len(M), function(m) run_one(s, X, lam),
                     mc.cores = ncores, mc.set.seed = TRUE)
    arr      <- simplify2array(sims)          # 5 models x 3 metrics x M
    mean_mat <- apply(arr, c(1, 2), mean, na.rm = TRUE)
    
    results <- rbind(results, data.frame(
      n = n, p = p, s = s, model = rownames(mean_mat),
      PESR = mean_mat[, "PESR"], FDR = mean_mat[, "FDR"], TPR = mean_mat[, "TPR"],
      row.names = NULL))
    cat(sprintf("[%s]  s = %2d done in %5.1f s\n", tag, s,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  
  results[c("PESR", "FDR", "TPR")] <- round(results[c("PESR", "FDR", "TPR")], 3)
  out_rds <- file.path(script_dir, sprintf("results/phase_transition_%s.rds", tag))
  out_csv <- file.path(script_dir, sprintf("results/phase_transition_%s.csv", tag))
  saveRDS(results, out_rds)
  write.csv(results, out_csv, row.names = FALSE)
  cat("Saved ->", normalizePath(out_csv), "\n\n")
  invisible(results)
}

## ------------------------------ run both cases -----------------------------
res_A <- run_case(n = 100, p = 100, smax = 13, seed = 1, tag = "n100_p100")
res_B <- run_case(n = 250, p = 50,  smax = 20, seed = 2, tag = "n250_p50")

cat("All cases done.\n")