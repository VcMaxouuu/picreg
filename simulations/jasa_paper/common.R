# =============================================================================
# Shared helpers for the JASA paper simulations.
#
# Sourced by simu1_phase_transition.R, simu2_weak_signals.R and
# simu3_real_data.R.
#
# A case uses a fixed standardized Gaussian design X (n x p) with Toeplitz
# correlation Sigma_jk = rho^|j-k| (rho = 0: independent columns). Each Monte
# Carlo replicate:
#   - draws beta with s nonzero entries at uniformly random positions,
#     random signs and magnitudes |beta_j| ~ U(amp[1], amp[2])
#     (amp[1] = amp[2]: fixed magnitude),
#   - generates y from the family with intercept beta0,
#   - fits every selector on the SAME dataset and records PESR, FDR, TPR.
#
# Selectors:
#   pic_l1   : PIC, lasso penalty, at the pivotal detection boundary
#   pic_scad : PIC, scad penalty,  at the pivotal detection boundary
#   pic_l0   : (Gaussian only) forward selection on log(RSS_s) + lambda s, with
#              lambda the (1-alpha)-quantile of log(RSS_0/RSS_1) under the null
#   glmnet   : cv.glmnet (10 folds) at lambda.min
#   bic      : forward selection, BIC  = -2 loglik + s log n
#   ebic     : forward selection, EBIC = BIC + 2 gamma log choose(p, s)
# =============================================================================

library(glmnet)
library(parallel)
library(picreg)

## ---------------------- parameters common to all cases ---------------------
M <- 200 # Monte Carlo replicates per setting
alpha <- 0.05 # level for the PIC detection boundaries
gamma_ebic <- 0.5 # EBIC gamma
maxsteps <- 40 # cap on forward-selection steps
n_simu_pdb <- 2000L # Monte Carlo draws for the detection boundaries (M)
ncores <- max(1, detectCores() - 1) # set to 1 on Windows

all_methods <- c("pic_l1", "pic_scad", "pic_l0", "glmnet", "bic", "ebic")

## ------------------------------ design -------------------------------------
# Same convention as picreg's internal standardization (population sd).
standardize_design <- function(X) {
  X <- sweep(X, 2, colMeans(X), "-")
  sX <- sqrt(colMeans(X^2))
  sX[sX == 0] <- 1
  sweep(X, 2, sX, "/")
}

make_design <- function(n, p, rho = 0) {
  X <- matrix(rnorm(n * p), n, p)
  if (rho != 0) X <- X %*% chol(rho^abs(outer(seq_len(p), seq_len(p), "-")))
  standardize_design(X)
}

## ------------------------- data generation ---------------------------------
# Inverse noise sd at the baseline mean mu0, used to put the signal on the
# same scale across families:
#   gaussian : 1
#   binomial : 1 / sqrt(mu0 (1 - mu0))   (= 2 at mu0 = 0.5)
#   poisson  : 1 / sqrt(mu0)             (~ 0.45 at mu0 = 5)
sigma_eff <- function(family, mu0) {
  switch(family,
         gaussian = 1,
         binomial = 1 / sqrt(mu0 * (1 - mu0)),
         poisson  = 1 / sqrt(mu0)
  )
}

draw_beta <- function(p, s, amp) {
  beta <- numeric(p)
  if (s > 0) {
    idx <- sample.int(p, s)
    beta[idx] <- runif(s, amp[1], amp[2]) * sample(c(-1, 1), s, replace = TRUE)
  }
  beta
}

# cfg: list(family, beta0, amp)
generate_data <- function(X, s, cfg) {
  beta <- draw_beta(ncol(X), s, cfg$amp)
  eta <- cfg$beta0 + as.numeric(X %*% beta)
  n <- nrow(X)
  y <- switch(cfg$family,
              gaussian = eta + rnorm(n),
              binomial = rbinom(n, 1, plogis(eta)),
              poisson  = rpois(n, exp(eta))
  )
  list(y = y, support = which(beta != 0))
}

## --------------------------- forward selection -----------------------------
# Generic greedy forward selection shared by BIC, EBIC and PIC-l0: the path
# (order in which variables enter) does not depend on the penalty, so it is
# computed once and each criterion stops at its first non-improving step.
#
#   step(S) : returns list(j, dev) for the best addition to S, or NULL;
#             dev is -2 loglik (up to a constant) of the model S + j.
#   dev0    : -2 loglik of the intercept-only model.
#   pens    : named list of functions s -> penalty on the -2 loglik scale.
forward_select <- function(step, dev0, pens, maxsteps) {
  S <- integer(0)
  D <- dev0
  active <- setNames(rep(TRUE, length(pens)), names(pens))
  sel <- setNames(vector("list", length(pens)), names(pens))
  while (any(active) && length(S) < maxsteps) {
    st <- step(S)
    if (is.null(st) || !is.finite(st$dev)) break
    k <- length(S)
    for (cr in names(pens)[active]) {
      if (st$dev + pens[[cr]](k + 1) >= D + pens[[cr]](k)) {
        active[cr] <- FALSE
        sel[[cr]] <- S
      }
    }
    S <- c(S, st$j)
    D <- st$dev
  }
  for (cr in names(pens)[active]) sel[[cr]] <- S
  sel
}

# Exact OLS forward selection via Gram-Schmidt updates (X centered, intercept
# always in). Keeps the residual and the design orthogonalized against the
# current model, so each step costs O(np).
gaussian_stepper <- function(X, y) {
  n <- nrow(X)
  r <- y - mean(y)
  Xt <- X
  rss <- sum(r^2)
  function(S) {
    nrm2 <- colSums(Xt^2)
    score <- drop(crossprod(Xt, r))^2 / nrm2
    score[S] <- -Inf
    score[nrm2 < 1e-8 * n] <- -Inf
    j <- which.max(score)
    if (!is.finite(score[j])) {
      return(NULL)
    }
    q <- Xt[, j] / sqrt(nrm2[j])
    r <<- r - q * sum(q * r)
    Xt <<- Xt - tcrossprod(q, drop(crossprod(q, Xt)))
    rss <<- sum(r^2)
    list(j = j, dev = n * log(rss / n))
  }
}

glm_stepper <- function(X, y, fam) {
  p <- ncol(X)
  function(S) {
    best <- NULL
    for (j in setdiff(seq_len(p), S)) {
      fit <- tryCatch(
        suppressWarnings(glm.fit(cbind(1, X[, c(S, j), drop = FALSE]), y,
                                 family = fam
        )),
        error = function(e) NULL
      )
      if (is.null(fit) || !is.finite(fit$deviance)) next
      if (is.null(best) || fit$deviance < best$dev) {
        best <- list(j = j, dev = fit$deviance)
      }
    }
    best
  }
}

# Returns list(bic = ..., ebic = ..., [pic_l0 = ...]) of selected supports.
fit_forward <- function(X, y, family, lam_l0 = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  pens <- list(
    bic  = function(s) s * log(n),
    ebic = function(s) s * log(n) + 2 * gamma_ebic * lchoose(p, s)
  )
  if (family == "gaussian") {
    # PIC-l0: log(RSS_s) + lam s  <=>  n log(RSS_s / n) + n lam s
    if (!is.null(lam_l0)) pens$pic_l0 <- function(s) n * lam_l0 * s
    step <- gaussian_stepper(X, y)
    dev0 <- n * log(sum((y - mean(y))^2) / n)
    ms <- min(maxsteps, n - 2)
  } else {
    fam <- switch(family,
                  binomial = binomial(),
                  poisson = poisson()
    )
    step <- glm_stepper(X, y, fam)
    dev0 <- suppressWarnings(glm.fit(matrix(1, n, 1), y, family = fam)$deviance)
    ms <- maxsteps
  }
  forward_select(step, dev0, pens, ms)
}

# PIC-l0 detection boundary: (1 - alpha)-quantile of log(RSS_0 / RSS_1) under
# the null y ~ N(mu, sigma^2 I). Pivotal, so mu = 0, sigma = 1 suffices. With X
# centered, RSS_1 = RSS_0 - max_j (x_j' r)^2 / ||x_j||^2.
lambda_pic_l0 <- function(X, n_simu, alpha) {
  n <- nrow(X)
  Y <- matrix(rnorm(n * n_simu), n, n_simu)
  R <- sweep(Y, 2, colMeans(Y), "-")
  rss0 <- colSums(R^2)
  gain <- apply(crossprod(X, R)^2 / colSums(X^2), 2, max)
  unname(quantile(log(rss0 / (rss0 - gain)), 1 - alpha))
}

## ----------------------------- other selectors -----------------------------
fit_glmnet <- function(X, y, family) {
  cvfit <- tryCatch(
    suppressWarnings(cv.glmnet(X, y, family = family, nfolds = 10)),
    error = function(e) NULL
  )
  if (is.null(cvfit)) {
    return(integer(0))
  }
  b <- as.numeric(coef(cvfit, s = "lambda.min"))[-1]
  which(b != 0)
}

fit_pic <- function(X, y, family, penalty, lam) {
  fit <- tryCatch(
    suppressWarnings(pic(X, y, family = family, penalty = penalty, lambda = lam)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(integer(0))
  }
  which(fit$beta != 0)
}

## ----------------------------- metrics -------------------------------------
metrics <- function(Shat, Strue) {
  Shat <- unique(Shat)
  Strue <- unique(Strue)
  tp <- length(intersect(Shat, Strue))
  fp <- length(setdiff(Shat, Strue))
  nsel <- length(Shat)
  ntrue <- length(Strue)
  c(
    PESR = as.numeric(setequal(Shat, Strue)),
    FDR = if (nsel == 0) 0 else fp / nsel,
    TPR = if (ntrue == 0) 1 else tp / ntrue,
    SIZE = nsel
  )
}

## ----------------------------- one design ----------------------------------
# Detection boundaries, computed once per design X (lambda depends on X only).
design_lambdas <- function(X, family, methods) {
  list(
    lam = lambda_pdb(X,
                     family = family, method = "mc_exact",
                     n_simu = n_simu_pdb, alpha = alpha
    )$value,
    lam_l0 = if (family == "gaussian" && "pic_l0" %in% methods) {
      lambda_pic_l0(X, n_simu_pdb, alpha)
    }
  )
}

# One Monte Carlo replicate; all selectors fitted on the SAME dataset.
run_one <- function(s, X, cfg, lams, methods) {
  fam <- cfg$family
  dat <- generate_data(X, s, cfg)
  fw <- if (any(c("bic", "ebic", "pic_l0") %in% methods)) {
    fit_forward(X, dat$y, fam, lams$lam_l0)
  }
  sel <- list(
    pic_l1 = if ("pic_l1" %in% methods) fit_pic(X, dat$y, fam, "lasso", lams$lam),
    pic_scad = if ("pic_scad" %in% methods) fit_pic(X, dat$y, fam, "scad", lams$lam),
    pic_l0 = if ("pic_l0" %in% methods) fw$pic_l0,
    glmnet = if ("glmnet" %in% methods) fit_glmnet(X, dat$y, fam),
    bic = if ("bic" %in% methods) fw$bic,
    ebic = if ("ebic" %in% methods) fw$ebic
  )
  sel <- sel[methods[methods %in% names(sel)]]
  sel <- sel[!vapply(sel, is.null, logical(1))]
  do.call(rbind, lapply(sel, metrics, Strue = dat$support))
}

# M replicates at sparsity s. Every replicate gets its own seed, so results do
# not depend on the number of cores. Returns a data.frame with one row per
# method: mean and Monte Carlo standard error of PESR, FDR, TPR and SIZE
# (number of selected variables).
run_reps <- function(s, X, cfg, lams, methods, seed, tag = "") {
  sims <- mclapply(seq_len(M), function(m) {
    set.seed(seed * 1e6 + s * 1e3 + m)
    run_one(s, X, cfg, lams, methods)
  }, mc.cores = ncores)
  failed <- vapply(sims, inherits, logical(1), "try-error")
  if (any(failed)) {
    stop(sprintf(
      "[%s] s = %d: %d replicates failed:\n%s",
      tag, s, sum(failed), sims[[which(failed)[1]]]
    ))
  }
  arr <- simplify2array(sims) # methods x 4 metrics x M
  mn <- apply(arr, c(1, 2), mean)
  se <- apply(arr, c(1, 2), sd) / sqrt(M)
  data.frame(
    model = rownames(mn),
    PESR = mn[, "PESR"], FDR = mn[, "FDR"], TPR = mn[, "TPR"], SIZE = mn[, "SIZE"],
    PESR_se = se[, "PESR"], FDR_se = se[, "FDR"], TPR_se = se[, "TPR"],
    SIZE_se = se[, "SIZE"],
    row.names = NULL
  )
}

## ------------------------------ output -------------------------------------
save_results <- function(results, out_dir, name) {
  num <- vapply(results, is.numeric, logical(1)) &
    names(results) %in% c("PESR", "FDR", "TPR", "SIZE",
                          "PESR_se", "FDR_se", "TPR_se", "SIZE_se")
  results[num] <- round(results[num], 3)
  dir.create(file.path(out_dir, "results"), showWarnings = FALSE)
  saveRDS(results, file.path(out_dir, sprintf("results/%s.rds", name)))
  out_csv <- file.path(out_dir, sprintf("results/%s.csv", name))
  write.csv(results, out_csv, row.names = FALSE)
  cat("Saved ->", normalizePath(out_csv), "\n\n")
}
