# =============================================================================
# Simulation 1 -- phase transition of exact support recovery.
#
# Large, FIXED coefficient magnitudes and an independent design, so that a
# drop in PESR can only come from s being too large. p is fixed and n varies:
# a larger n pushes the transition to larger s.
#
#   design  : p = 100, independent standard Gaussian columns (rho = 0),
#             standardized, fixed per (family, n)
#   signal  : s = 0, ..., 25 nonzero coefficients at random positions with
#             random signs and fixed magnitude
#               gaussian : |beta_j| = 3,    beta0 = 0,      eps ~ N(0, 1)
#               binomial : |beta_j| = 2,    beta0 = 0
#               poisson  : |beta_j| = 0.5,  beta0 = log 5
#   n       : gaussian {75, 150}; binomial 400; poisson 200
#   methods : PIC-l1, PIC-scad, PIC-l0 (Gaussian only), cv.glmnet, BIC, EBIC
#   metrics : PESR, FDR, TPR, model size (mean and MC standard error over M
#             replicates; M is set in common.R)
#
# Output: results/simu1_<family>_n<n>_p<p>.csv, plotted by simu1_plot.R.
# =============================================================================

## locate this folder (Rscript, source(), or RStudio)
script_dir <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) {
    return(dirname(normalizePath(f[1])))
  }
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) {
      return(dirname(normalizePath(of)))
    }
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  }
  getwd()
})
source(file.path(script_dir, "common.R"))

## ------------------------------ settings -----------------------------------
p <- 100
smax <- 25
run_families <- c("binomial") #)gaussian", "poisson", "binomial") # subset to run fewer

settings <- list(
  gaussian = list(
    cfg = list(family = "gaussian", beta0 = 0, amp = c(3, 3)),
    n_grid = c(75, 150), seed = 1
  ),
  binomial = list(
    cfg = list(family = "binomial", beta0 = 0, amp = c(2, 2)),
    n_grid = 400, seed = 2
  ),
  poisson = list(
    cfg = list(family = "poisson", beta0 = log(5), amp = c(0.5, 0.5)),
    n_grid = 200, seed = 3
  )
)

## --------------------- one (family, n) phase transition --------------------
run_phase_transition <- function(cfg, n, seed, tag) {
  methods <- if (cfg$family == "gaussian") all_methods else setdiff(all_methods, "pic_l0")
  set.seed(seed)
  X <- make_design(n, p, rho = 0)
  lams <- design_lambdas(X, cfg$family, methods)
  cat(sprintf(
    "[%s]  lambda = %.5f%s\n", tag, lams$lam,
    if (is.null(lams$lam_l0)) "" else sprintf("  lambda_l0 = %.5f", lams$lam_l0)
  ))
  
  results <- do.call(rbind, lapply(0:smax, function(s) {
    t0 <- Sys.time()
    r <- run_reps(s, X, cfg, lams, methods, seed = seed, tag = tag)
    cat(sprintf(
      "[%s]  s = %2d done in %5.1f s\n", tag, s,
      as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ))
    cbind(family = cfg$family, n = n, p = p, s = s, r)
  }))
  save_results(results, script_dir, sprintf("simu1_%s", tag))
  invisible(results)
}

## ------------------------------ run ----------------------------------------
for (fam in run_families) {
  st <- settings[[fam]]
  for (k in seq_along(st$n_grid)) {
    n <- st$n_grid[k]
    run_phase_transition(st$cfg, n,
                         seed = 100 * st$seed + k,
                         tag = sprintf("%s_n%d_p%d", fam, n, p)
    )
  }
}
cat("Simulation 1 done.\n")
