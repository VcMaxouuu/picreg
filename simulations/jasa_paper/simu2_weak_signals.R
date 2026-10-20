# =============================================================================
# Simulation 2 -- weak signals and correlated designs, at fixed sparsity.
#
# s is fixed below the capacity of PIC-l1, so that the difficulty comes from
# signal amplitudes close to the detection threshold and from correlation,
# not from s.
#
#   design  : p = 100, Toeplitz Sigma_jk = rho^|j-k|, rho in {0, 0.5, 0.8},
#             standardized, fixed per (family, n, rho). rho = 0 isolates the
#             effect of weak signals; rho > 0 adds that of correlation.
#   signal  : s = 5, random positions and signs,
#             |beta_j| = sigma_eff * U(0.5, 1.5), i.e.
#               gaussian (linear)   : U(0.5, 1.5)  beta0 = 0, eps ~ N(0, 1)
#               binomial (logistic) : U(1, 3)      mu0 = 0.5, beta0 = 0
#   null    : the same designs are also run with s = 0 (pure noise), to check
#             the calibration P(S_hat != empty set) ~ alpha for every rho
#   n       : gaussian 100; binomial 300 (add values to n_grid for
#             supplementary settings)
#   lambda  : PIC detection boundary computed on each design
#   methods : PIC-l1, PIC-scad, cv.glmnet, BIC, EBIC
#   metrics : TPR, FDR, model size at s = 5; P(S_hat != empty) at s = 0
#             (PESR also stored), M replicates
#
# Output (results/):
#   simu2_results.csv   one row per (family, n, rho, s, method): mean and MC se
#   simu2_table.csv     formatted table: one block per model, one row per
#                       method; per rho: TPR, FDR, Size (s = 5) and the s = 0
#                       column = proportion of replicates selecting >= 1 variable
#   simu2_table.tex     same table in LaTeX (booktabs)
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
s <- 5
rho_grid <- c(0, 0.5, 0.8)
methods <- c("pic_l1", "pic_scad", "glmnet", "bic", "ebic")
run_families <- c("gaussian", "binomial") # subset to run fewer

settings <- list(
  gaussian = list(family = "gaussian", mu0 = 0, beta0 = 0, n_grid = 100, seed = 1),
  binomial = list(family = "binomial", mu0 = 0.5, beta0 = 0, n_grid = 300, seed = 2)
)

## ------------------------------ run ----------------------------------------
# For each design: the null case s = 0 and the signal case s.
results <- list()
for (fam in run_families) {
  st <- settings[[fam]]
  cfg <- list(
    family = st$family, beta0 = st$beta0,
    amp = sigma_eff(st$family, st$mu0) * c(0.5, 1.5)
  )
  for (k in seq_along(st$n_grid)) {
    n <- st$n_grid[k]
    for (r in seq_along(rho_grid)) {
      rho <- rho_grid[r]
      seed <- 1000 * st$seed + 10 * k + r
      tag <- sprintf("%s_n%d_rho%.1f", fam, n, rho)
      set.seed(seed)
      X <- make_design(n, p, rho)
      lams <- design_lambdas(X, fam, methods)
      for (s_run in c(0, s)) {
        t0 <- Sys.time()
        res <- run_reps(s_run, X, cfg, lams, methods, seed = seed, tag = tag)
        cat(sprintf(
          "[%s]  s = %d  |beta| in [%.2f, %.2f]  lambda = %.4f  done in %.1f s\n",
          tag, s_run, cfg$amp[1], cfg$amp[2], lams$lam,
          as.numeric(difftime(Sys.time(), t0, units = "secs"))
        ))
        results[[paste(tag, s_run)]] <- cbind(
          family = fam, n = n, p = p, s = s_run,
          rho = rho, res
        )
      }
    }
  }
}
results <- do.call(rbind, results)
rownames(results) <- NULL
save_results(results, script_dir, "simu2_results")

## ------------------------------ table --------------------------------------
# Per rho: TPR, FDR and mean model size at s = 5, and the s = 0 column, i.e.
# the proportion of null replicates in which at least one variable is
# selected (= mean FDR at s = 0; the target for PIC is alpha = 0.05).
family_labels <- c(gaussian = "Linear regression", binomial = "Logistic regression")
method_labels <- c(
  pic_l1 = "PIC-l1", pic_scad = "PIC-SCAD", glmnet = "GLMNet",
  bic = "BIC", ebic = "EBIC"
)
method_tex <- c(
  pic_l1 = "PIC:$\\ell_1$", pic_scad = "PIC:SCAD", glmnet = "GLMNet",
  bic = "BIC", ebic = "EBIC"
)

blocks <- unique(results[c("family", "n")])
# Empty cell (method not run for that setting).
fmt <- function(x, digits = 2, empty = "-") {
  if (length(x)) formatC(x, format = "f", digits = digits) else empty
}

# Cells of one (block, method, rho): TPR, FDR, Size at s, null rate at s = 0.
cells_of <- function(d, m, rho) {
  x <- d[d$model == m & d$rho == rho & d$s == s, ]
  x0 <- d[d$model == m & d$rho == rho & d$s == 0, ]
  c(
    TPR = fmt(x$TPR), FDR = fmt(x$FDR), Size = fmt(x$SIZE, digits = 1),
    s0 = fmt(x0$FDR)
  )
}

# CSV: one block per (family, n), one row per method.
tab <- do.call(rbind, lapply(seq_len(nrow(blocks)), function(b) {
  d <- results[results$family == blocks$family[b] & results$n == blocks$n[b], ]
  do.call(rbind, lapply(methods, function(m) {
    row <- data.frame(
      Model = sprintf("%s (n = %d)", family_labels[[blocks$family[b]]], blocks$n[b]),
      Method = method_labels[[m]]
    )
    for (rho in rho_grid) {
      cl <- cells_of(d, m, rho)
      row[[sprintf("rho=%.1f TPR", rho)]] <- cl[["TPR"]]
      row[[sprintf("rho=%.1f FDR", rho)]] <- cl[["FDR"]]
      row[[sprintf("rho=%.1f Size", rho)]] <- cl[["Size"]]
      row[[sprintf("rho=%.1f s=0", rho)]] <- cl[["s0"]]
    }
    row
  }))
}))
out_tab <- file.path(script_dir, "results/simu2_table.csv")
write.csv(tab, out_tab, row.names = FALSE)
print(tab, row.names = FALSE)
cat("Saved ->", normalizePath(out_tab), "\n")

# LaTeX (booktabs).
nr <- length(rho_grid)
nc <- 4 # columns per rho
tex <- c(
  sprintf("\\begin{tabular}{l%s}", paste(rep(strrep("c", nc), nr), collapse = "")),
  "\\toprule",
  paste0(" & ", paste(sprintf("\\multicolumn{%d}{c}{$\\rho = %s$}", nc, rho_grid),
                      collapse = " & "
  ), " \\\\"),
  paste(sprintf("\\cmidrule(lr){%d-%d}", 2 + nc * (seq_len(nr) - 1), 1 + nc * seq_len(nr)),
        collapse = " "
  ),
  paste0("Method & ", paste(rep("TPR & FDR & Size & $s{=}0$", nr), collapse = " & "), " \\\\")
)
for (b in seq_len(nrow(blocks))) {
  d <- results[results$family == blocks$family[b] & results$n == blocks$n[b], ]
  tex <- c(
    tex, "\\midrule",
    sprintf(
      "\\multicolumn{%d}{l}{\\textit{%s ($n = %d$)}} \\\\",
      1 + nc * nr, family_labels[[blocks$family[b]]], blocks$n[b]
    )
  )
  for (m in methods) {
    cells <- unlist(lapply(rho_grid, function(rho) {
      cl <- cells_of(d, m, rho)
      ifelse(cl == "-", "--", cl)
    }))
    tex <- c(tex, paste0(method_tex[[m]], " & ", paste(cells, collapse = " & "), " \\\\"))
  }
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")
out_tex <- file.path(script_dir, "results/simu2_table.tex")
writeLines(tex, out_tex)
cat("Saved ->", normalizePath(out_tex), "\n")
cat("Simulation 2 done.\n")
