# Compare R vs C++ implementations of the PDB Monte Carlo selectors.
# Run from the project root (picreg/).

devtools::load_all()

# ---- problem size ---------------------------------------------------------
n      <- 200
p      <- 100
n_simu <- 5000L
alpha  <- 0.05

X <- scale(matrix(rnorm(n * p), n, p))   # standardised design

# ---- bench helper ---------------------------------------------------------
bench <- function(label, expr_call) {
  t0  <- Sys.time()
  res <- eval(expr_call)
  dt  <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf("  %-4s  lambda = %.6f   time = %.3f s",
              label, res$value, dt))
  invisible(list(lambda = res$value, time = dt))
}

cat(sprintf("Setup: n = %d, p = %d, n_simu = %d, alpha = %.3f\n\n",
            n, p, n_simu, alpha))

# ---- mc_exact across families --------------------------------------------
cat("=== mc_exact ===\n")

for (fam_name in c("gaussian", "binomial", "poisson", "exponential", "gumbel", "cox")) {
  fam <- picreg:::get_family(fam_name)
  cat(sprintf("[%s]\n", fam_name))
  c_res <- bench("C++", quote(lambda_pdb(X, fam, n_simu = n_simu, alpha = alpha)))
}