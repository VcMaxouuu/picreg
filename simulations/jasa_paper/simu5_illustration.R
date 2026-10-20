# =============================================================================
# Simulation 5 -- pivotal vs non-pivotal calibration (Section 4.1).
#
# An information criterion pre-sets lambda, without knowing the nuisance
# parameter. This illustration contrasts, for Poisson data, the zero-
# thresholding statistic of PIC (weighted score loss, canonical link) with the
# one of the canonical GLM (negative log-likelihood, log link, phi = id):
#   PIC       : Lambda = ||X^T (y - ybar 1)||_inf / (n sqrt(ybar))
#   canonical : Lambda = ||X^T (y - ybar 1)||_inf / n
# Both lambdas are pre-set once, with the same procedure: the (1 - alpha)-
# quantile of the statistic simulated under H0 at the reference mean
# mu_ref = 5 (picreg's mc_exact default). The PIC statistic is
# (asymptotically) pivotal, so its null law and its lambda do not depend on
# mu0; the canonical statistic scales like sqrt(mu0), so its pre-set lambda is
# too large for mu0 < mu_ref (missed signals) and too small for mu0 > mu_ref
# (false detections). One dataset is analysed for each baseline mean mu0, the
# values of mu0 differing from mu_ref.
#
# Note: estimating mu0 and recalibrating the canonical lambda on each dataset
# would amount, at the null fit, to normalizing the gradient by sqrt(ybar),
# i.e. to the PIC transformation; this is why it is not compared here.
#
#   design   : p = 100, n = 200, one fixed design with orthogonal, centered
#              and standardized columns (X^T X = n I, Condition 1), so that the
#              inputs are decoupled at the null fit, as for the canonical model
#              X = I_n
#   null     : y_i ~ Poisson(mu0); M_null draws per mu0 (null laws)
#   signal   : s = 3 nonzero coefficients, random positions and signs,
#              |beta_j| = c_sig / sqrt(mu0) (same signal-to-noise ratio for
#              every mu0), beta0 = log(mu0)
#   no fit   : only the gradient at the null fit is compared with lambda; an
#              input above lambda makes the origin lose its local optimality
#   lambda   : one pre-set value per method, Monte Carlo (n_simu_pdb draws)
#              at mu_ref
#
# Output: results/simu5_illustration.rds (null statistics, gradients of the
#         analysed datasets, lambdas), plotted by simu5_plot.R.
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
n <- 200
s <- 3
mu_grid <- c(2, 10, 50) # must differ from mu_ref
mu_ref <- 5 # reference mean used to pre-set both lambdas
c_sig <- 0.5 # |beta_j| = c_sig / sqrt(mu0)
M_null <- 1e4 # null draws per mu0 (null laws)

## ------------------------------ design -------------------------------------
# Orthogonal, centered, standardized columns: QR of a centered Gaussian matrix.
make_orthogonal_design <- function(n, p) {
    Q <- qr.Q(qr(make_design(n, p, rho = 0)))
    sqrt(n) * Q
}

## --------------------------- statistics ------------------------------------
# Gradient at the null fit, divided by n, for every column of Y (X centered).
grad_pic <- function(X, Y) {
  yb <- colMeans(Y)
  sweep(crossprod(X, Y), 2, nrow(X) * sqrt(pmax(yb, 1e-300)), "/") * (yb > 0)
}
grad_can <- function(X, Y) crossprod(X, Y) / nrow(X)

null_stats <- function(X, mu, m) {
  Y <- matrix(rpois(nrow(X) * m, mu), nrow(X), m)
  list(
    pic = apply(abs(grad_pic(X, Y)), 2, max),
    can = apply(abs(grad_can(X, Y)), 2, max)
  )
}

## ----------------------------- design and lambdas --------------------------
X <- make_orthogonal_design(n, p)
stopifnot(max(abs(colMeans(X))) < 1e-10, max(abs(crossprod(X) / n - diag(p))) < 1e-10)

# PIC: picreg's own calibration (mc_exact simulates at mu = 5 = mu_ref).
lam_pic <- lambda_pdb(X, family = "poisson", method = "mc_exact",
                      n_simu = n_simu_pdb, alpha = alpha)$value
# Canonical: same Monte Carlo at mu_ref, with the canonical statistic.
lam_can <- unname(quantile(null_stats(X, mu_ref, n_simu_pdb)$can, 1 - alpha))
cat(sprintf("lambda_PIC = %.4f   lambda_can = %.4f (mu_ref = %g)\n", lam_pic, lam_can, mu_ref))

## ------------------------------ null laws ----------------------------------
nulls <- do.call(rbind, lapply(mu_grid, function(mu) {
  st <- null_stats(X, mu, M_null)
  rbind(
    data.frame(mu0 = mu, method = "pic", stat = st$pic),
    data.frame(mu0 = mu, method = "can", stat = st$can)
  )
}))

## ------------------------ one dataset per mu0 -----------------------------
# Same support and signs for the three means, for comparability.
supp <- sort(sample.int(p, s))
signs <- sample(c(-1, 1), s, replace = TRUE)
gen_y <- function(mu, support, sgn) {
  beta <- numeric(p)
  beta[support] <- sgn * c_sig / sqrt(mu)
  rpois(n, mu * exp(drop(X %*% beta)))
}
examples <- do.call(rbind, lapply(mu_grid, function(mu) {
  y <- gen_y(mu, supp, signs)
  rbind(
    data.frame(mu0 = mu, method = "pic", j = seq_len(p), grad = abs(drop(grad_pic(X, cbind(y)))),
               signal = seq_len(p) %in% supp, lambda = lam_pic),
    data.frame(mu0 = mu, method = "can", j = seq_len(p), grad = abs(drop(grad_can(X, cbind(y)))),
               signal = seq_len(p) %in% supp, lambda = lam_can)
  )
}))

## --------------------------- analysis -------------------------------------
# For each dataset: lambda, and numbers of signal and noise inputs above it.
analysis <- do.call(rbind, lapply(split(examples, list(examples$method, examples$mu0)), function(d) {
  above <- d$grad > d$lambda
  data.frame(mu0 = d$mu0[1], method = d$method[1], lambda = d$lambda[1],
             signal_above = sum(above & d$signal), noise_above = sum(above & !d$signal))
}))
print(analysis[order(analysis$mu0, analysis$method), ], row.names = FALSE, digits = 4)

## ------------------------------ output -------------------------------------
dir.create(file.path(script_dir, "results"), showWarnings = FALSE)
saveRDS(list(nulls = nulls, examples = examples, analysis = analysis,
             lambda = c(pic = lam_pic, can = lam_can), mu_ref = mu_ref, n = n, p = p),
        file.path(script_dir, "results/simu5_illustration.rds"))
cat("Saved ->", normalizePath(file.path(script_dir, "results/simu5_illustration.rds")), "\n")
