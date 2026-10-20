# =============================================================================
# Simulation 4 -- sensitivity of the PDB calibration to the nuisance parameter
# (Section 3.1; answers R1, point 2a).
#
# For the one-parameter exponential family, the null statistic Lambda is only
# asymptotically pivotal: at finite n its law still depends on the null mean
# mu0. In practice lambda_PDB is computed by Monte Carlo at an arbitrary
# reference mean mu_ref (picreg's mc_exact uses mu_ref = 0.5 for binomial and
# mu_ref = 5 for poisson), or by the Gaussian approximation (mc_gaussian).
# The Gaussian family (MSE loss with phi = sqrt) is included as a benchmark:
# Lambda is then exactly pivotal (location-scale family), so every Monte Carlo
# calibration matches the reference up to Monte Carlo error, whatever n.
# This simulation separates the two sources of error:
#   - Monte Carlo error   : spread of lambda over R repetitions of a
#                           calibration with M = n_simu_pdb = 2000 draws (shrinks with M),
#   - nuisance mismatch   : offset between a calibration at mu_ref and the
#                           reference calibration at mu0 (shrinks with n).
#
#   design      : p = 100, independent standard Gaussian columns, standardized,
#                 fixed per (family, n)
#   n           : 20, 50, 100, 200, 500
#   gaussian    : reference sd 1,      reference sd in {0.1, 10, 100} (intercept 0;
#                 Lambda is invariant to both the intercept and the sd)
#   binomial    : reference mu0 = 0.2, mu_ref in {0.05, 0.35, 0.5}
#   poisson     : reference mu0 = 1,   mu_ref in {0.1, 5, 50}
#   calibrations: reference (MC at the reference value mu0, M_reference draws), MC at each mu_ref and
#                 Gaussian approximation (R repetitions of n_simu_pdb draws
#                 each), Bonferroni closed form
#   level       : P_{mu0}(Lambda > lambda), i.e. the probability that PIC
#                 selects at least one variable under H0, averaged over the R
#                 calibrations and estimated with the reference draws
#
# Output: results/simu4_lambdas.{rds,csv} (one row per calibration draw) and
#         results/simu4_levels.csv (summary), plotted by simu4_plot.R.
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
n_grid <- c(20, 50, 100, 200, 500)
R_rep <- 100L # repetitions of each calibration (shows the Monte Carlo error)
M_reference <- 1e5 # draws for the reference quantile and for the levels
chunk <- 1e4 # draws per matrix product (memory)

# mu0 / mu_ref: the nuisance parameter, i.e. the null mean (binomial,
# poisson) or the noise sd (gaussian).
settings <- list(
    gaussian = list(mu0 = 1, mu_ref = c(0.1, 10, 100), seed = 40),
    binomial = list(mu0 = 0.2, mu_ref = c(0.05, 0.35, 0.5), seed = 41),
    poisson = list(mu0 = 1, mu_ref = c(0.1, 5, 50), seed = 42)
)

## --------------------------- null statistic --------------------------------
# Lambda = ||X^T (y - ybar 1)||_inf / (n sqrt(V)), with V = mean((y - ybar)^2)
# for gaussian (MSE loss, phi = sqrt), V = ybar (1 - ybar) for binomial and
# V = ybar for poisson; X is centered, so X^T (y - ybar 1) = X^T y. This is the
# zero-thresholding statistic of Table 1 and coincides with picreg's residuals
# at the null intercept. When ybar is on the boundary (all y_i equal)
# the MLE does not exist and PIC returns beta = 0: we set Lambda = 0.
draw_null <- function(n, m, family, mu) {
    switch(family,
        gaussian = matrix(rnorm(n * m, 0, mu), n, m),
        binomial = matrix(rbinom(n * m, 1, mu), n, m),
        poisson  = matrix(rpois(n * m, mu), n, m)
    )
}

null_stat <- function(X, Y, family) {
    n <- nrow(X)
    yb <- colMeans(Y)
    v <- switch(family,
        gaussian = colMeans(sweep(Y, 2, yb)^2),
        binomial = yb * (1 - yb),
        poisson  = yb
    )
    G <- apply(abs(crossprod(X, Y)), 2, max) / n
    ifelse(v > 0, G / sqrt(pmax(v, 1e-300)), 0)
}

sim_stats <- function(X, family, mu, m) {
    nb <- ceiling(m / chunk)
    unlist(lapply(seq_len(nb), function(b) {
        mb <- min(chunk, m - (b - 1) * chunk)
        null_stat(X, draw_null(nrow(X), mb, family, mu), family)
    }))
}

# mclapply with one seed per task; stops on the first failed task.
par_tasks <- function(n_tasks, seed0, fun) {
    out <- mclapply(seq_len(n_tasks), function(t) {
        set.seed(seed0 + t)
        fun(t)
    }, mc.cores = ncores)
    failed <- vapply(out, inherits, logical(1), "try-error")
    if (any(failed)) stop(out[[which(failed)[1]]])
    out
}

mc_lambda <- function(X, family, mu, m) {
    unname(quantile(sim_stats(X, family, mu, m), 1 - alpha))
}

## -------------------------- one (family, n) --------------------------------
run_case <- function(family, n) {
    st <- settings[[family]]
    seed0 <- st$seed * 1e7 + n * 1e4
    set.seed(seed0)
    X <- make_design(n, p, rho = 0)

    # reference draws, in parallel chunks
    n_chunks <- ceiling(M_reference / chunk)
    reference_stats <- unlist(par_tasks(n_chunks, seed0 + 5000, function(b) {
        sim_stats(X, family, st$mu0, min(chunk, M_reference - (b - 1) * chunk))
    }))
    lam_reference <- unname(quantile(reference_stats, 1 - alpha))
    level <- function(lam) mean(reference_stats > lam)

    # R repetitions of every calibration, in parallel
    cal_names <- c(sprintf("mc_%g", st$mu_ref), "gaussian")
    tasks <- expand.grid(rep = seq_len(R_rep), k = seq_along(cal_names))
    vals <- unlist(par_tasks(nrow(tasks), seed0, function(t) {
        k <- tasks$k[t]
        if (cal_names[k] == "gaussian") {
            lambda_pdb(X,
                family = family, method = "mc_gaussian",
                n_simu = n_simu_pdb, alpha = alpha
            )$value
        } else {
            mc_lambda(X, family, st$mu_ref[k], n_simu_pdb)
        }
    }))
    cal <- split(vals, factor(cal_names[tasks$k], levels = cal_names))
    lam_bonf <- lambda_pdb(X, family = family, method = "analytical", alpha = alpha)$value

    lambdas <- rbind(
        data.frame(family, n, calibration = "reference", rep = NA_integer_, lambda = lam_reference),
        do.call(rbind, lapply(names(cal), function(k) {
            data.frame(family, n, calibration = k, rep = seq_len(R_rep), lambda = cal[[k]])
        })),
        data.frame(family, n, calibration = "bonferroni", rep = NA_integer_, lambda = lam_bonf)
    )
    lambdas$mu0 <- st$mu0
    lambdas$sqrt_n_lambda <- sqrt(n) * lambdas$lambda

    # level of each calibration: mean over the R calibrations of
    # P_{mu0}(Lambda > lambda), with its Monte Carlo standard error
    levels <- do.call(rbind, lapply(split(lambdas$lambda, lambdas$calibration), function(l) {
        lv <- vapply(l, level, numeric(1))
        se_rep <- if (length(lv) > 1) sd(lv) / sqrt(length(lv)) else 0
        c(level = mean(lv), level_se = sqrt(se_rep^2 + mean(lv) * (1 - mean(lv)) / M_reference))
    }))
    levels <- data.frame(
        family, n,
        mu0 = st$mu0, calibration = rownames(levels),
        lambda_median = vapply(rownames(levels), function(k) {
            median(lambdas$lambda[lambdas$calibration == k])
        }, numeric(1)),
        levels, row.names = NULL
    )
    levels$frac_degenerate <- mean(reference_stats == 0)
    cat(sprintf(
        "[%s n = %4d]  reference sqrt(n) lambda = %.3f  levels: %s\n", family, n,
        sqrt(n) * lam_reference,
        paste(sprintf("%s %.3f", levels$calibration, levels$level), collapse = ", ")
    ))
    list(lambdas = lambdas, levels = levels)
}

## ------------------------ sanity check vs picreg ---------------------------
# picreg's mc_exact simulates at sd = 1 (gaussian), mu = 0.5 (binomial) and
# mu = 5 (poisson); our
# statistic at the same mu must give the same quantile up to Monte Carlo error.
local({
    set.seed(1)
    X <- make_design(200, p)
    for (fam in names(settings)) {
        mu_pkg <- c(gaussian = 1, binomial = 0.5, poisson = 5)[[fam]]
        a <- lambda_pdb(X, family = fam, method = "mc_exact", n_simu = 2e4, alpha = alpha)$value
        b <- mc_lambda(X, fam, mu_pkg, 2e4)
        cat(sprintf("[check %s] picreg mc_exact = %.5f, this script = %.5f\n", fam, a, b))
        stopifnot(abs(a / b - 1) < 0.02)
    }
})

## --------------------------------- run -------------------------------------
cases <- expand.grid(n = n_grid, family = names(settings), stringsAsFactors = FALSE)
out <- lapply(seq_len(nrow(cases)), function(i) run_case(cases$family[i], cases$n[i]))

lambdas <- do.call(rbind, lapply(out, `[[`, "lambdas"))
levels <- do.call(rbind, lapply(out, `[[`, "levels"))

dir.create(file.path(script_dir, "results"), showWarnings = FALSE)
saveRDS(lambdas, file.path(script_dir, "results/simu4_lambdas.rds"))
write.csv(lambdas, file.path(script_dir, "results/simu4_lambdas.csv"), row.names = FALSE)
levels$lambda_median <- signif(levels$lambda_median, 5)
levels[c("level", "level_se", "frac_degenerate")] <-
    round(levels[c("level", "level_se", "frac_degenerate")], 4)
write.csv(levels, file.path(script_dir, "results/simu4_levels.csv"), row.names = FALSE)
cat("Saved ->", normalizePath(file.path(script_dir, "results/simu4_levels.csv")), "\n")
