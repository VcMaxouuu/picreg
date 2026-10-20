# =============================================================================
# Simulation 5 figure, 2 x 4: rows PIC (top) and canonical GLM (bottom).
#   column 1   : null law of the zero-thresholding statistic for each mu0, with
#                the lambda pre-set at mu_ref (vertical line)
#   columns 2-4: one dataset with s = 3 signals for each mu0: absolute
#                gradient at the null fit for each input (signals in black),
#                with the same pre-set lambda (horizontal line); circles mark
#                the inputs above lambda. The panels of a row share their y-axis.
# Reads results/simu5_illustration.rds written by simu5_illustration.R.
# =============================================================================
library(ggplot2)
library(patchwork)

script_dir <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) return(dirname(normalizePath(f[1])))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  getwd()
})
res_dir <- file.path(script_dir, "results")
res <- readRDS(file.path(res_dir, "simu5_illustration.rds"))

mu_grid <- sort(unique(res$nulls$mu0))
row_lab <- c(pic = "PIC (weighted score loss)", can = "Canonical GLM (log-likelihood)")
mu_lty  <- setNames(c("solid", "dashed", "dotted")[seq_along(mu_grid)], mu_grid)

base_theme <- theme_classic(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        legend.position = "none")

null_panel <- function(m) {
  d <- res$nulls[res$nulls$method == m, ]
  d$mu0 <- factor(d$mu0, levels = mu_grid)
  ggplot(d, aes(stat, linetype = mu0)) +
    geom_density(linewidth = 0.5, adjust = 1.2, key_glyph = "path") +
    geom_vline(xintercept = res$lambda[[m]], colour = "black", linewidth = 0.4) +
    scale_linetype_manual(values = mu_lty, name = expression(mu[0])) +
    labs(x = expression(Lambda), y = row_lab[[m]],
         title = if (m == "pic") "Null law of the statistic" else NULL) +
    base_theme + theme(legend.position = c(0.8, 0.75), legend.key.width = unit(0.8, "cm"))
}

# common y-axis within each row (same lambda, same scale for every mu0)
ylims <- sapply(c("pic", "can"), function(m) {
  e <- res$examples[res$examples$method == m, ]
  c(0, max(e$grad, e$lambda))
}, simplify = FALSE)

grad_panel <- function(m, mu) {
  d <- res$examples[res$examples$method == m & res$examples$mu0 == mu, ]
  lam <- d$lambda[1]
  d$above <- d$grad > lam
  ggplot(d, aes(j, grad)) +
    geom_hline(yintercept = lam, linetype = "dashed", linewidth = 0.4) +
    geom_point(data = d[!d$signal, ], colour = "grey60", size = 0.8) +
    geom_point(data = d[d$signal, ], colour = "black", size = 2, shape = 17) +
    geom_point(data = d[d$above, ], shape = 1, size = 3.6, stroke = 0.5) +
    coord_cartesian(ylim = ylims[[m]]) +
    labs(x = "input  j", y = if (mu == mu_grid[1]) "|gradient|" else NULL,
         title = if (m == "pic") bquote(mu[0] == .(mu)) else NULL) +
    base_theme
}

rows <- lapply(c("pic", "can"), function(m) {
  wrap_plots(c(list(null_panel(m)), lapply(mu_grid, function(mu) grad_panel(m, mu))), nrow = 1)
})
fig <- rows[[1]] / rows[[2]]

out <- file.path(res_dir, "simu5_illustration.pdf")
ggsave(out, fig, width = 13, height = 6)
cat("Saved ->", normalizePath(out), "\n")
