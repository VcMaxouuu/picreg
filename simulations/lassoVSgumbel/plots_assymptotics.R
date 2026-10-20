suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# Must match assymptotic_conv.R; sigma is constant over the run and is
# therefore not stored in the CSV.
SIGMA <- 5
EULER_GAMMA <- 0.5772156649015329

# Cramer-Rao bounds for a Gumbel(mu, sigma) sample. The oracle MLE attains
# them, so they are the reference the Exp-Lasso is measured against:
#   n * Var(mu_hat)    -> sigma^2 (1 + 6(1-gamma)^2 / pi^2)
#   n * Var(sigma_hat) -> 6 sigma^2 / pi^2
CRAMER_RAO_LOC <- SIGMA^2 * (1 + 6 * (1 - EULER_GAMMA)^2 / pi^2)
CRAMER_RAO_SCALE <- 6 * SIGMA^2 / pi^2

script_directory <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg))
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    path <- rstudioapi::getActiveDocumentContext()$path
    if (length(path) == 1L && nzchar(path)) return(dirname(path))
  }
  getwd()
}

script_dir <- script_directory()
results_dir <- Sys.getenv("PICREG_RESULTS_DIR", unset = script_dir)
plots_dir <- Sys.getenv("PICREG_PLOTS_DIR", unset = results_dir)
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

results_file <- file.path(results_dir, "asymptotic_convergence_gumbel.csv")
if (!file.exists(results_file)) stop("Results file not found: ", results_file)

results <- read.csv(results_file, stringsAsFactors = FALSE)

required <- c("n", "pesr", "mean_loc_hat", "mean_scale_hat",
              "mean_mle_loc_hat", "mean_mle_scale_hat",
              "var_loc_hat", "var_scale_hat",
              "var_mle_loc_hat", "var_mle_scale_hat")
missing <- setdiff(required, names(results))
if (length(missing))
  stop("Missing columns: ", paste(missing, collapse = ", "))

results <- results %>%
  arrange(n) %>%
  mutate(
    # Rescaled by sqrt(n): the lemma claims the exp-Lasso and the oracle MLE
    # are asymptotically EQUAL, i.e. the gap is o_P(n^(-1/2)) and therefore
    # still vanishes after this rescaling. The raw gap would only show that
    # both estimators are consistent.
    loc_diff   = sqrt(n) * (mean_loc_hat - mean_mle_loc_hat),
    scale_diff = sqrt(n) * (mean_scale_hat - mean_mle_scale_hat)
  )

metric_levels <- c(
  "PESR",
  "sqrt(n) * (hat(mu)[Exp-Lasso] - hat(mu)[MLE])",
  "sqrt(n) * (hat(sigma)[Exp-Lasso] - hat(sigma)[MLE])",
  "n %.% Var(hat(mu))",
  "n %.% Var(hat(sigma))"
)

# Row 1 carries a single curve; the variance panels of row 2 show the
# Exp-Lasso against the oracle MLE, so the series is a mapped aesthetic.
series_levels <- c("Exp-Lasso", "Oracle MLE")

results_long <- bind_rows(
  results %>% transmute(n, metric = metric_levels[1], series = series_levels[1],
                        value = pesr),
  results %>% transmute(n, metric = metric_levels[2], series = series_levels[1],
                        value = loc_diff),
  results %>% transmute(n, metric = metric_levels[3], series = series_levels[1],
                        value = scale_diff),
  # n * Var, not Var: a horizontal line then means the parametric sqrt(n)
  # rate, and the height of the plateau IS the asymptotic variance.
  results %>% transmute(n, metric = metric_levels[4], series = series_levels[1],
                        value = n * var_loc_hat),
  results %>% transmute(n, metric = metric_levels[4], series = series_levels[2],
                        value = n * var_mle_loc_hat),
  results %>% transmute(n, metric = metric_levels[5], series = series_levels[1],
                        value = n * var_scale_hat),
  results %>% transmute(n, metric = metric_levels[5], series = series_levels[2],
                        value = n * var_mle_scale_hat)
) %>%
  mutate(metric = factor(metric, levels = metric_levels),
         series = factor(series, levels = series_levels))

targets <- data.frame(
  metric = factor(metric_levels[2:5], levels = metric_levels),
  target = c(0, 0, CRAMER_RAO_LOC, CRAMER_RAO_SCALE)
)

publication_theme <- theme_classic(base_size = 9.5, base_family = "serif") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 9, face = "plain", margin = margin(b = 4)),
    axis.title.x = element_text(size = 9.5, margin = margin(t = 5)),
    axis.text = element_text(size = 8, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = grid::unit(2.5, "pt"),
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    panel.grid = element_blank(),
    panel.spacing.x = grid::unit(1.1, "lines"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8.5),
    legend.background = element_blank(),
    plot.margin = margin(t = 10, r = 12, b = 12, l = 10)
  )

convergence_plot <- ggplot(
  results_long,
  aes(x = n, y = value, linetype = series, group = series)
) +
  geom_hline(
    data = targets, aes(yintercept = target), inherit.aes = FALSE,
    linetype = "dotted", linewidth = 0.45, colour = "black"
  ) +
  geom_line(colour = "black", linewidth = 0.55, na.rm = TRUE) +
  facet_wrap(vars(metric), ncol = 3, scales = "free_y",
             labeller = label_parsed) +
  scale_linetype_manual(values = c("Exp-Lasso" = "solid",
                                   "Oracle MLE" = "dashed")) +
  guides(linetype = guide_legend(keywidth = grid::unit(1.4, "cm"))) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.06))) +
  labs(
    x = expression(n), y = NULL,
  ) +
  publication_theme

output_file <- file.path(plots_dir, "asymptotic_convergence.pdf")

# Two facet rows now: three panels on top, the two variance panels below.
ggsave(output_file, convergence_plot, width = 10, height = 6.4, units = "in")
message("Saved: ", normalizePath(output_file))