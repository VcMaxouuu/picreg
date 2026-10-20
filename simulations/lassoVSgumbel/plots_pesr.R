# ============================================================
# Figures for the Exp-Lasso vs square-root Lasso experiment.
#
# Reads the CSV written by simulation.R and produces:
#   1. PESR and prediction norm ||X (beta_hat - beta_0)||_2 / sqrt(n) vs s,
#      stacked as two rows
#   2. Gumbel location / scale estimates vs s (Exp-Lasso only)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ============================================================
# Paths
# ============================================================

script_directory <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]))))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    path <- rstudioapi::getActiveDocumentContext()$path
    if (length(path) == 1L && nzchar(path)) return(dirname(path))
  }

  getwd()
}

script_dir <- script_directory()
results_dir <- Sys.getenv("PICREG_RESULTS_DIR", unset = script_dir)
plots_dir <- Sys.getenv("PICREG_PLOTS_DIR", unset = results_dir)

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

file_paths <- sort(list.files(
  results_dir,
  pattern = "^PESR_Gumbel.*\\.csv$",
  full.names = TRUE
))

if (!length(file_paths)) {
  stop(
    "No simulation CSV files found in ", results_dir,
    ". Run simulation.R first or set PICREG_RESULTS_DIR.",
    call. = FALSE
  )
}

# ============================================================
# Data loading and validation
# ============================================================

results <- bind_rows(lapply(file_paths, function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE)
}))

required_columns <- c(
  "n", "p", "log_p_over_n", "sigma", "s", "n_sim",
  "fitted_family", "pesr", "pesr_se", "prediction_norm",
  "loc_hat", "scale_hat"
)

missing_columns <- setdiff(required_columns, names(results))
if (length(missing_columns)) {
  stop(
    "The CSV files use an obsolete or incomplete simulation format. ",
    "Run the current simulation.R again. Missing columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

unknown_families <- setdiff(unique(results$fitted_family),
                            c("gumbel", "gaussian"))
if (length(unknown_families)) {
  stop("Unknown fitted families: ", paste(unknown_families, collapse = ", "),
       call. = FALSE)
}

facet_log_ratio <- if ("log_p_over_n_target" %in% names(results)) {
  results$log_p_over_n_target
} else {
  results$log_p_over_n
}

results <- results %>%
  mutate(
    method = factor(
      fitted_family,
      levels = c("gumbel", "gaussian"),
      labels = c("exp-Lasso", "sqrt-Lasso")
    ),
    panel = paste0(
      "list(n == ", n, ", log(p/n) == ", formatC(facet_log_ratio, format = "f", digits = 1),
      ")"
    ),
    true_loc = 0,
    true_scale = sigma
  )

results$panel <- factor(results$panel, levels = unique(results$panel))

method_linetypes <- c("exp-Lasso" = "solid", "sqrt-Lasso" = "dashed")
n_panel_columns <- length(unique(results$n))
x_breaks <- pretty(range(results$s, finite = TRUE), n = 5)

# One facet row per sigma. Size the canvas from the actual grid rather than
# hard-coding it, so a single-sigma run does not produce a near-empty page.
n_panel_rows <- ceiling(nlevels(results$panel) / n_panel_columns)
figure_width <- 11
figure_height <- 2.2 + 3.4 * n_panel_rows

publication_theme <- theme_classic(base_size = 9.5, base_family = "serif") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "plain", margin = margin(b = 4)),
    axis.title.x = element_text(size = 9.5, margin = margin(t = 5)),
    axis.title.y = element_text(size = 9.5, margin = margin(r = 5)),
    axis.text = element_text(size = 8, colour = "black"),
    axis.text.y.right = element_text(colour = "grey35"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = grid::unit(2.5, "pt"),
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    panel.grid = element_blank(),
    panel.spacing.x = grid::unit(0.6, "lines"),
    panel.spacing.y = grid::unit(0.5, "lines"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8.5),
    legend.key.width = grid::unit(1.2, "cm"),
    legend.background = element_blank(),
    plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0),
    plot.margin = margin(t = 10, r = 10, b = 15, l = 10)
  )

# ============================================================
# 1. PESR and prediction norm, stacked as two facet rows
#
# Row 1: PESR. Row 2: ||X (beta_hat - beta_0)||_2 / sqrt(n). Both curves
# are black; the method is carried by the line type alone, shown in the
# legend.
# ============================================================

# Plotmath labels for the two facet rows, parsed by label_parsed below.
metric_levels <- c(
  "PESR",
  # plotmath collapses group('||', ...) to a SINGLE bar; nesting two
  # group('|', ...) is what actually produces the double-bar norm.
  "group('|', group('|', X * (hat(beta) - beta[0]), '|'), '|')[2] / sqrt(n)"
)

results_long <- bind_rows(
  results %>% mutate(metric = metric_levels[1L], value = pesr),
  results %>% mutate(metric = metric_levels[2L], value = prediction_norm)
) %>%
  mutate(metric = factor(metric, levels = metric_levels))

if (!any(is.finite(results_long$value))) {
  stop("Neither PESR nor the prediction norm contains a finite value.",
       call. = FALSE)
}

main_plot <- ggplot(
  results_long,
  aes(x = s, y = value, linetype = method, group = method)
) +
  geom_line(colour = "black", linewidth = 0.55, na.rm = TRUE) +
  facet_grid(
    rows = vars(metric), cols = vars(panel),
    scales = "free_y", switch = "y", labeller = label_parsed,
    axes = "all_x"
  ) +
  scale_x_continuous(breaks = x_breaks,
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.04))) +
  scale_linetype_manual(values = method_linetypes) +
  labs(x = expression(s^"*"), y = NULL) +
  # The method is carried by the line type; a wide key is needed for the
  # dashed pattern to read at this size.
  guides(linetype = guide_legend(
    keywidth = grid::unit(1.4, "cm"),
    override.aes = list(colour = "black", linewidth = 0.55)
  )) +
  publication_theme +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, margin = margin(r = 4))
  )

main_file <- file.path(plots_dir, "PESR_error_plot.pdf")

ggsave(
  filename = main_file,
  plot = main_plot,
  # Two metric rows (PESR, prediction norm) for every sigma row.
  width = figure_width,
  height = 2.0 + 2.8 * 2 * n_panel_rows,
  units = "in"
)

# ============================================================
# 2. Gumbel location and scale estimates
#
# Only the Exp-Lasso fit estimates these; the sqrt-Lasso rows are
# NA by construction and are dropped here.
# ============================================================

gumbel_results <- dplyr::filter(results, fitted_family == "gumbel")

gumbel_parameter_plot <- function(mean_column, target_column, y_label, title) {
  ggplot(gumbel_results, aes(x = s)) +
    geom_line(
      aes(y = .data[[target_column]]),
      colour = "black", linetype = "dotted", linewidth = 0.5, na.rm = TRUE
    ) +
    geom_line(
      aes(y = .data[[mean_column]]),
      colour = "black", linewidth = 0.7, na.rm = TRUE
    ) +
    geom_point(
      aes(y = .data[[mean_column]]),
      colour = "black", size = 0.9, na.rm = TRUE
    ) +
    facet_wrap(
      vars(panel), ncol = n_panel_columns,
      labeller = label_parsed, scales = "free_y"
    ) +
    scale_x_continuous(breaks = x_breaks,
                       expand = expansion(mult = c(0.01, 0.01))) +
    labs(
      x = expression(s^"*"), y = y_label, title = title,
      caption = "Solid line: Monte Carlo mean over repetitions. Dotted line: true value."
    ) +
    publication_theme +
    theme(legend.position = "none")
}

gumbel_plots <- list(
  gumbel_parameter_plot(
    "loc_hat",
    "true_loc",
    expression(Monte~Carlo~mean~of~hat(mu)),
    "Gumbel location estimator (true value 0)"
  ),
  gumbel_parameter_plot(
    "scale_hat",
    "true_scale",
    expression(Monte~Carlo~mean~of~hat(sigma)),
    "Gumbel scale estimator (true value sigma)"
  )
)

gumbel_file <- file.path(plots_dir, "Gumbel_LocationScale.pdf")

pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(gumbel_file, width = figure_width, height = figure_height,
           onefile = TRUE)
for (plot in gumbel_plots) print(plot)
grDevices::dev.off()

message("Saved: ", normalizePath(main_file))
message("Saved: ", normalizePath(gumbel_file))
