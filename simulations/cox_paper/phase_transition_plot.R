# =============================================================================
# Combined figure, two independent rows:
#   row 1: n=100, p=100  -> x-axis 0..15, Cox+BH dropped (d <= p)
#   row 2: n=250, p=50   -> x-axis 0..25, Cox+BH kept
# Shared y-axis [0,1]; single legend at the bottom. The two rows share nothing on x.
# =============================================================================
library(ggplot2)
library(patchwork)
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

## ----- model display order, labels, linetypes, colours --------------------
model_levels <- c("pic_lasso", "pic_scad", "glmnet", "bic", "bh")
model_labels <- c("PIC-lasso", "PIC-scad", "cv.glmnet", "BIC forward", "Cox + BH")
model_ltypes <- c("PIC-lasso"   = "solid",
                  "PIC-scad"    = "longdash",
                  "cv.glmnet"   = "dashed",
                  "BIC forward" = "dotted",
                  "Cox + BH"    = "dotdash")
model_colour <- c("PIC-lasso"   = "black",
                  "PIC-scad"    = "black",
                  "cv.glmnet"   = "grey65",
                  "BIC forward" = "grey65",
                  "Cox + BH"    = "grey65")
metric_levels <- c("PESR", "TPR", "FDR")

## ----- read one case to long form ------------------------------------------
read_long <- function(csv_path, drop_bh) {
  res <- read.csv(csv_path)
  long <- do.call(rbind, lapply(metric_levels, function(mt)
    data.frame(s = res$s, model = res$model, metric = mt, value = res[[mt]])))
  if (drop_bh) long <- long[long$model != "bh", , drop = FALSE]
  long$model  <- factor(long$model,  levels = model_levels, labels = model_labels)
  long$metric <- factor(long$metric, levels = metric_levels)
  long
}

dat_A <- read_long(file.path(script_dir, "phase_transition_n100_p100.csv"), drop_bh = TRUE)
dat_B <- read_long(file.path(script_dir, "phase_transition_n250_p50.csv"),  drop_bh = FALSE)

## ----- one row (a 1x3 facet) with its own x range --------------------------
make_row <- function(dat, xmax, ylab, show_strip, show_xtitle, show_legend) {
  ggplot(dat, aes(s, value, linetype = model, colour = model)) +
    geom_line(linewidth = 0.7, key_glyph = "path") +
    facet_wrap(~ metric, nrow = 1) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, xmax, 5)) +
    scale_linetype_manual(values = model_ltypes, name = NULL, drop = FALSE) +
    scale_colour_manual(values = model_colour, name = NULL, drop = FALSE) +
    labs(x = if (show_xtitle) "sparsity  s" else NULL, y = ylab) +
    theme_classic(base_size = 11) +
    theme(strip.background = element_blank(),
          strip.text       = if (show_strip) element_text(size = 11) else element_blank(),
          panel.spacing    = unit(1, "lines"),
          axis.title.y     = element_text(size = 11),
          legend.position  = if (show_legend) "bottom" else "none")
}

row_top <- make_row(dat_A, xmax = 15, ylab = "n = 100,  p = 100",
                    show_strip = TRUE,  show_xtitle = FALSE, show_legend = FALSE)
row_bot <- make_row(dat_B, xmax = 25, ylab = "n = 250,  p = 50",
                    show_strip = FALSE, show_xtitle = TRUE,  show_legend = TRUE)

## ----- stack, collect a single legend at the bottom ------------------------
fig <- (row_top / row_bot) &
  theme(legend.key.width = unit(1.2, "cm"),
        legend.text      = element_text(size = 9))

ggsave(file.path(script_dir, "results/phase_transition_plot.pdf"),
       fig, width = 11.0, height = 6.5)
cat("Saved ->", normalizePath(file.path(script_dir, "results/phase_transition_plot.pdf")), "\n")