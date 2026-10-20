# =============================================================================
# Simulation 1 figures, 2 x 3 each (columns: PESR, TPR, FDR):
#   simu1_gaussian.pdf : rows n = 75, n = 150   (p = 100)
#   simu1_glm.pdf      : rows binomial (n = 400), poisson (n = 200), p = 100
# Every method is drawn on every panel. Reads the CSVs written by
# simu1_phase_transition.R.
#
# Styling:
#   - greyscale only: methods are told apart by line type and marker, fixed
#     per method and identical in every figure; PIC methods are drawn in black
#     with thicker lines and filled markers, competitors in grey with thinner
#     lines and open markers;
#   - theme_bw with boxed strips: metrics on top, settings on the right;
#   - markers staggered per method (method k at s = k - 1, k - 1 + K, ...,
#     K = number of methods), so no two methods put a marker at the same s;
#   - vector PDF via cairo_pdf (embedded fonts, plotmath glyphs), plus a
#     600 dpi PNG with the same name.
# =============================================================================
library(ggplot2)

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

## ----- method order, labels and encodings (fixed per method) ---------------
model_levels <- c("pic_l1", "pic_scad", "pic_l0", "glmnet", "bic", "ebic")
model_labels <- list(expression("PIC:" * "\u2113"[1]), "PIC:SCAD", expression("PIC:" * "\u2113"[0]),
                     "GLMNet", "BIC", "EBIC")
model_ltypes <- c(pic_l1 = "solid", pic_scad = "44",   pic_l0 = "13",
                  glmnet = "1343",  bic      = "82",   ebic   = "2262")
model_shapes <- c(pic_l1 = 16, pic_scad = 15, pic_l0 = 17,    # filled
                  glmnet = 4,  bic = 1,       ebic = 5)       # open
# PIC methods: black, thicker lines; competitors: grey, thinner lines.
model_colour <- c(pic_l1 = "black",  pic_scad = "black",  pic_l0 = "black",
                  glmnet = "grey55", bic      = "grey55", ebic   = "grey55")
model_lwd    <- c(pic_l1 = 0.75, pic_scad = 0.75, pic_l0 = 0.75,
                  glmnet = 0.45, bic      = 0.45, ebic   = 0.45)
metric_levels <- c("PESR", "TPR", "FDR")

## ----- read cases to long form ---------------------------------------------
# cases: named character vector, row label (plotmath) -> csv tag
read_long <- function(cases) {
  do.call(rbind, lapply(names(cases), function(lab) {
    res <- read.csv(file.path(res_dir, sprintf("simu1_%s.csv", cases[[lab]])))
    do.call(rbind, lapply(metric_levels, function(mt)
      data.frame(case = lab, s = res$s, model = res$model, metric = mt,
                 value = res[[mt]])))
  })) |>
    transform(case   = factor(case, levels = names(cases)),
              model  = factor(model, levels = model_levels),
              metric = factor(metric, levels = metric_levels))
}

# Staggered markers: the k-th method present gets markers at
# s = k - 1, k - 1 + K, k - 1 + 2K, ..., with K the number of methods.
marker_points <- function(dat) {
  present <- model_levels[model_levels %in% dat$model]
  K <- length(present)
  offset <- setNames(seq_along(present) - 1L, present)
  dat[(dat$s - offset[as.character(dat$model)]) %% K == 0, ]
}

## ----- 2 x 3 grid -----------------------------------------------------------
# Metrics in boxed strips on top, settings (n, family) in boxed strips on the
# right, attached to the panels (theme_bw strips).
make_fig <- function(dat) {
  present <- model_levels[model_levels %in% dat$model]
  lab     <- model_labels[match(present, model_levels)]
  smax    <- max(dat$s)
  ggplot(dat, aes(s, value, linetype = model, shape = model, colour = model)) +
    geom_line(aes(linewidth = model)) +
    geom_point(data = marker_points(dat), size = 1.6, stroke = 0.5, fill = "white") +
    facet_grid(case ~ metric, labeller = label_parsed) +
    scale_x_continuous(breaks = seq(0, smax, 5)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_linetype_manual(values = model_ltypes, breaks = present, labels = lab, name = NULL) +
    scale_shape_manual(values = model_shapes,    breaks = present, labels = lab, name = NULL) +
    scale_colour_manual(values = model_colour,   breaks = present, labels = lab, name = NULL) +
    scale_linewidth_manual(values = model_lwd,   breaks = present, labels = lab, name = NULL) +
    labs(x = expression("Sparsity " * italic(s)), y = NULL) +
    guides(linetype = guide_legend(nrow = 1)) +
    theme_bw(base_size = 9, base_family = "Helvetica") +
    theme(
      panel.grid.minor  = element_blank(),
      strip.text        = element_text(size = 8.5),
      legend.position   = "bottom",
      legend.key        = element_blank(),
      legend.key.width  = unit(2.6, "lines"),
      legend.text       = element_text(size = 8, margin = margin(l = 2, r = 5)),
      legend.margin     = margin(t = -2)
    )
}

# JASA text width is about 6.5 in. Each figure is written twice: a vector PDF
# (resolution-independent, preferred for submission) and a 600 dpi PNG.
save_fig <- function(fig, file, width = 6.5, height = 4.4, dpi = 600) {
  out <- file.path(res_dir, file)
  ggsave(out, fig, width = width, height = height, device = cairo_pdf)
  out_png <- sub("\\.pdf$", ".png", out)
  grDevices::png(out_png, width = width, height = height, units = "in", res = dpi,
                 type = "cairo")
  print(fig)
  grDevices::dev.off()
  cat("Saved ->", normalizePath(out), "and", basename(out_png), sprintf("(%d dpi)", dpi), "\n")
}

## ----- Gaussian: rows n = 75, n = 150 --------------------------------------
dat_gauss <- read_long(c("italic(n) == 75*','~~italic(p) == 100"  = "gaussian_n75_p100",
                         "italic(n) == 150*','~~italic(p) == 100" = "gaussian_n150_p100"))
save_fig(make_fig(dat_gauss), "simu1_gaussian.pdf")

## ----- GLM: rows binomial, poisson -----------------------------------------
dat_glm <- read_long(c("Logistic*','~italic(n) == 400" = "binomial_n400_p100",
                       "Poisson*','~italic(n) == 200"  = "poisson_n200_p100"))
save_fig(make_fig(dat_glm), "simu1_glm.pdf")
