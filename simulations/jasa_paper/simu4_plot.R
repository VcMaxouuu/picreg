# =============================================================================
# Simulation 4 figure: ratio lambda / lambda_alpha of each calibration to the
# reference calibration, against n, one panel per family.
#   line at 1 : reference calibration at the reference nuisance value, for every n
#   boxes     : R repetitions of the Monte Carlo calibration at each mu_ref and
#               of the Gaussian approximation; the width is the Monte Carlo
#               error, the offset from one the error due to the unknown
#               nuisance parameter (MC) or to the Gaussian approximation
# Reads results/simu4_lambdas.rds written by simu4_calibration.R.
#
# Styling matches simu1_plot.R (publication-ready, JASA text width): greyscale,
# theme_bw with boxed strips on top, vector PDF via cairo_pdf plus a 600 dpi PNG.
# =============================================================================
library(ggplot2)

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
res_dir <- file.path(script_dir, "results")

lam <- readRDS(file.path(res_dir, "simu4_lambdas.rds"))
lam$n_f <- factor(lam$n, levels = sort(unique(lam$n)))

## ----- panels: families present in the results -----------------------------
fam_info <- data.frame(
    family = c("gaussian", "binomial", "poisson"),
    name   = c("Gaussian", "Logistic", "Poisson"),
    nuis   = c("sd", "mean", "mean")
)
fam_info <- fam_info[fam_info$family %in% lam$family, ]
fams <- fam_info$family
mu0 <- tapply(lam$mu0, lam$family, unique)
panel_labels <- sprintf("%s (reference %s %g)", fam_info$name, fam_info$nuis, mu0[fams])
lam <- lam[lam$family %in% fams, ]
lam$panel <- factor(lam$family, levels = fams, labels = panel_labels)

## ----- calibration order and labels (references sorted by value) -----------
cal_of <- function(fam) {
    refs <- unique(lam$calibration[lam$family == fam & grepl("^mc_", lam$calibration)])
    refs[order(as.numeric(sub("mc_", "", refs)))]
}
boxes <- lam[grepl("^mc_", lam$calibration) | lam$calibration == "gaussian", ]
boxes$slot <- ave(seq_len(nrow(boxes)), boxes$family, FUN = function(i) {
    k <- boxes$calibration[i]
    ifelse(k == "gaussian", 4L, match(k, cal_of(boxes$family[i[1]])))
})
boxes$slot <- factor(boxes$slot,
    levels = 1:4,
    labels = c(
        "MC, low value", "MC, middle value",
        "MC, high value", "Gaussian approximation"
    )
)
reference <- lam[lam$calibration == "reference", c("family", "n", "lambda")]
boxes$ratio <- boxes$lambda /
    reference$lambda[match(paste(boxes$family, boxes$n), paste(reference$family, reference$n))]

# Reference values used for the MC calibrations, written in each panel.
ref_note <- do.call(rbind, lapply(seq_along(fams), function(k) {
    data.frame(
        panel = factor(panel_labels[k], levels = panel_labels),
        label = paste0(
            "other ", fam_info$nuis[k], "s: ",
            paste(sub("mc_", "", cal_of(fams[k])), collapse = ", ")
        )
    )
}))

## ----- figure ---------------------------------------------------------------
fig <- ggplot(boxes, aes(n_f, ratio)) +
    geom_hline(yintercept = 1, linetype = "22", linewidth = 0.5, colour = "black") +
    geom_boxplot(aes(fill = slot),
        position = position_dodge(width = 0.8), width = 0.7,
        linewidth = 0.3, colour = "black",
        outlier.size = 0.4, outlier.shape = 16, outlier.colour = "grey40"
    ) +
    geom_text(
        data = ref_note, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
        hjust = 1.04, vjust = 1.6, size = 2.6, colour = "grey20", family = "Helvetica"
    ) +
    facet_wrap(~panel) +
    scale_fill_manual(values = c("grey30", "grey60", "grey85", "white"), name = NULL) +
    scale_y_continuous(breaks = scales::breaks_width(0.05)) +
    labs(
        x = expression("Sample size " * italic(n)),
        y = expression(hat(lambda)[alpha]^{PDB} / lambda[alpha])
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    theme_bw(base_size = 9, base_family = "Helvetica") +
    theme(
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text         = element_text(size = 8.5),
        legend.position    = "bottom",
        legend.key         = element_blank(),
        legend.key.width   = unit(1.2, "lines"),
        legend.text        = element_text(size = 8, margin = margin(l = 2, r = 8)),
        legend.margin      = margin(t = -2)
    )

# JASA text width is about 6.5 in. Written twice: a vector PDF
# (resolution-independent, preferred for submission) and a 600 dpi PNG.
# Greek letters of plotmath (lambda, alpha) in a serif symbol font, so that
# alpha is not mistaken for an 'a' (Helvetica's alpha looks like one).
fig_w <- 6.5
fig_h <- 3.6
dpi <- 600
symbol_font <- grDevices::cairoSymbolFont("Times New Roman")
out <- file.path(res_dir, "simu4_calibration.pdf")
pdf_device <- function(...) grDevices::cairo_pdf(..., symbolfamily = symbol_font)
ggsave(out, fig, width = fig_w, height = fig_h, device = pdf_device)
