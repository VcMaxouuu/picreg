# =============================================================================
# Real-data results table: median #selected and median C-index per model
# =============================================================================
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

## ----- model display names (same as phase_transition_plot.R) ----------------
model_levels <- c("pic_lasso", "pic_scad", "glmnet", "bic", "bh")
model_labels <- c("PIC-lasso", "PIC-scad", "cv.glmnet", "BIC forward", "Cox + BH")
names(model_labels) <- model_levels

## ----- dataset display names ------------------------------------------------
dataset_files <- c(lung = "results_lung.csv", pbc = "results_pbc.csv",
                   gbsg = "results_gbsg.csv", lymphoma = "results_lymphoma.csv",
                   nki70 = "results_nki70.csv", mcl = "results_mcl.csv")

## ----- build table ----------------------------------------------------------
rows <- lapply(names(dataset_files), function(ds) {
  d <- read.csv(file.path(script_dir, "results", dataset_files[ds]))
  row <- data.frame(Dataset = ds, stringsAsFactors = FALSE)
  for (m in model_levels) {
    col_nsel <- paste0(m, "_nselected")
    col_ci   <- paste0(m, "_cindex")
    med_nsel <- median(d[[col_nsel]])
    med_ci   <- median(d[[col_ci]])
    row[[paste0(model_labels[m], " (nsel)")]]   <- med_nsel
    row[[paste0(model_labels[m], " (C-index)")]] <- round(med_ci, 2)
  }
  row
})

tab <- do.call(rbind, rows)
print(tab, row.names = FALSE)

## ----- save to CSV ----------------------------------------------------------
out_path <- file.path(script_dir, "results/real_data_table.csv")
write.csv(tab, out_path, row.names = FALSE)
cat("\nSaved ->", normalizePath(out_path), "\n")
