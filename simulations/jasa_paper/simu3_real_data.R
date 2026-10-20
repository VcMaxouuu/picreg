# =============================================================================
# Simulation 3 -- real data: predictive performance vs model size.
#
# Seven datasets (4 regression, 3 binary classification). For each dataset,
# n_splits random 70/30 train/test splits. On each split:
#   - covariates are standardized with the TRAIN mean / population sd
#     (columns constant on the train set are dropped),
#   - every selector picks a support on the train set,
#   - an UNPENALIZED model (lm for regression, logistic glm for
#     classification) is refitted on the train set restricted to that support,
#   - its test performance is recorded: MSE (regression) or accuracy at the
#     0.5 threshold (classification), together with the model size |S|.
#     For regression the relative MSE is also recorded:
#       rel_mse = MSE_test / mean((y_test - mean(y_train))^2),
#     i.e. the MSE relative to the intercept-only model fitted on the train
#     set (1 = no better than the mean; unit-free, comparable across datasets).
#
# Models: PIC-l1, PIC-scad, PIC-l0 (regression only), cv.glmnet, BIC, EBIC.
# Selectors, forward selection and PIC detection boundaries are the ones of
# common.R, recomputed on each training set.
#
# Output (results/):
#   simu3_<dataset>.csv       one row per (split, model): size, perf, rel_mse
#   simu3_summary.csv         mean and 1 SE (sd / sqrt(n_splits)) per cell
#   simu3_table.csv           formatted "mean ± se" table (perf, rel_mse, size)
# =============================================================================

## locate this folder (Rscript, source(), or RStudio)
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
source(file.path(script_dir, "common.R"))

## ------------------------------ parameters ---------------------------------
n_splits   <- 100               # random train/test splits per dataset
prop_train <- 0.70              # train fraction
data_dir   <- file.path(script_dir, "data")
dir.create(data_dir, showWarnings = FALSE)

## ------------------------------ loaders ------------------------------------
# Each loader returns list(X = numeric matrix, y = numeric vector, task).
# task is "regression" (gaussian) or "classification" (binomial, y in {0, 1}).

# Mean imputation of numeric columns, dummy coding of the others.
clean_features <- function(df) {
  for (j in seq_along(df)) {
    if (is.numeric(df[[j]])) {
      df[[j]][is.na(df[[j]])] <- mean(df[[j]], na.rm = TRUE)
    } else {
      v <- as.character(df[[j]])
      v[is.na(v)] <- names(which.max(table(v)))      # most frequent category
      df[[j]] <- factor(v)
    }
  }
  X <- model.matrix(~ ., data = df)
  X[, colnames(X) != "(Intercept)", drop = FALSE]
}

to_binary <- function(y) {
  if (is.numeric(y) && all(y %in% c(0, 1))) return(as.numeric(y))
  as.numeric(factor(y)) - 1
}

# Generic CSV loader. target = column name, or NULL for the last column.
# For classification, `positive` is the label coded 1 (others 0).
# Index columns written by pandas / write.csv ("X", "Unnamed...") are dropped.
load_csv <- function(file, target = NULL, drop = character(0), task,
                     positive = NULL) {
  path <- file.path(data_dir, file)
  if (!file.exists(path)) stop("Missing data file: ", path, call. = FALSE)
  d <- read.csv(path, check.names = TRUE)
  idx_cols <- names(d) == "X" | startsWith(names(d), "Unnamed")
  d <- d[, !idx_cols & !(names(d) %in% drop), drop = FALSE]
  if (is.null(target)) target <- names(d)[ncol(d)]
  y <- d[[target]]
  ok <- !is.na(y)
  X <- clean_features(d[ok, names(d) != target, drop = FALSE])
  y <- if (task == "regression") as.numeric(y[ok])
  else if (!is.null(positive)) as.numeric(y[ok] == positive)
  else to_binary(y[ok])
  list(X = X, y = y, task = task)
}

# Download a UCI file once and cache it in data/.
uci_file <- function(url, file) {
  path <- file.path(data_dir, file)
  if (!file.exists(path)) download.file(url, path, quiet = TRUE)
  path
}
uci <- "https://archive.ics.uci.edu/ml/machine-learning-databases"

# Communities & Crime (UCI 183): drop the 5 non-predictive columns, keep
# features with < 10% missing, mean-impute. Target: ViolentCrimesPerPop.
load_communities <- function() {
  d <- read.csv(uci_file(file.path(uci, "communities/communities.data"),
                         "communities.data"),
                header = FALSE, na.strings = "?")
  y <- d[[ncol(d)]]
  X <- d[, 6:(ncol(d) - 1)]          # cols 1-5: state, county, community, name, fold
  X <- X[, colMeans(is.na(X)) < 0.1]
  list(X = clean_features(X), y = y, task = "regression")
}

# Ionosphere (UCI 52): drop Attribute2 (constant), class g = 1, b = 0.
load_ionosphere <- function() {
  d <- read.csv(uci_file(file.path(uci, "ionosphere/ionosphere.data"),
                         "ionosphere.data"), header = FALSE)
  y <- as.numeric(d[[ncol(d)]] == "g")
  X <- as.matrix(d[, setdiff(seq_len(ncol(d) - 1), 2)])
  list(X = X, y = y, task = "classification")
}

# Breast Cancer Wisconsin Diagnostic (= sklearn load_breast_cancer):
# id, diagnosis, 30 features. sklearn coding: malignant = 0, benign = 1.
load_breast_cancer <- function() {
  d <- read.csv(uci_file(file.path(uci, "breast-cancer-wisconsin/wdbc.data"),
                         "wdbc.data"), header = FALSE)
  y <- as.numeric(d[[2]] == "B")
  X <- as.matrix(d[, 3:ncol(d)])
  list(X = X, y = y, task = "classification")
}

# Diabetes (Efron et al. 2004, package 'lars'): x2 = 442 x 64 design with the
# 10 main effects, 45 two-way interactions and 9 squares (sex is binary, so it
# has no square). Target: disease progression after one year.
load_diabetes <- function() {
  e <- new.env()
  data("diabetes", package = "lars", envir = e)
  list(X = unclass(e$diabetes$x2), y = as.numeric(e$diabetes$y),
       task = "regression")
}

# CSV datasets (in data/), same preprocessing as the original Python code:
#   prostate   : ProstateCancer.csv, drop "train", target "lpsa"
#   riboflavin : riboflavin.csv, target "y"
#   sonar      : sonar.csv, target "Class", R (rock) = 1, M (mine) = 0
datasets <- list(
  prostate      = function() load_csv("ProstateCancer.csv", target = "lpsa",
                                      drop = "train", task = "regression"),
  diabetes      = load_diabetes,
  communities   = load_communities,
  riboflavin    = function() load_csv("riboflavin.csv", target = "y",
                                      task = "regression"),
  breast_cancer = load_breast_cancer,
  ionosphere    = load_ionosphere,
  sonar         = function() load_csv("sonar.csv", target = "Class",
                                      task = "classification", positive = "R")
)

## ------------------------- one train/test split ----------------------------
# Standardize with train statistics; drop columns constant on the train set.
standardize_split <- function(Xtr, Xte) {
  mu <- colMeans(Xtr)
  sd <- sqrt(colMeans(sweep(Xtr, 2, mu, "-")^2))
  keep <- sd > 1e-10
  list(tr = sweep(sweep(Xtr[, keep, drop = FALSE], 2, mu[keep], "-"), 2, sd[keep], "/"),
       te = sweep(sweep(Xte[, keep, drop = FALSE], 2, mu[keep], "-"), 2, sd[keep], "/"))
}

# Unpenalized refit on the selected support, evaluated on the test set.
# Rank-deficient refits (|S| >= n_train) get NA coefficients set to 0.
# Returns c(perf, rel_mse): perf = MSE (regression) or accuracy
# (classification); rel_mse = MSE / MSE of the train-mean predictor
# (regression only, NA for classification).
refit_perf <- function(Xtr, ytr, Xte, yte, S, family) {
  Ztr <- cbind(1, Xtr[, S, drop = FALSE])
  Zte <- cbind(1, Xte[, S, drop = FALSE])
  if (family == "gaussian") {
    b <- lm.fit(Ztr, ytr)$coefficients
    b[is.na(b)] <- 0
    mse <- mean((yte - Zte %*% b)^2)
    c(perf = mse, rel_mse = mse / mean((yte - mean(ytr))^2))
  } else {
    b <- suppressWarnings(glm.fit(Ztr, ytr, family = binomial())$coefficients)
    b[is.na(b)] <- 0
    c(perf = mean((Zte %*% b > 0) == (yte == 1)),  # P(y = 1) > 0.5  <=>  eta > 0
      rel_mse = NA_real_)
  }
}

run_split <- function(ds, family) {
  n   <- length(ds$y)
  tr  <- sort(sample.int(n, round(prop_train * n)))
  te  <- setdiff(seq_len(n), tr)
  Z   <- standardize_split(ds$X[tr, , drop = FALSE], ds$X[te, , drop = FALSE])
  ytr <- ds$y[tr]; yte <- ds$y[te]
  
  lam    <- lambda_pdb(Z$tr, family = family, method = "mc_exact",
                       n_simu = n_simu_pdb, alpha = alpha)$value
  lam_l0 <- if (family == "gaussian") lambda_pic_l0(Z$tr, n_simu_pdb, alpha)
  fw     <- fit_forward(Z$tr, ytr, family, lam_l0)
  
  sel <- list(
    pic_l1   = fit_pic(Z$tr, ytr, family, "lasso", lam),
    pic_scad = fit_pic(Z$tr, ytr, family, "scad",  lam),
    pic_l0   = fw$pic_l0,
    glmnet   = fit_glmnet(Z$tr, ytr, family),
    bic      = fw$bic,
    ebic     = fw$ebic
  )
  sel <- sel[!vapply(sel, is.null, logical(1))]
  pf  <- vapply(sel, function(S) refit_perf(Z$tr, ytr, Z$te, yte, S, family),
                numeric(2))
  data.frame(
    model   = names(sel),
    size    = vapply(sel, length, numeric(1)),
    perf    = pf["perf", ],
    rel_mse = pf["rel_mse", ],
    row.names = NULL)
}

## ------------------------- run one dataset ---------------------------------
# Every split gets its own seed, so results do not depend on the number of cores.
run_dataset <- function(name, seed) {
  ds     <- datasets[[name]]()
  family <- if (ds$task == "regression") "gaussian" else "binomial"
  cat(sprintf("[%s] n=%d p=%d task=%s\n", name, nrow(ds$X), ncol(ds$X), ds$task))
  
  t0   <- Sys.time()
  rows <- mclapply(seq_len(n_splits), function(m) {
    set.seed(seed * 1e4 + m)
    cbind(split = m, run_split(ds, family))
  }, mc.cores = ncores)
  failed <- vapply(rows, inherits, logical(1), "try-error")
  if (any(failed)) stop(sprintf("[%s] %d splits failed:\n%s",
                                name, sum(failed), rows[[which(failed)[1]]]))
  res <- cbind(dataset = name, task = ds$task, n = nrow(ds$X), p = ncol(ds$X),
               do.call(rbind, rows))
  
  out_csv <- file.path(script_dir, sprintf("results/simu3_%s.csv", name))
  write.csv(res, out_csv, row.names = FALSE)
  cat(sprintf("[%s] done in %.1f min -> %s\n\n", name,
              as.numeric(difftime(Sys.time(), t0, units = "mins")),
              normalizePath(out_csv)))
  invisible(res)
}

## ------------------------------ summary ------------------------------------
se <- function(x) sd(x) / sqrt(length(x))

summarise_results <- function(all_res) {
  key <- unique(all_res[c("dataset", "task", "n", "p", "model")])
  summ <- do.call(rbind, lapply(seq_len(nrow(key)), function(i) {
    r <- merge(key[i, ], all_res)
    cbind(key[i, ],
          perf_mean = mean(r$perf), perf_se = se(r$perf),
          rel_mse_mean = mean(r$rel_mse), rel_mse_se = se(r$rel_mse),
          size_mean = mean(r$size), size_se = se(r$size))
  }))
  rownames(summ) <- NULL
  summ
}

# Wide table: one row per dataset, "mean ± se" for perf, relative MSE
# (regression only) and size per model.
format_table <- function(summ) {
  models <- c("pic_scad", "pic_l1", "pic_l0", "glmnet", "bic", "ebic")
  labels <- c("PIC-scad", "PIC-l1", "PIC-l0", "cv.glmnet", "BIC", "EBIC")
  dsets  <- unique(summ$dataset)
  tab <- do.call(rbind, lapply(dsets, function(d) {
    s   <- summ[summ$dataset == d, ]
    row <- data.frame(Dataset = d, Task = s$task[1], n = s$n[1], p = s$p[1],
                      Metric = if (s$task[1] == "regression") "test MSE" else "test accuracy")
    for (k in seq_along(models)) {
      r <- s[s$model == models[k], ]
      row[[paste(labels[k], "perf")]] <- if (nrow(r)) sprintf("%.3f ± %.3f", r$perf_mean, r$perf_se) else "–"
      row[[paste(labels[k], "rel MSE")]] <- if (nrow(r) && !is.na(r$rel_mse_mean))
        sprintf("%.3f ± %.3f", r$rel_mse_mean, r$rel_mse_se) else "–"
      row[[paste(labels[k], "size")]] <- if (nrow(r)) sprintf("%.1f ± %.1f", r$size_mean, r$size_se) else "–"
    }
    row
  }))
  tab
}

## ------------------------------ run ----------------------------------------
all_res <- do.call(rbind, lapply(seq_along(datasets), function(k)
  run_dataset(names(datasets)[k], seed = k)))

summ <- summarise_results(all_res)
write.csv(summ, file.path(script_dir, "results/simu3_summary.csv"), row.names = FALSE)

tab <- format_table(summ)
out_tab <- file.path(script_dir, "results/simu3_table.csv")
write.csv(tab, out_tab, row.names = FALSE, fileEncoding = "UTF-8")
print(tab, row.names = FALSE)
cat("\nSaved ->", normalizePath(out_tab), "\n")
