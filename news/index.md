# Changelog

## picreg 0.1.4

### New features

- [`predict()`](https://rdrr.io/r/stats/predict.html) for Cox fits now
  accepts `type = "survival"`, which returns subject-specific survival
  curves (a list with `time` and `survival` components). This replaces
  the standalone `predict_survival_function()`, which has been removed.

### Breaking changes

- `predict_survival_function()` has been removed. Use
  `predict(object, newx, type = "survival")` instead.

## picreg 0.1.3

CRAN release: 2026-06-04

### Bug fixes

- Resolved a C++ One Definition Rule (ODR) violation reported under
  link-time optimization (LTO): the file-local helper structs in
  `fista.cpp` and `fista_cox.cpp` now have internal linkage (anonymous
  namespaces).

### Minor improvements

- [`coef()`](https://rdrr.io/r/stats/coef.html) now returns a one-column
  **sparse matrix** instead of a two-column data frame. Zero
  coefficients are printed as `.`, which keeps the display compact in
  high dimensions.

- Added a [`summary()`](https://rdrr.io/r/base/summary.html) method for
  `pic` fits.

- `intercept = FALSE` is now permitted for every family. Previously it
  raised an error for non-Cox families when `standardize = TRUE`. When
  `intercept = FALSE`, the data are assumed to be centered.

- Standardized the input-validation error messages (consistent
  back-quoting of argument names and `call. = FALSE`).
