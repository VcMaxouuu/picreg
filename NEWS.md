# picreg 0.2.0

## Breaking changes

* The `relax` argument of `pic()` has been removed. To debias the
  coefficients, refit on `fit$selected` with the estimator of your choice
  (e.g. `lm()` or `glm()`).

* For Cox fits, `predict(type = "response")` now returns the relative risk
  `exp(eta)`, as in `glmnet`, instead of the linear predictor, and `coef()`
  no longer has an `(Intercept)` row.

* For SCAD and MCP, the solver now follows Section 3.3 of the paper: the
  l1 path runs down to the target `lambda`, and the SCAD / MCP stage starts
  from that l1 solution (it previously started from the l1 solution at the
  penultimate path value). Selected supports may differ slightly from
  previous versions.

* The default number of Monte Carlo draws is now 2000 everywhere, as in the
  paper: `lambda_pdb()`, `phase_transition()` and `pdb_asymptotic()` used
  5000.

* With `intercept = FALSE` and `standardize = TRUE`, the `(Intercept)` row
  of `coef()` now absorbs the centering of the columns on the original
  scale, so that `coef()` reproduces `predict()`. It is zero when `X` is
  centered.

* `fit$beta` is now a plain numeric vector (it was a one-column matrix).

## New features

* `pic()` accepts a `pic.lambda_pdb` object as `lambda`, typically
  `fit$lambda_pdb`, to reuse a calibration across penalties while keeping
  it on the new fit for `pdb_summary()` and `plot()`. Its family and design
  dimensions are checked.

* New `lambda_seed` argument of `pic()` and `seed` argument of
  `lambda_pdb()` for a reproducible PDB calibration. Only the calibration is
  seeded: the global random stream is restored afterwards. Default `NULL`.

* Binomial responses may be logical (`TRUE` mapped to 1) or a two-level
  factor (second level mapped to 1). `assess()` maps new responses with the
  training levels.

* Fits store the matched `call`, the total number of FISTA iterations
  `n_iter` and a `converged` flag, and `pic()` warns when the final FISTA
  step reaches `maxit`.

* `lambda_pdb()` results also store `seed`, `family`, `n` and `p`.

* `phase_transition()`:
  * failed fits are excluded from the averages, counted in a new `n_failed`
    column of `results` and reported in a single warning (they were silently
    counted as failed recoveries);
  * with `parallel = TRUE`, the caller's `future` plan is restored on exit
    (it was reset to `sequential`), and `workers` defaults to
    `future::availableCores(omit = 1)`, which honors the core limits set by
    CRAN or a scheduler.

## Bug fixes

* `feature_effects_on_survival()` failed on Cox fits with
  `standardize = FALSE`.

* A constant column of `X` triggered a spurious "does not appear
  standardized" warning. Constant columns are also detected robustly to
  rounding.

* Arguments of `pic()` and `lambda_pdb()` are validated with explicit
  messages. Previously, e.g., `lambda = c(0.1, 0.2)` gave a cryptic error,
  `lambda = -1` produced `NaN`s and `lambda_alpha = 0` an infinite `lambda`.

* `predict()` reports a mismatch between the number of columns of `newx`
  and of the training `X`.

* `print()` of a `phase_transition()` result listed no curve, and `print()`
  of a family descriptor missed its final newline.

* Gumbel family:
  * the null intercept and scale, used for the warm start, `lambda_max` and
    the `"mc_exact"` calibration, are now the joint location-scale
    maximum-likelihood estimates (previously a moment-based intercept
    followed by a conditional scale estimate);
  * the profiled scale is solved by a bracketed, safeguarded Newton method
    without clipping the exponentials, and the loss and gradient share a
    single solve;
  * fitting is numerically stable under large shifts of the response, and a
    constant response raises an explicit error;
  * `assess()` reports the deviance with the scale profiled at its
    maximum-likelihood estimate instead of a moment estimate.

## Performance

* The Monte Carlo PDB calibration processes the draws in blocks, which
  bounds its memory (it allocated two `n x n_simu` matrices, about 640 MB
  for `n = 20000`) without changing the result for a given seed. The Cox
  permutation calibration stores the simulated events as bytes.

* `concordance_index()` is computed in C++ in O(n log n) time and O(n)
  memory instead of O(n^2) matrices, and `cox_partial_log_likelihood()` no
  longer loops over the unique times.

* The Cox solver evaluates the loss and its gradient in a single pass over
  the risk sets and no longer records unused per-iteration diagnostics.

* The per-column value grids used by `feature_effects_on_survival()` are
  only computed for Cox fits (they were computed and discarded for every
  family, which dominated preprocessing for wide designs).
  `feature_effects_on_survival()` predicts all values in one call.

## Dependencies

* `future` and `future.apply` move from Imports to Suggests: they are only
  needed for `phase_transition(parallel = TRUE)`. `parallel` is no longer
  used.

## Documentation

* Fixed the documentation of `lambda_pdb()` (probability of selecting a
  variable under the null, tightness of the Bonferroni bound), of the
  `selected` and `beta` components of a fit, and of the value returned by
  `phase_transition()`, plus several typos in the vignette (formulas of the
  Gaussian and exponential objectives, `PBD`, section references).

* The vignette comparison with `glmnet` uses random signs for the non-zero
  coefficients, as stated in the text.

# picreg 0.1.4

## New features

* `predict()` for Cox fits now accepts `type = "survival"`, which returns
  subject-specific survival curves (a list with `time` and `survival`
  components). This replaces the standalone `predict_survival_function()`,
  which has been removed.

## Breaking changes

* `predict_survival_function()` has been removed. Use
  `predict(object, newx, type = "survival")` instead.

# picreg 0.1.3

## Bug fixes

* Resolved a C++ One Definition Rule (ODR) violation reported under
  link-time optimization (LTO): the file-local helper structs in
  `fista.cpp` and `fista_cox.cpp` now have internal linkage (anonymous
  namespaces).

## Minor improvements

* `coef()` now returns a one-column **sparse matrix** instead of a two-column data frame. 
  Zero coefficients are printed as `.`, which keeps the display compact in high dimensions.
  
* Added a `summary()` method for `pic` fits.

* `intercept = FALSE` is now permitted for every family. Previously it raised
  an error for non-Cox families when `standardize = TRUE`. When
  `intercept = FALSE`, the data are assumed to be centered.

* Standardized the input-validation error messages (consistent back-quoting
  of argument names and `call. = FALSE`).
