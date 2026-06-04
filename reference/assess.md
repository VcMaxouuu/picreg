# Assess performance of a `pic` fit.

Given a test set, reports a small set of family-appropriate predictive
metrics. Optionally appends support-recovery diagnostics when the true
active set is known.

## Usage

``` r
assess(object, ...)

# S3 method for class 'pic'
assess(object, newx, newy, true_features = NULL, ...)
```

## Arguments

- object:

  A fitted `pic` object.

- ...:

  Unused; present for S3 method consistency.

- newx:

  Numeric design matrix at which predictions are evaluated, with the
  same columns as the training data.

- newy:

  Response on the new observations. Numeric vector for all families
  except Cox; for Cox, a two-column matrix `(time, event)` matching the
  training response format.

- true_features:

  Optional integer or character vector listing the indices (or names) of
  the true active variables. When supplied, support-recovery metrics are
  appended.

## Value

A two-column `data.frame` with columns `metric` (character) and `value`
(numeric).

## Details

The metrics depend on the family:

- `"gaussian"`:

  MSE, MAE, R-squared.

- `"binomial"`:

  accuracy, AUC, binomial deviance.

- `"poisson"`:

  MSE, MAE, Poisson deviance.

- `"exponential"`:

  MSE, MAE, Exponential deviance.

- `"gumbel"`:

  MSE, MAE, deviance computed from the per-sample negative
  log-likelihood with a moment estimate of the scale parameter
  \\\hat\sigma = \mathrm{sd}(y - \hat\eta)\sqrt{6}/\pi\\.

- `"cox"`:

  Harrell's C-index and the Breslow partial log-likelihood (negative,
  normalized by `n`).

When `true_features` is non-`NULL`, four support-recovery metrics are
appended (independent of `newx` / `newy`): `exact_recovery`, `tpr`
(true-positive rate / sensitivity), `fdr` (false-discovery rate) and
`f1` (harmonic mean of precision and recall). Names are accepted in
addition to integer positions when the fit carries column names.

## Examples

``` r
data(QuickStartExample)
X <- QuickStartExample$X
y <- QuickStartExample$y
fit <- pic(X, y, family = "gaussian", penalty = "scad")
assess(fit, X, y, true_features = paste0("gene_", 1:5))
#>          metric     value
#>             MSE 1.0551556
#>             MAE 0.8266425
#>              R2 0.7627548
#>  exact_recovery 1.0000000
#>             tpr 1.0000000
#>             fdr 0.0000000
#>              f1 1.0000000
```
