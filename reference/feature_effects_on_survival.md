# Effect of one feature on the Cox survival curve.

Builds the family of survival curves obtained by varying a single
covariate while holding the others at their column means. The returned
object has the same shape as
[`predict_survival_function()`](https://vcmaxouuu.github.io/picreg/reference/predict_survival_function.md),
so it composes directly with
[`plot_survival_curves()`](https://vcmaxouuu.github.io/picreg/reference/plot_survival_curves.md).

## Usage

``` r
feature_effects_on_survival(object, idx, values = NULL)
```

## Arguments

- object:

  A fitted `pic.cox` object.

- idx:

  Index of the feature to vary. Either an integer column position in the
  training `X` or, if `X` carried column names, the variable name. The
  feature must lie in the model's selected support; otherwise the curve
  would be flat in `v` and the call is rejected.

- values:

  Optional numeric vector of values to evaluate. When `NULL` (default),
  the cached grid described in *Details* is used.

## Value

A list with components `time` (length `K`) and `survival` (matrix
`K x length(values)`, one column per evaluated value). Column names of
`survival` are formatted as `"<feature_name> = <value>"` and are picked
up automatically by
[`plot_survival_curves()`](https://vcmaxouuu.github.io/picreg/reference/plot_survival_curves.md)
for the legend.

## Details

Both the per-column mean row and a default grid of representative values
(used when `values = NULL`) are cached on the fit by
[`pic()`](https://vcmaxouuu.github.io/picreg/reference/pic.md) under
`attr(fit, "preproc")`, so the training design matrix does not need to
be passed back in. The cached default grid uses the unique values of the
column when there are at most five of them (handy for ordinal /
categorical covariates), and the four equispaced empirical quantiles
(`0`, `1/3`, `2/3`, `1`) otherwise.

## See also

[`predict_survival_function()`](https://vcmaxouuu.github.io/picreg/reference/predict_survival_function.md),
[`plot_survival_curves()`](https://vcmaxouuu.github.io/picreg/reference/plot_survival_curves.md).
