# Validate and optionally standardize the design matrix.

Standardizes columns to zero mean and unit variance (population sd; `n`
divisor) when `standardize_X` is `TRUE`. When `X_mean` and `X_std` are
supplied they are reapplied without recomputing — used at prediction
time on new data. Column names of `X` (or column names of the source
data frame) are captured and returned as `feature_names`; they are
propagated through
[`pic()`](https://vcmaxouuu.github.io/picreg/reference/pic.md) to
annotate `fit$selected`, [`coef()`](https://rdrr.io/r/stats/coef.html),
and the coefficient plot.

## Usage

``` r
check_X(X, standardize_X = TRUE, X_mean = NULL, X_std = NULL)
```

## Arguments

- X:

  A numeric matrix or coercible (data.frame).

- standardize_X:

  Logical; standardize columns of `X`.

- X_mean:

  Optional pre-computed column means.

- X_std:

  Optional pre-computed column standard deviations.

## Value

A list with components `X`, `X_mean`, `X_std`, `feature_names` (the
column names of `X` if any, otherwise `NULL`).
