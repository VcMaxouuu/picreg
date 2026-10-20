# Coefficients of a fitted pic model.

Returns the fitted coefficients as a one-column **sparse matrix** of
class `"dgCMatrix"` (from the Matrix package), mirroring the output of
[`coef()`](https://rdrr.io/r/stats/coef.html) for glmnet. The first row
is the intercept, labeled `"(Intercept)"`; the remaining rows are the
predictors. Cox fits have no intercept row, the Cox model being
invariant to a shift of the linear predictor. Row names are taken from
the column names of `X` when available (matrix `colnames(X)` or
data-frame column names) and otherwise default to `V1, ..., Vp`. Zero
coefficients are stored implicitly and printed as `"."`, which keeps the
display compact in high dimensions.

## Usage

``` r
# S3 method for class 'pic'
coef(object, standardized = FALSE, ...)
```

## Arguments

- object:

  A fitted `pic` object.

- standardized:

  Logical; if `TRUE`, return the coefficients on the standardized scale
  used internally during fitting. Default `FALSE` (return coefficients
  on the original scale of `X`).

- ...:

  Unused; present for S3 method consistency.

## Value

A sparse `(p + 1)` by `1` matrix of class `"dgCMatrix"`, with row names
`c("(Intercept)", <variables>)` and column name `"coefficient"` (`p` by
`1`, without the intercept row, for Cox fits). Use
[`as.numeric()`](https://rdrr.io/r/base/numeric.html) to obtain a plain
numeric vector, or `which(coef(fit) != 0)` to list the selected entries.

## Details

Internally the model is fitted on a standardized design matrix, so the
raw coefficients live on the standardized scale. By default this method
rescales them back to the **original scale** of `X` — the values to plug
into the un-standardized design for prediction — via \$\$beta\\orig =
beta / s \quad intercept\\orig = intercept - sum(m \* beta\\orig)\$\$
where `m` and `s` are the column mean and standard deviation. The
intercept absorbs the centering even for fits with `intercept = FALSE`,
so that the coefficients reproduce
[`predict.pic()`](https://vcmaxouuu.github.io/picreg/reference/predict.pic.md).
Pass `standardized = TRUE` to skip the rescaling and return the raw fit
values (identical to `fit$beta`, with intercept `0` when none was
fitted).
