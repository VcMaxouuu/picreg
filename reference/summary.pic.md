# Summarize a fitted pic model.

Returns a structured summary of the fit: the family, penalty, selected
\\\hat\lambda\\, the problem dimensions, the number of selected
variables, the intercept, and a table of the **non-zero** coefficients
ordered by decreasing absolute value. Coefficients are returned on the
original scale of `X` by default (pass `standardized = TRUE` for the
internal standardized scale, matching `fit$beta`).

## Usage

``` r
# S3 method for class 'pic'
summary(object, standardized = FALSE, ...)
```

## Arguments

- object:

  A fitted `pic` object.

- standardized:

  Logical; if `TRUE`, summarize the standardized coefficients used
  internally. Default `FALSE` (original scale of `X`).

- ...:

  Unused; present for S3 method consistency.

## Value

An object of class `"summary.pic"`: a list with elements `family`,
`penalty`, `lambda`, `n_samples`, `n_features`, `df`, `intercept`,
`standardized`, and `coefficients` (a two-column `data.frame` of the
non-zero coefficients).

## Examples

``` r
data(QuickStartExample)
fit <- pic(QuickStartExample$X, QuickStartExample$y)
summary(fit)
#> pic fit summary
#>   family    : gaussian
#>   penalty   : lasso
#>   lambda    : 0.3099
#>   dimensions: n = 100, p = 30
#>   selected  : 5 / 30
#>   intercept : 0.006378
#> 
#>   Non-zero coefficients (original scale):
#>  variable coefficient
#>    gene_3     -0.7544
#>    gene_4     -0.6537
#>    gene_1     -0.4065
#>    gene_5      0.1422
#>    gene_2     -0.0502
```
