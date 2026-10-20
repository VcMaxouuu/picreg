# Linear predictor / response prediction for a pic fit.

Linear predictor / response prediction for a pic fit.

## Usage

``` r
# S3 method for class 'pic'
predict(object, newx, type = c("response", "link", "class", "survival"), ...)
```

## Arguments

- object:

  A fitted `pic` object.

- newx:

  Matrix of new values at which predictions are to be made.

- type:

  `"link"` (linear predictor), `"response"` (default; family `g`-link
  applied, and for Cox the relative risk \\\exp(\eta)\\, as in glmnet),
  `"class"` (binomial only; `0`/`1` labels, where `1` stands for `TRUE`
  or the second factor level of the training response), or `"survival"`
  (Cox only; returns a list with components `time` and `survival`).

- ...:

  Unused; present for S3 method consistency.

## Value

A numeric vector for types `"link"`, `"response"`, and `"class"`. For
`type = "survival"`, a list with `time` (length `K`) and `survival`
(matrix `K x m`, one column per subject).
