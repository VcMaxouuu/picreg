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
  applied), `"class"` (binomial only), or `"survival"` (Cox only;
  returns a list with components `time` and `survival`).

- ...:

  Unused; present for S3 method consistency.

## Value

A numeric vector for types `"link"`, `"response"`, and `"class"`. For
`type = "survival"`, a list with `time` (length `K`) and `survival`
(matrix `K x m`, one column per subject).
