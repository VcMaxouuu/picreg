# Linear predictor / response prediction for a pic fit.

Linear predictor / response prediction for a pic fit.

## Usage

``` r
# S3 method for class 'pic'
predict(object, newx, type = c("response", "link", "class"), ...)
```

## Arguments

- object:

  A fitted `pic` object.

- newx:

  Matrix of new values at which predictions are to be made.

- type:

  `"link"` (linear predictor) or `"response"` (default; family `g`-link
  applied).

- ...:

  Unused; present for S3 method consistency.

## Value

A numeric vector.
