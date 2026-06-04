# Horizontal lollipop plot of the non-zero coefficients of a pic fit.

One row per selected variable, sorted by descending absolute coefficient
value (largest at the top). Each variable is drawn as a horizontal
segment from zero to its fitted value.

## Usage

``` r
# S3 method for class 'pic'
plot(x, standardized = TRUE, max_features = NULL, ...)
```

## Arguments

- x:

  A fitted `pic` object.

- standardized:

  Logical; if `TRUE`, plot the standardized coefficients used internally
  during fitting. If `FALSE`, use the coefficients on the original scale
  of `X`. Default `TRUE`.

- max_features:

  Optional cap on the number of features displayed (the strongest are
  kept).

- ...:

  Additional graphical parameters forwarded to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html) for
  the empty frame.

## Value

Invisibly returns the plotted (named) coefficient vector.
