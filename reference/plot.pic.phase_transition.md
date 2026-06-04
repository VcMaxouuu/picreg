# Phase-transition plot for a `pic.phase_transition` object.

Plots the chosen recovery metric as a function of the sparsity level
`s`. Curves are distinguished by line type (grayscale ramp beyond five
curves). When multiple `(n, p)` configurations are compared across
several penalties, panels are laid out in a grid with at most three
columns per row.

## Usage

``` r
# S3 method for class 'pic.phase_transition'
plot(x, metric = c("exact_recovery", "tpr", "fdr"), ...)
```

## Arguments

- x:

  An object returned by
  [`phase_transition()`](https://vcmaxouuu.github.io/picreg/reference/phase_transition.md).

- metric:

  One of `"exact_recovery"` (default), `"tpr"`, `"fdr"`.

- ...:

  Additional graphical parameters forwarded to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns `x`.
