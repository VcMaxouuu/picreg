# Plot subject-specific Cox survival curves.

Step-line visualization of the output of
[`predict_survival_function()`](https://vcmaxouuu.github.io/picreg/reference/predict_survival_function.md)
or
[`feature_effects_on_survival()`](https://vcmaxouuu.github.io/picreg/reference/feature_effects_on_survival.md):
one survival curve per subject (or per feature value) on a common time
grid. To keep curves distinguishable, each curve is drawn with its own
line type and a sparse set of marker glyphs is overlaid at regular time
points.

## Usage

``` r
plot_survival_curves(
  sf,
  subjects = NULL,
  max_subjects = 10L,
  labels = NULL,
  n_marks = 8L,
  main = "Individual survival curves",
  ...
)
```

## Arguments

- sf:

  A list as returned by
  [`predict_survival_function()`](https://vcmaxouuu.github.io/picreg/reference/predict_survival_function.md)
  or
  [`feature_effects_on_survival()`](https://vcmaxouuu.github.io/picreg/reference/feature_effects_on_survival.md),
  with components `time` (length `K`) and `survival` (matrix `K x m`,
  one column per curve).

- subjects:

  Optional integer vector selecting which columns of `sf$survival` to
  plot. Defaults to all curves, capped to the first `max_subjects` for
  legibility.

- max_subjects:

  Maximum number of curves drawn when `subjects` is `NULL`. Default 10.

- labels:

  Optional character vector of labels used in the legend, one per
  plotted curve. Defaults to the column names of `sf$survival` when set,
  otherwise `"subject 1"`, `"subject 2"`, etc.

- n_marks:

  Number of marker glyphs overlaid on each curve to help distinguish
  them. Default `8`. Set to `0` to disable.

- main:

  Plot title.

- ...:

  Additional graphical parameters forwarded to
  [`graphics::matplot()`](https://rdrr.io/r/graphics/matplot.html).

## Value

Invisibly returns `NULL`.
