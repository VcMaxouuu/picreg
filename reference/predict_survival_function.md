# Survival curves for new data from a fitted Cox pic model.

Survival curves for new data from a fitted Cox pic model.

## Usage

``` r
predict_survival_function(object, newx)
```

## Arguments

- object:

  A fitted `pic.cox` model.

- newx:

  New design matrix (rows = subjects).

## Value

A list with `time` (length `K`) and `survival` (matrix `K x m`, one
column per subject).
