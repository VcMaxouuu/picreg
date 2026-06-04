# Plot Cox baseline cumulative hazard and baseline survival.

Two-panel step plot for a fitted `pic.cox` model: the Breslow baseline
cumulative hazard \\H_0(t)\\ on top and the baseline survival \\S_0(t) =
\exp(-H_0(t))\\ below.

## Usage

``` r
plot_baseline(model)
```

## Arguments

- model:

  A fitted `pic.cox` object.

## Value

Invisibly returns `NULL`.
