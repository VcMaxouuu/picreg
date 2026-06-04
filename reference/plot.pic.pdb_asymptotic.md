# Plot of the PDB asymptotic behavior.

Multi-panel histogram comparison of the simulated null gradient-norm
statistic under the `"mc_exact"` (light grey fill) and `"mc_gaussian"`
(dashed outline) selectors, one panel per `n` in `n_grid`. Two vertical
lines are added per panel:

- \\\hat\lambda^{\rm PDB}\_\alpha\\ - solid navy, the empirical (1 -
  alpha) quantile of `"mc_exact"` (what pic would actually use).

- \\\hat\lambda\_{analytical}\\ - dashed black, the Bonferroni
  closed-form bound.

## Usage

``` r
# S3 method for class 'pic.pdb_asymptotic'
plot(x, breaks = 40L, ...)
```

## Arguments

- x:

  An object returned by
  [`pdb_asymptotic()`](https://vcmaxouuu.github.io/picreg/reference/pdb_asymptotic.md).

- breaks:

  Number of histogram bins (default 40).

- ...:

  Unused; present for S3 method consistency.

## Value

Invisibly returns `x`.
