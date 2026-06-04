# Plot methods for pic fits and diagnostics.

Entry points:

## Details

- `plot(fit)` - lollipop plot of the non-zero coefficients.

- `plot(fit$lambda_pdb)` - histogram of the PDB null distribution with a
  vertical line at the selected \\\hat\lambda\\.

- [`plot_baseline()`](https://vcmaxouuu.github.io/picreg/reference/plot_baseline.md) -
  Cox-only: baseline cumulative hazard and baseline survival, two
  panels.

- `plot(phase_transition(...))` - recovery curves.

- `plot(pdb_asymptotic(...))` - null-distribution convergence panels.
