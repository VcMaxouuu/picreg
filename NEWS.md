# picreg 0.1.0

Initial release.

## Features

- `pic(X, y, family, penalty)` — top-level fit. Selects the regularisation
  parameter via the Pivotal Detection Boundary (PDB) and runs a 4-step
  warm-started FISTA path (three L1 warm-ups + one final step under the
  requested penalty).
- Six families dispatched by string: `"gaussian"`, `"binomial"`,
  `"poisson"`, `"exponential"`, `"gumbel"`, `"cox"`.
- Three penalties: `"lasso"`, `"scad"` (concavity `scad_a`), `"mcp"`
  (concavity `mcp_gamma`).
- `lambda_pdb()` — Pivotal Detection Boundary selector with three
  methods: `"mc_exact"`, `"mc_gaussian"`, `"analytical"`. Family-aware
  Monte Carlo for non-Cox families; dedicated permutation-based pivot
  for Cox.
- S3 methods on fitted objects: `predict()`, `coef()`, `print()`,
  `plot()`.
- Cox utilities: `concordance_index()`, `cox_partial_log_likelihood()`,
  `baseline_functions()`, `predict_survival_function()`, `plot_baseline()`.
- Empirical diagnostics: `phase_transition()` (support-recovery curves),
  `pdb_asymptotic()` (null-distribution convergence).
- C++ backend via Rcpp + RcppArmadillo for FISTA, family/penalty
  registries, and PDB Monte Carlo.
- Three vignettes: introduction, Cox regression, diagnostics.
