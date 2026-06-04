# Package index

## Model fitting

The single high-level entry point.

- [`pic()`](https://vcmaxouuu.github.io/picreg/reference/pic.md) :
  Sparse linear regression using the Pivotal Information Criterion.

## Methods for fitted models

S3 methods shared across all families.

- [`summary(`*`<pic>`*`)`](https://vcmaxouuu.github.io/picreg/reference/summary.pic.md)
  : Summarize a fitted pic model.

- [`coef(`*`<pic>`*`)`](https://vcmaxouuu.github.io/picreg/reference/coef.pic.md)
  : Coefficients of a fitted pic model.

- [`predict(`*`<pic>`*`)`](https://vcmaxouuu.github.io/picreg/reference/predict.pic.md)
  : Linear predictor / response prediction for a pic fit.

- [`plot(`*`<pic>`*`)`](https://vcmaxouuu.github.io/picreg/reference/plot.pic.md)
  : Horizontal lollipop plot of the non-zero coefficients of a pic fit.

- [`assess()`](https://vcmaxouuu.github.io/picreg/reference/assess.md) :

  Assess performance of a `pic` fit.

## Families and penalties

- [`pic_families`](https://vcmaxouuu.github.io/picreg/reference/pic_families.md)
  : GLM families for pic — descriptor layer.
- [`pic_penalties`](https://vcmaxouuu.github.io/picreg/reference/pic_penalties.md)
  : Sparsity-inducing penalties for pic.

## Choosing the penalty level

The Pivotal Detection Boundary selector and its diagnostics.

- [`lambda_pdb()`](https://vcmaxouuu.github.io/picreg/reference/lambda_pdb.md)
  : Pivotal Detection Boundary regularization selector
- [`pdb_summary()`](https://vcmaxouuu.github.io/picreg/reference/pdb_summary.md)
  : Summary of the PDB lambda selector.
- [`pdb_asymptotic()`](https://vcmaxouuu.github.io/picreg/reference/pdb_asymptotic.md)
  : Asymptotic behavior of the PDB null distribution.
- [`plot(`*`<pic.lambda_pdb>`*`)`](https://vcmaxouuu.github.io/picreg/reference/plot.pic.lambda_pdb.md)
  : Plot the PDB null distribution.
- [`plot(`*`<pic.pdb_asymptotic>`*`)`](https://vcmaxouuu.github.io/picreg/reference/plot.pic.pdb_asymptotic.md)
  : Plot of the PDB asymptotic behavior.

## Variable-selection diagnostics

- [`phase_transition()`](https://vcmaxouuu.github.io/picreg/reference/phase_transition.md)
  : Phase-transition analysis of support recovery.

- [`plot(`*`<pic.phase_transition>`*`)`](https://vcmaxouuu.github.io/picreg/reference/plot.pic.phase_transition.md)
  :

  Phase-transition plot for a `pic.phase_transition` object.

## Survival utilities

Helpers for the Cox proportional-hazards family.

- [`predict_survival_function()`](https://vcmaxouuu.github.io/picreg/reference/predict_survival_function.md)
  : Survival curves for new data from a fitted Cox pic model.

- [`feature_effects_on_survival()`](https://vcmaxouuu.github.io/picreg/reference/feature_effects_on_survival.md)
  : Effect of one feature on the Cox survival curve.

- [`plot_baseline()`](https://vcmaxouuu.github.io/picreg/reference/plot_baseline.md)
  : Plot Cox baseline cumulative hazard and baseline survival.

- [`plot_survival_curves()`](https://vcmaxouuu.github.io/picreg/reference/plot_survival_curves.md)
  : Plot subject-specific Cox survival curves.

- [`baseline_functions()`](https://vcmaxouuu.github.io/picreg/reference/baseline_functions.md)
  : Breslow baseline cumulative hazard and survival.

- [`concordance_index()`](https://vcmaxouuu.github.io/picreg/reference/concordance_index.md)
  : Harrell's concordance index.

- [`cox_partial_log_likelihood()`](https://vcmaxouuu.github.io/picreg/reference/cox_partial_log_likelihood.md)
  :

  Breslow partial log-likelihood (negative, normalized by `n`).

## Datasets

- [`QuickStartExample`](https://vcmaxouuu.github.io/picreg/reference/QuickStartExample.md)
  : Small Gaussian dataset for the introductory vignette.
- [`BinomialExample`](https://vcmaxouuu.github.io/picreg/reference/BinomialExample.md)
  : Small binary-classification dataset for the Binomial section of the
  vignette.
- [`CoxExample`](https://vcmaxouuu.github.io/picreg/reference/CoxExample.md)
  : Small survival dataset for the Cox section of the vignette.

## Internal

Overview pages and print methods.

- [`pic_methods`](https://vcmaxouuu.github.io/picreg/reference/pic_methods.md)
  : S3 methods for fitted pic objects.

- [`pic_plots`](https://vcmaxouuu.github.io/picreg/reference/pic_plots.md)
  : Plot methods for pic fits and diagnostics.

- [`pic_survival`](https://vcmaxouuu.github.io/picreg/reference/pic_survival.md)
  : Survival utilities for the Cox family.

- [`print(`*`<pic.family>`*`)`](https://vcmaxouuu.github.io/picreg/reference/print.pic.family.md)
  : Pretty-print a pic family descriptor.

- [`print(`*`<pic.pdb_asymptotic>`*`)`](https://vcmaxouuu.github.io/picreg/reference/print.pic.pdb_asymptotic.md)
  : Print PDB asymptotic diagnostic.

- [`print(`*`<pic.phase_transition>`*`)`](https://vcmaxouuu.github.io/picreg/reference/print.pic.phase_transition.md)
  : Print phase-transition analysis.

- [`print(`*`<summary.pic>`*`)`](https://vcmaxouuu.github.io/picreg/reference/print.summary.pic.md)
  :

  Print a `summary.pic` object.
