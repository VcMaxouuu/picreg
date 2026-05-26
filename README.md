# picreg

**Variable selection using the Pivotal Information Criterion.**

Sparse regression and classification via the Pivotal Information
Criterion (PIC), an alternative to BIC, cross-validation, and
Lasso-based tuning. The regularisation parameter is selected from a
pivotal null-distribution statistic, eliminating the need for
cross-validation and yielding sharper support recovery.

Provides FISTA optimisation for the L1, SCAD, and MCP penalties across
six response distributions:

| Family      | `family =`      | Response       |
| ----------- | --------------- | -------------- |
| Gaussian    | `"gaussian"`    | continuous     |
| Binomial    | `"binomial"`    | 0/1 binary     |
| Poisson     | `"poisson"`     | count          |
| Exponential | `"exponential"` | positive cont. |
| Gumbel      | `"gumbel"`      | continuous     |
| Cox PH      | `"cox"`         | (time, event)  |

Under standard sparsity assumptions, the selector achieves a phase
transition for exact support recovery, analogous to results in
compressed sensing.

## Installation

The development version can be installed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("VcMaxouuu/picreg")
```

## Quick start

```r
library(picreg)

data(QuickStartExample)
fit <- pic(QuickStartExample$X, QuickStartExample$y)

fit$selected   # names of selected variables
fit$lambda     # PDB-selected lambda
coef(fit)      # coefficients on the original scale
```

## Documentation

The full walk-through - fitting across all six families, predicting,
visualising, choosing penalties, and running diagnostics
(`phase_transition()`, `pdb_asymptotic()`) - lives in the package
vignette:

```r
vignette("vignette", package = "picreg")
```

## Reference

Sardy, van Cutsem, and van de Geer. *The Pivotal Information
Criterion.* <https://arxiv.org/abs/2603.04172>
(<doi:10.48550/arXiv.2603.04172>)
