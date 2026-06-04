# picreg

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/picreg)](https://CRAN.R-project.org/package=picreg)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/picreg)](https://CRAN.R-project.org/package=picreg)
[![License: GPL-2](https://img.shields.io/badge/license-GPL--2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
<!-- badges: end -->

**Variable selection using the Pivotal Information Criterion.**

Sparse regression and classification via the Pivotal Information
Criterion (PIC), an alternative to BIC, cross-validation, and
Lasso-based tuning. The regularization parameter is selected from a
pivotal null-distribution statistic, eliminating the need for
cross-validation and yielding sharper support recovery.

Provides FISTA optimization for the L1, SCAD, and MCP penalties across
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

Install the released version from CRAN:

```r
install.packages("picreg")
```

Or the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("VcMaxouuu/picreg")
```

## Quick start

```r
library(picreg)

data(QuickStartExample)
fit <- pic(QuickStartExample$X, QuickStartExample$y)

fit$selected   # names of the selected variables
fit$lambda     # PDB-selected lambda (no cross-validation)
summary(fit)   # family, penalty, lambda, and non-zero coefficients
coef(fit)      # coefficients (sparse matrix, original scale of X)
predict(fit, newx = QuickStartExample$X[1:5, ])
```

The design deliberately mirrors `glmnet`: a single `pic()` fitting
function returning an object equipped with `print()`, `summary()`,
`coef()`, `predict()`, `plot()`, and `assess()` methods that behave
consistently across all six families.

## Documentation

The full walk-through — fitting across all six families, predicting,
visualizing, choosing penalties, and running diagnostics
(`phase_transition()`, `pdb_asymptotic()`) — lives in the package
vignette:

```r
vignette("vignette", package = "picreg")
```

## Reference

Sardy, van Cutsem, and van de Geer. *The Pivotal Information
Criterion.* <https://arxiv.org/abs/2603.04172>
(<doi:10.48550/arXiv.2603.04172>)
