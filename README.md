# picreg

**Variable Selection using the Pivotal Information Criterion.**

Sparse regression and classification in high dimensions with **automatic
regularisation parameter selection** via the Pivotal Detection Boundary
(PDB) method. FISTA optimisation, L1 / SCAD / MCP penalties across six
families (Gaussian, Binomial, Poisson, Exponential, Gumbel, Cox).

The PDB selector chooses `lambda` as a quantile of a pivotal
null-distribution statistic — no cross-validation loop required.

## Installation

```r
# from a local clone
install.packages("picreg", repos = NULL, type = "source")
```

## Quick start

```r
library(picreg)

set.seed(1)
n <- 100; p <- 100; s <- 5
X <- matrix(rnorm(n * p), n, p)
beta <- numeric(p); beta[sample.int(p, s)] <- 3
y <- as.numeric(X %*% beta + rnorm(n))

fit <- pic(X, y, family = "gaussian", penalty = "scad")
coef(fit)[which(coef(fit) != 0)]
fit$lambda        # PDB-selected lambda
fit$selected      # indices of non-zero coefficients
```

## Supported families and penalties

| Family        | `family =`       | Response       |
| ------------- | ---------------- | -------------- |
| Gaussian      | `"gaussian"`     | continuous     |
| Binomial      | `"binomial"`     | 0/1 binary     |
| Poisson       | `"poisson"`      | count          |
| Exponential   | `"exponential"`  | positive cont. |
| Gumbel        | `"gumbel"`       | continuous     |
| Cox PH        | `"cox"`          | (time, event)  |

| Penalty | `penalty =`  | Tuning parameter            |
| ------- | ------------ | --------------------------- |
| Lasso   | `"lasso"`    | —                           |
| SCAD    | `"scad"`     | `scad_a = 3.7` (default)    |
| MCP     | `"mcp"`      | `mcp_gamma = 3.0` (default) |

## Diagnostics

Two empirical-check tools are exported alongside the core fit:

- `phase_transition()` — support-recovery probability as a function of
  the sparsity level `s`.
- `pdb_asymptotic()` — convergence of the family-exact null distribution
  to the Gaussian approximation as `n` grows.

See `vignette("picr-diagnostics", package = "picreg")`.

## Reference

See <doi:10.48550/arXiv.2603.04172>.
