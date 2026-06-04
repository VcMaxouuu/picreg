# Pivotal Detection Boundary regularization selector

Computes the data-driven regularization parameter \\\hat\lambda^{\rm
PDB}\_\alpha\\ using the Pivotal Detection Boundary (PDB) principle. The
selected value is defined as the empirical \\(1 - \alpha)\\ quantile of
a null-distribution gradient statistic \$\$\hat\lambda^{\rm PDB}\_\alpha
=q\_{1-\alpha}\left(\left\\\nabla \ell_0\right\\\_\infty \right),\$\$
where \\\ell_0\\ denotes the loss evaluated under the null model.

## Usage

``` r
lambda_pdb(
  X,
  family,
  n_simu = 5000L,
  alpha = 0.05,
  method = c("mc_exact", "mc_gaussian", "analytical")
)
```

## Arguments

- X:

  Numeric design matrix. Columns should typically be standardized to
  zero mean and unit variance.

- family:

  A PIC family specification. Can be either a family object or a
  character string accepted by
  [`get_family()`](https://vcmaxouuu.github.io/picreg/reference/get_family.md).

- n_simu:

  Number of Monte Carlo simulations used by the stochastic methods.
  Ignored when `method = "analytical"`.

- alpha:

  Nominal tail probability used to define the quantile level.

- method:

  One of "mc_exact", "mc_gaussian", or "analytical".

## Value

An object of class `"pic.lambda_pdb"`.

- value:

  Selected regularization parameter \\\hat\lambda^{\rm PDB}\_\alpha\\.

- statistics:

  Simulated null statistics used to estimate the quantile. `NULL` for
  the analytical method.

- method:

  Estimation method used.

- alpha:

  Quantile level parameter.

- n_simu:

  Number of Monte Carlo simulations.

## Details

Under the null \\\beta = 0\\, the gradient of the loss carries only
sampling noise. The smallest \\\lambda\\ large enough to dominate this
noise is the natural threshold separating signal from noise, i.e. the
value above which a coefficient should be kept rather than shrunk to
zero. Calibrating \\\hat\lambda\\ this way has two consequences. First,
the quantile depends only on the design matrix \\X\\, the family, and
the level \\\alpha\\ (not on the response \\y\\), so cross-validation is
no longer needed. Second, and more importantly, it leads to sharper
support recovery than prediction-error-based selectors such as
cross-validated lasso: by targeting the noise level of the gradient
directly, PDB controls the inclusion of noise variables and more
reliably identifies the true non-zero coefficients. Computing the
quantile requires only the distribution of \\\\\nabla
\ell_0\\\_\infty\\, which `lambda_pdb()` estimates through one of the
three methods below.

### Details on `method` option

The empirical quantile is obtained using one of the three following
methods:

- `"mc_exact"`:

  Family-aware Monte Carlo. For each of the `n_simu` draws, a response
  vector is sampled under the null model (\\\beta = 0\\) of the chosen
  family; the family-specific gradient of the loss at \\\beta = 0\\ is
  evaluated on the fixed design \\X\\, and its supremum norm is
  recorded. The empirical \\(1 - \alpha)\\ quantile of the `n_simu`
  recorded norms is returned. Most accurate but slowest of the three.

- `"mc_gaussian"`:

  Monte Carlo under the Gaussian approximation of the null gradient. A
  central-limit argument gives \\\nabla \ell_0 \approx \mathcal{N}(0,\\
  c(n)\\ \Sigma_X / n)\\ with \\\Sigma_X = X^\top X / n\\. Each of the
  `n_simu` draws samples directly from this Gaussian and records its
  supremum norm — no family-specific evaluation needed. Family-agnostic
  and noticeably faster than `"mc_exact"`; valid in the regime where the
  CLT kicks in (moderate to large \\n\\).

- `"analytical"`:

  Closed-form Bonferroni bound on the Gaussian tail. Combining a union
  bound over the \\p\\ coordinates of the gradient with the standard
  Gaussian tail bound gives \$\$\hat\lambda\_\alpha^{\rm analytical} =
  \Phi^{-1}\\\left(1 - \alpha / (2p)\right)\\ \sqrt{c(n) / n}.\$\$
  Deterministic and \\O(1)\\ — no simulation. Conservative (slightly
  over-estimates the true quantile) and gets looser as \\p\\ grows, but
  useful when speed matters or when Monte Carlo is overkill.

The `"mc_gaussian"` and `"analytical"` methods use a variance scaling
factor \\c(n)\\ depending on the family:

- Gaussian / Binomial / Poisson / Exponential: \\c(n) = 1\\

- Gumbel: \\c(n) = \exp(2(\gamma + 1))\\, with \\\gamma\\ the
  Euler-Mascheroni constant.

- Cox: \\c(n) = 1 / (4 \log n)\\

### Computational cost

Both Monte Carlo methods are dominated by a single \\p \times n\_{\rm
simu}\\ matrix product \\X^\top R\\, where \\R\\ stacks the simulated
residuals. This product is dispatched to BLAS as a single \\\tt gemm\\,
which is essentially as fast as a Monte Carlo selector can be on dense
designs. For very large \\n\\ or \\p\\, however, the constant becomes
large and `"mc_exact"` may be unnecessarily expensive:

- As \\n\\ grows, the central-limit approximation tightens and
  `"mc_gaussian"` gives essentially the same \\\hat\lambda\\ as
  `"mc_exact"` at a fraction of the cost (no family-specific residual
  generation, fewer dependencies on \\y\\-draws).

- `"analytical"` is \\O(1)\\ and useful as a quick upper bound —
  slightly conservative, but accurate enough for triage and for very
  high \\p\\ where the Bonferroni tail is tight.

Practical rule of thumb: prefer `"mc_exact"` for small to moderate
problems (default), `"mc_gaussian"` once \\n\\ grows and the design is
dense, and `"analytical"` when even the Monte Carlo cost is a concern.

### Details on `alpha` option

The level \\\alpha\\ is the nominal Type-I error of the test of
\\H_0\\\colon \beta = 0\\. By construction of the quantile,
\$\$\Pr\\\left(\left\\\nabla \ell_0\right\\\_\infty \> \hat\lambda^{\rm
PDB}\_\alpha \\\big\|\\ H_0\right) = \alpha,\$\$ so under the null model
no variable enters the active set with probability at most \\\alpha\\.
With the default \\\alpha = 0.05\\, this caps the false-discovery rate
at \\5\\\\ under the null: when the data carry no signal, picreg returns
the empty support \\95\\\\ of the time.

## Examples

``` r
X <- scale(matrix(rnorm(100 * 20), 100, 20))
lam <- lambda_pdb(
  X,
  family = "gaussian",
  method = "mc_exact"
)
#> Warning: lambda_pdb: design matrix X does not appear standardized; results may be unreliable.
print(lam)
#> pic lambda_pdb (method = mc_exact, alpha = 0.05)
#>   lambda = 0.295355
```
