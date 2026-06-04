# Phase-transition analysis of support recovery.

For a fixed `(n, p)` configuration, varies the sparsity level `s` from
`0` to `s_max` and, for each `s`, estimates by Monte Carlo (`m`
replications) the probability that
[`pic()`](https://vcmaxouuu.github.io/picreg/reference/pic.md) recovers
exactly the support of the true coefficient vector. If several penalties
are supplied, one result is produced for each `(n, p, penalty)`
combination.

## Usage

``` r
phase_transition(
  n,
  p,
  type = c("gaussian", "binomial", "poisson", "exponential", "gumbel", "cox"),
  s_max,
  m = 50,
  penalty = "lasso",
  beta_value = 3,
  lambda_method = "mc_exact",
  lambda_alpha = 0.05,
  lambda_n_simu = 5000L,
  verbose = FALSE,
  parallel = FALSE,
  workers = parallel::detectCores() - 1L
)
```

## Arguments

- n:

  Integer or integer vector — number of observations per configuration.

- p:

  Integer or integer vector of the same length as `n` — number of
  features.

- type:

  Family name: `"gaussian"`, `"binomial"`, `"poisson"`, `"exponential"`,
  `"gumbel"`, or `"cox"`.

- s_max:

  Largest sparsity level evaluated. Must satisfy `s_max < min(p)`.

- m:

  Number of Monte Carlo replications per `s`.

- penalty:

  One or more penalties for
  [`pic()`](https://vcmaxouuu.github.io/picreg/reference/pic.md):
  `"lasso"`, `"scad"`, or `"mcp"`.

- beta_value:

  Magnitude of the non-zero coefficients used to generate the response.
  The sign is fixed to `+`.

- lambda_method:

  Passed to
  [`pic()`](https://vcmaxouuu.github.io/picreg/reference/pic.md).
  Default `"mc_exact"`.

- lambda_alpha:

  Nominal level for the PDB selector.

- lambda_n_simu:

  Monte Carlo size for the PDB selector.

- verbose:

  Logical; if `TRUE`, prints a one-line progress message per
  `(n, p, penalty, s)`.

- parallel:

  Logical; if `TRUE`, distribute the `m` Monte Carlo replications of
  each `(n, p, penalty, s)` cell across multiple R processes via the
  `future` framework
  ([`future::multisession`](https://future.futureverse.org/reference/multisession.html)
  plan).

- workers:

  Integer; number of background R processes to use when
  `parallel = TRUE`. Ignored otherwise.

## Value

An object of class `c("pic.phase_transition", "pic.diagnostic")`.

- s_grid:

  `0:s_max`.

- exact_recovery, tpr, fdr:

  Matrices of shape `(length(n) * length(penalty), s_max + 1)` - one row
  per `(n, p, penalty)` curve, one column per sparsity level

- curve_n, curve_p, curve_penalty:

  Curve descriptors aligned with the rows of the metric matrices.

- config:

  A list of all configuration arguments for downstream plotting /
  reporting.

- call:

  The call.

## Details

At each Monte Carlo replicate a design matrix is drawn from a standard
Gaussian, `s` features are sampled uniformly at random, the response is
generated under the chosen `type`, and one `pic` fit is run for each
requested penalty. The selected support is compared to the truth and
three metrics are stored:

- `exact_recovery`:

  1 if the selected set equals the true set, 0 otherwise.

- `tpr`:

  \\\left\|\hat{S}\cap S\right\| / \|S\|\\ - true positive rate. For
  `s = 0`, this is set to 1 when no true feature exists.

- `fdr`:

  \\\left\|\hat{S}\setminus S\right\| / \max{(\|\hat{S}\|, 1)}\\ - false
  discovery rate.

### Details on data sampling

At each replicate, the design \\X\\ is drawn iid \\\mathcal{N}(0, 1)\\,
the true support \\S\\ is sampled uniformly at random, and \\\beta_j =
\\ `beta_value` for \\j \in S\\, \\\beta_j = 0\\ otherwise. The linear
predictor is \\\eta = X\beta\\. The response \\y\\ is then drawn
conditionally on \\\eta\\ according to the requested family:

- `"gaussian"`:

  \\y_i = \eta_i + \varepsilon_i\\ with \\\varepsilon_i \sim
  \mathcal{N}(0, 1)\\.

- `"binomial"`:

  \\y_i \sim \mathrm{Bernoulli}\\\left(\sigma(\eta_i)\right)\\, where
  \\\sigma(z) = 1 / (1 + e^{-z})\\ is the logistic function.

- `"poisson"`:

  \\y_i \sim \mathrm{Poisson}\\\left(e^{\eta_i}\right)\\.

- `"exponential"`:

  \\y_i \sim \mathrm{Exp}\\\left(\mathrm{rate} = e^{\eta_i}\right)\\.

- `"gumbel"`:

  \\y_i = \eta_i + \varepsilon_i\\ with \\\varepsilon_i \sim
  \mathrm{Gumbel}(0, 1)\\, drawn as \\-\log(-\log U_i)\\ for \\U_i \sim
  \mathcal{U}(0, 1)\\.

- `"cox"`:

  Event times \\T_i \sim \mathrm{Exp}\\\left(e^{\eta_i}\right)\\ and
  independent censoring times \\C_i \sim \mathrm{Exp}(1)\\. The response
  is the 2-column matrix \\\bigl(\min(T_i, C_i),\\ \mathbf{1}\\T_i \le
  C_i\\\bigr)\\.

## See also

[`plot.pic.phase_transition()`](https://vcmaxouuu.github.io/picreg/reference/plot.pic.phase_transition.md)
for visualization.

## Examples

``` r
# \donttest{
pt <- phase_transition(n = 50, p = 100, type = "gaussian",
                       s_max = 8, m = 20,
                       penalty = c("lasso", "scad"), 
                       parallel = TRUE)
plot(pt)

# }
```
