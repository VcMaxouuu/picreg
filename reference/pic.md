# Sparse linear regression using the Pivotal Information Criterion.

Fits sparse regression, classification, and survival models built on a
linear predictor, using a family-specific loss whose gradient at the
null is (asymptotically) pivotal, combined with a sparsity-inducing
penalty. Supported families are Gaussian, binomial, Poisson,
exponential, Gumbel, and Cox proportional hazards.

## Usage

``` r
pic(
  X,
  y,
  family = c("gaussian", "binomial", "poisson", "exponential", "gumbel", "cox"),
  penalty = "lasso",
  standardize = TRUE,
  intercept = TRUE,
  lambda_method = c("mc_exact", "mc_gaussian", "analytical"),
  relax = FALSE,
  mcp_gamma = 3,
  scad_a = 3.7,
  lambda = NULL,
  lambda_alpha = 0.05,
  lambda_n_simu = 2000,
  tol = 1e-05,
  path_length = 10L,
  maxit = 1000
)
```

## Arguments

- X:

  Numeric design matrix of shape `(n, p)`, where rows are observations
  and columns are candidate predictors. `X` may be a matrix or a data
  frame coercible to a numeric matrix. It must contain at least two
  columns and no `NA`, `NaN`, or infinite values.

- y:

  Response variable with length `n`. For Gaussian and Gumbel models, `y`
  must be numeric. For Binomial models, `y` must contain only `0` and
  `1`. For Poisson and Exponential models, `y` must be non-negative. For
  `family = "cox"`, `y` must be a two-column numeric matrix
  `(time, event)`, where `time` is non-negative and `event` is coded as
  `0` or `1`.

- family:

  A character name determining the response distribution and the
  associated loss function used during optimization. See **Details on
  family-specific losses** below for the exact objective functions
  associated with each family. Default `"gaussian"`.

- penalty:

  Penalty name: one of `"lasso"`, `"scad"`, or `"mcp"` (lowercase).
  Default `"lasso"`. See
  [pic_penalties](https://vcmaxouuu.github.io/picreg/reference/pic_penalties.md)
  for the precise form of each penalty and the role of `scad_a` /
  `mcp_gamma`.

- standardize:

  Whether to standardize columns of `X` to zero mean and unit variance
  prior to computing the PDB regularization parameter and fitting the
  model. Strongly recommended for proper PDB calibration. When `TRUE`,
  the optimization is performed on the standardized design matrix and
  the stored fitted coefficients therefore correspond to the
  standardized scale. Use
  [`coef.pic()`](https://vcmaxouuu.github.io/picreg/reference/coef.pic.md)
  to recover coefficients on the original scale of `X`. Default is
  `TRUE`. Standardization centers `X`; for non-Cox families the
  resulting column-mean shift is absorbed by the intercept, so combining
  `standardize = TRUE` with `intercept = FALSE` assumes the data are
  already centered (see `intercept`).

- intercept:

  Whether to fit an unpenalized intercept. Forced to `FALSE` for the Cox
  family. We assume the data has been centered if intercept = FALSE.

- lambda_method:

  One of `"mc_exact"`, `"mc_gaussian"`, `"analytical"`. See
  [`lambda_pdb()`](https://vcmaxouuu.github.io/picreg/reference/lambda_pdb.md)
  for details. The default `"mc_exact"` is a safe choice for small to
  moderate designs; for very large \\n\\ or \\p\\, `"mc_gaussian"`
  matches `"mc_exact"` closely at a fraction of the cost (it amortises
  only one BLAS \\\tt gemm\\, no family-specific residual draws), and
  `"analytical"` is closed-form.

- relax:

  If `TRUE`, run an unpenalized refit on the selected support after the
  regularization path (debiasing).

- mcp_gamma:

  MCP concavity parameter when `penalty = "mcp"`. Default 3.0.

- scad_a:

  SCAD concavity parameter when `penalty = "scad"`. Default 3.7.

- lambda:

  Optional user-supplied regularization parameter. When `NULL`,
  [`lambda_pdb()`](https://vcmaxouuu.github.io/picreg/reference/lambda_pdb.md)
  is used to select it automatically. Providing `lambda` avoids
  recomputing the PDB calibration, which is useful for example when
  fitting several penalties (`"lasso"`, `"scad"`, `"mcp"`) on the same
  dataset, since the PDB choice of \\\lambda\\ does not depend on the
  penalty itself.

- lambda_alpha:

  Nominal level for PDB. See
  [`lambda_pdb()`](https://vcmaxouuu.github.io/picreg/reference/lambda_pdb.md).

- lambda_n_simu:

  Number of Monte Carlo draws for PDB. See
  [`lambda_pdb()`](https://vcmaxouuu.github.io/picreg/reference/lambda_pdb.md).

- tol:

  FISTA convergence tolerance at the final path step.

- path_length:

  Number of points in the warm-start regularization path running from
  \\\lambda\_{\max}\\ down to \\\lambda^{\rm PDB}\\. Default 10.

- maxit:

  Maximum number of iterations for FISTA if convergence is not yet
  reached.

## Value

An object of class `c("pic.<family>", "pic")`.

- beta:

  Fitted coefficients (length `p`).

- intercept:

  Fitted intercept or `NULL`.

- df:

  The number of nonzero coefficients.

- selected:

  Indices of the nonzero coefficients.

- lambda:

  \\\lambda\\ used during training (PDB or user-supplied).

- family:

  The family object used.

- penalty:

  The penalty object used.

- lambda_pdb:

  Result from
  [`lambda_pdb()`](https://vcmaxouuu.github.io/picreg/reference/lambda_pdb.md)
  if \\\lambda\\ is not user-supplied.

- n_samples:

  The number of observations in the training set.

- n_features:

  The number of variables in the training set.

For `family = "cox"`, the fit additionally carries:

- baseline_cumulative_hazard:

  Estimated Breslow baseline cumulative hazard function.

- baseline_survival:

  Estimated baseline survival function.

- unique_times:

  Sorted unique event times used in the baseline estimation.

- censoring_rate:

  Percentage of censored observations in the training set.

## Details

Minimises the objective function \$\$\hat\beta = \arg\min\_\beta\\
L(\beta) \\+\\ \lambda\_\alpha^{\rm PDB}\\\mathrm{pen}(\beta)
\quad\text{with}\quad L(\beta) =
\phi\\\left(\ell_n\\\left(g(\beta)\right)\right),\$\$ where the
regularization parameter is selected automatically using the Pivotal
Detection Boundary (PDB) principle implemented in
[`lambda_pdb()`](https://vcmaxouuu.github.io/picreg/reference/lambda_pdb.md).
The PDB choice is not just a substitute for cross-validation: by
calibrating \\\lambda\\ against a pivotal null-distribution statistic
rather than out-of-sample prediction error, it yields sharper support
recovery - fewer false positives on noise variables and more accurate
identification of the true non-zero coefficients. Here, \\\phi\\ and
\\g\\ are family-specific transformations chosen so that the gradient at
the null has a distribution free of the nuisance parameter, which is
precisely what makes the use of \\\lambda^{\rm PDB}\_\alpha\\ valid.
Optimization is performed using a warm-started FISTA regularization
path.

### Details on family-specific losses

- `"gaussian"`:

  Uses the mean squared error for \\\ell_n\\, with \\\phi(\cdot) =
  \sqrt{\cdot}\\ and \\g(\cdot) = \mathrm{Id}\\, giving the square-root
  lasso loss \$\$L(\beta) = \sqrt{\frac{1}{n}\sum\_{i=1}^n \left(y_i -
  (\beta_0 + \mathbf{x}\_i^\top\beta)\right)^2}.\$\$ This makes the
  gradient at \\\beta = 0\\ scale-free, so the PDB threshold needs no
  estimate of the noise standard deviation.

- `"binomial"`:

  A variance-stabilized transformation of the Bernoulli likelihood. With
  \\\phi(\cdot) = \mathrm{Id}\\ and the logistic link \\\theta_i =
  g(\eta_i) = (1 + e^{-\eta_i})^{-1}\\, the loss is \$\$L(\beta) =
  \frac{1}{n}\sum\_{i=1}^n \left(2 y_i \sqrt{\frac{1 -
  \theta_i}{\theta_i}} + 2 (1 - y_i) \sqrt{\frac{\theta_i}{1 -
  \theta_i}}\right).\$\$ The classical logistic link is preserved; only
  the loss itself is modified to obtain a pivotal gradient at the null.

- `"poisson"`:

  The same pivotalization applied to the Poisson likelihood. With
  \\\phi(\cdot) = \mathrm{Id}\\ and the canonical log link \\\theta_i =
  g(\eta_i) = e^{\eta_i}\\, the loss is \$\$L(\beta) =
  \frac{1}{n}\sum\_{i=1}^n \left(\frac{2 y_i}{\sqrt{\theta_i}} + 2
  \sqrt{\theta_i}\right).\$\$ The canonical log link is kept unchanged;
  pivotality is obtained through the loss rather than through the link.

- `"exponential"`:

  Uses the standard Exponential negative log-likelihood directly (no
  transformation). With \\\phi(\cdot) = \mathrm{Id}\\ and \\\theta_i =
  g(\eta_i) = e^{\eta_i}\\, the loss is \$\$L(\beta) =
  \frac{1}{n}\sum\_{i=1}^n \left(\log\theta_i +
  \frac{y_i}{\theta_i}\right).\$\$

- `"gumbel"`:

  A location-scale model with extreme-value noise. The base
  log-likelihood is \$\$\ell_n(\theta, \sigma) = \log(\sigma) +
  \frac{1}{n}\sum\_{i=1}^n \left(z_i + e^{-z_i}\right), \qquad z_i =
  \frac{y_i - \theta_i}{\sigma}.\$\$ With \\\phi(\cdot) = \exp(\cdot)\\
  and the identity link \\\theta_i = g(\eta_i) = \eta_i\\, the objective
  is \\L(\beta, \sigma) = \exp(\ell_n(\theta, \sigma))\\. The scale
  \\\sigma\\ is re-estimated internally by maximum likelihood at every
  iteration; the user does not supply it.

- `"cox"`:

  A square-root-transformed Cox partial log-likelihood. With
  \\\phi(\cdot) = \sqrt{\cdot}\\ and the identity link \\\theta_i =
  g(\eta_i) = \eta_i\\, the objective is \$\$L(\beta) =
  \sqrt{-\frac{1}{n}\left( \sum\_{i=1}^n \delta_i \eta_i - \sum\_{i=1}^n
  \delta_i \log\left(\sum\_{j \in R_i} e^{\eta_j}\right)\right)},\$\$
  where \\\delta_i\\ is the event indicator and \\R_i\\ the risk set at
  event time \\t_i\\. The Breslow approximation is used for tied event
  times, and the model is always fitted without an intercept.

## Examples

``` r
data(QuickStartExample)
X <- QuickStartExample$X
y <- QuickStartExample$y
fit <- pic(X, y, family = "gaussian", penalty = "scad")
fit$beta
#>             [,1]
#>  [1,] -0.6371864
#>  [2,]  0.0000000
#>  [3,]  0.0000000
#>  [4,]  0.0000000
#>  [5,] -0.2308360
#>  [6,]  0.0000000
#>  [7,]  0.0000000
#>  [8,]  0.0000000
#>  [9,]  0.0000000
#> [10,] -1.1585369
#> [11,]  0.0000000
#> [12,]  0.0000000
#> [13,]  0.0000000
#> [14,]  0.0000000
#> [15,]  0.0000000
#> [16,]  0.0000000
#> [17,] -0.9759369
#> [18,]  0.0000000
#> [19,]  0.0000000
#> [20,]  0.0000000
#> [21,]  0.0000000
#> [22,]  0.0000000
#> [23,]  0.0000000
#> [24,]  0.0000000
#> [25,]  0.2894805
#> [26,]  0.0000000
#> [27,]  0.0000000
#> [28,]  0.0000000
#> [29,]  0.0000000
#> [30,]  0.0000000
fit$selected
#> [1] "gene_1" "gene_2" "gene_3" "gene_4" "gene_5"
```
