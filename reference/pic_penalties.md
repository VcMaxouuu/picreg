# Sparsity-inducing penalties for pic.

Three penalties are supported, identified by lowercase name to match the
C++ registry. Each penalty enters the
[`pic()`](https://vcmaxouuu.github.io/picreg/reference/pic.md) objective
as \$\$\mathrm{pen}(\beta) = \sum\_{j=1}^p p\_\lambda(\|\beta_j\|),\$\$
where \\p\_\lambda(\cdot)\\ depends on the penalty.

## Details

- `"lasso"`:

  L1 (soft-thresholding) penalty: \$\$p\_\lambda(\|t\|) = \lambda
  \|t\|.\$\$ Convex, gives the strongest shrinkage on large coefficients
  bias does not vanish as \\\|t\| \to \infty\\.

- `"scad"` (Smoothly Clipped Absolute Deviation, Fan & Li 2001):

  Non-convex penalty with concavity parameter `scad_a > 2` (default
  3.7): \$\$p\_\lambda'(\|t\|) = \lambda\\\left\\ \mathbf{1}\\\|t\| \le
  \lambda\\ + \frac{(a\lambda - \|t\|)\_+}{(a - 1)\lambda}
  \mathbf{1}\\\|t\| \> \lambda\\\right\\.\$\$ Behaves like the lasso for
  small \\\|t\|\\, then tapers off so large coefficients are barely
  penalized - yields nearly unbiased estimates on strong signals.

- `"mcp"` (Minimax Concave Penalty, Zhang 2010):

  Non-convex penalty with concavity parameter `mcp_gamma > 1` (default
  3.0): \$\$p\_\lambda'(\|t\|) = \left(\lambda -
  \frac{\|t\|}{\gamma}\right)\_+.\$\$ Similar motivation as SCAD but a
  smoother transition: starts at the lasso derivative for small
  \\\|t\|\\ and tapers linearly to zero at \\\|t\| = \gamma\lambda\\.

The actual evaluation and proximal operators live in C++
(`src/penalty_*.cpp`). Larger `scad_a` / `mcp_gamma` make the penalty
closer to the lasso; smaller values amplify the non-convexity (and the
bias reduction on strong signals).
