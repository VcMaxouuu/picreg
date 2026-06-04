# Summary of the PDB lambda selector.

Prints a formatted summary of the Pivotal Detection Boundary selector
used by [`pic()`](https://vcmaxouuu.github.io/picreg/reference/pic.md)
to choose \\\hat\lambda\\: method, nominal level, Monte Carlo size,
selected \\\hat\lambda\\, and a compact view of the null distribution
when Monte Carlo was run. For models fitted with `method = "analytical"`
or with a user-supplied `lambda`, only the selector metadata is shown.

## Usage

``` r
pdb_summary(model, digits = 4L)
```

## Arguments

- model:

  A fitted `pic` object.

- digits:

  Number of significant digits used in the distribution table.

## Value

Invisibly returns `NULL`. Called for its side effect (printing).

## Examples

``` r
data(QuickStartExample)
fit <- pic(QuickStartExample$X, QuickStartExample$y,
           family = "gaussian", penalty = "lasso")
pdb_summary(fit)
#> PDB lambda selector
#> -------------------
#>   method       : mc_exact
#>   alpha        : 0.05
#>   n_simu       : 2,000
#>   lambda_hat   : 0.3117
#> 
#>   Null distribution:
#>    min     q05     q25  median     q75     q95     max  
#> 0.1154  0.1687  0.2018  0.2271  0.2592  0.3117  0.4096  
#> 
#>   mean = 0.2326      sd = 0.0437
```
