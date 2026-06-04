# Plot the PDB null distribution.

Histogram of the simulated null gradient-norm statistics, with a
vertical line at the selected \\\hat\lambda\\.

## Usage

``` r
# S3 method for class 'pic.lambda_pdb'
plot(x, breaks = 40L, ...)
```

## Arguments

- x:

  A `pic.lambda_pdb` object (typically `fit$lambda_pdb`).

- breaks:

  Number of histogram bins (default 40).

- ...:

  Unused; present for S3 method consistency.

## Value

Invisibly returns `x`.
