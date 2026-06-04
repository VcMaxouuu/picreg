# Small Gaussian dataset for the introductory vignette.

A synthetic Gaussian regression dataset used to illustrate
[`pic()`](https://vcmaxouuu.github.io/picreg/reference/pic.md)
throughout the introductory vignette. It contains \\n = 100\\
observations of \\p = 30\\ predictors, of which only \\s = 5\\ carry
signal; the remaining 25 are pure noise.

## Usage

``` r
QuickStartExample
```

## Format

A list with two components:

- `X`:

  Numeric matrix of dimension \\100 \times 30\\ with column names
  `gene_*` and `noise_*`.

- `y`:

  Numeric vector of length \\100\\.

## Details

Column names are chosen to make the underlying support obvious at a
glance:

- `gene_1, ..., gene_5`: the five active variables, whose non-zero
  coefficients are drawn uniformly in \\\[0.5,\\ 1.5\]\\ with random
  sign.

- `noise_1, ..., noise_25`: the remaining inactive variables, with true
  coefficient \\0\\.

The columns are interleaved in random order; column names are the only
indicator of which features are part of the true support.

The response is generated as \\y = X\beta + \varepsilon\\ with
\\\varepsilon \sim \mathcal{N}(0, 1)\\.

## Examples

``` r
data(QuickStartExample)
fit <- pic(QuickStartExample$X, QuickStartExample$y)
fit$selected
#> [1] "gene_1" "gene_2" "gene_3" "gene_4" "gene_5"
```
