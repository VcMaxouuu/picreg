# Small binary-classification dataset for the Binomial section of the vignette.

A synthetic logistic-regression dataset used to illustrate
[`pic()`](https://vcmaxouuu.github.io/picreg/reference/pic.md) with the
Binomial family. It contains \\n = 300\\ observations of \\p = 50\\
predictors, of which \\s = 5\\ carry signal; the remaining 45 are noise.

## Usage

``` r
BinomialExample
```

## Format

A list with two components:

- `X`:

  Numeric matrix of dimension \\300 \times 50\\ with column names
  `gene_*` and `noise_*`.

- `y`:

  Integer vector of length \\300\\ containing \\0/1\\ class labels.

## Details

Column names follow the same convention as
[`QuickStartExample`](https://vcmaxouuu.github.io/picreg/reference/QuickStartExample.md):
active variables are labeled `gene_1, ..., gene_5` and inactive ones
`noise_1, ..., noise_45`, interleaved in random order.

The binary response is generated as \$\$Y_i \sim
\mathrm{Bernoulli}\\\left(\frac{1}{1 + e^{-X_i\beta}}\right),\$\$ with
non-zero coefficients drawn uniformly in \\\[1.5,\\ 3\]\\ with random
sign. The class balance is roughly \\45\\\\ of positives.

## Examples

``` r
data(BinomialExample)
fit <- pic(BinomialExample$X, BinomialExample$y, family = "binomial")
fit$selected
#> [1] "gene_1" "gene_2" "gene_3" "gene_4" "gene_5"
```
