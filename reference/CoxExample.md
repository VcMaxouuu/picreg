# Small survival dataset for the Cox section of the vignette.

A synthetic survival dataset used to illustrate
[`pic()`](https://vcmaxouuu.github.io/picreg/reference/pic.md) with the
Cox family. It contains \\n = 250\\ subjects observed on \\p = 50\\
covariates, of which \\s = 5\\ carry signal; the remaining 45 are noise.

## Usage

``` r
CoxExample
```

## Format

A list with two components:

- `X`:

  Numeric matrix of dimension \\250 \times 50\\ with column names
  `gene_*` and `noise_*`.

- `y`:

  Numeric matrix of dimension \\250 \times 2\\ with columns `time` and
  `event`.

## Details

Column names follow the same convention as
[`QuickStartExample`](https://vcmaxouuu.github.io/picreg/reference/QuickStartExample.md):
active variables are labeled `gene_1, ..., gene_5` and inactive ones
`noise_1, ..., noise_45`, interleaved in random order.

Event times are drawn from an exponential proportional-hazards model
\$\$T_i \sim \mathrm{Exp}\\\bigl(e^{X_i\beta}\bigr),\$\$ and independent
censoring times from \\C_i \sim \mathrm{Exp}(0.5)\\. The observed
response is the standard two-column \\(\min(T_i, C_i),\\ \mathbf{1}\\T_i
\le C_i\\)\\. The censoring rate is roughly \\40\\\\.

## Examples

``` r
data(CoxExample)
fit <- pic(CoxExample$X, CoxExample$y, family = "cox")
fit$selected
#> [1] "gene_1" "gene_2" "gene_3" "gene_4" "gene_5"
```
