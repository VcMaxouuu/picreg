# Validate response according to its distributional kind.

A binary response may also be logical (`FALSE`/`TRUE` mapped to `0`/`1`)
or a two-level factor (first level mapped to `0`, second to `1`).

## Usage

``` r
check_y(y, n_samples = NULL, y_kind = "continuous")
```

## Arguments

- y:

  Response vector or 2-column matrix (survival).

- n_samples:

  Expected number of observations.

- y_kind:

  One of `"continuous"`, `"positive"`, `"binary"`, `"survival"`.

## Value

Validated `y` (possibly coerced to matrix for survival).
