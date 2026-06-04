# Validate response according to its distributional kind.

Validate response according to its distributional kind.

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
