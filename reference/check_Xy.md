# Joint validation/preprocessing of `(X, y)`.

For survival data, rows of `X` and `y` are sorted by ascending time. For
a factor binary response, its two levels are returned as `y_levels`.

## Usage

``` r
check_Xy(
  X,
  y,
  y_kind,
  standardize_X = TRUE,
  X_mean = NULL,
  X_std = NULL,
  keep_feature_values = FALSE
)
```
