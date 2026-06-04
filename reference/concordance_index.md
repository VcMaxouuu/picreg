# Harrell's concordance index.

Harrell's concordance index.

## Usage

``` r
concordance_index(times, events, predictions)
```

## Arguments

- times:

  Survival times (length `n`).

- events:

  Event indicators (0 or 1; length `n`).

- predictions:

  Risk scores; higher = higher risk = shorter survival.

## Value

Numeric scalar in `[0, 1]`.
