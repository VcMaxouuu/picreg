# Harrell's concordance index.

A pair \\(i, j)\\ is comparable when \\i\\ has an event and either \\t_i
\< t_j\\, or \\t_i = t_j\\ and \\j\\ is censored. It is concordant when
\\i\\ has the higher risk score; tied scores count \\1/2\\. Computed in
\\O(n \log n)\\.

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

Numeric scalar in `[0, 1]` (`0` when no pair is comparable).
