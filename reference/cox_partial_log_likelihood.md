# Breslow partial log-likelihood (negative, normalized by `n`).

Breslow partial log-likelihood (negative, normalized by `n`).

## Usage

``` r
cox_partial_log_likelihood(times, events, predictions)
```

## Arguments

- times:

  Survival times (length `n`).

- events:

  Event indicators (0 or 1; length `n`).

- predictions:

  Risk scores; higher = higher risk = shorter survival.

## Value

A numeric scalar: the negative Breslow partial log-likelihood for the
Cox proportional-hazards model, normalized by the sample size `n`. Lower
values indicate a better fit.
