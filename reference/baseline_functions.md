# Breslow baseline cumulative hazard and survival.

Breslow baseline cumulative hazard and survival.

## Usage

``` r
baseline_functions(times, events, predictions)
```

## Arguments

- times:

  Survival times (length `n`).

- events:

  Event indicators (0 or 1; length `n`).

- predictions:

  Risk scores; higher = higher risk = shorter survival.

## Value

A data.frame with columns `time`, `cumulative_hazard`, `survival`.
