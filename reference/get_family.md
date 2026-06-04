# Resolve a family name into a pic family descriptor.

Accepts one of the six supported names: `"gaussian"`, `"binomial"`,
`"poisson"`, `"exponential"`, `"gumbel"`, `"cox"`. Anything else raises
an error — the build does not support user-defined families.
Already-built descriptors are returned unchanged (internal re-entry).

## Usage

``` r
get_family(family)
```

## Arguments

- family:

  A family name or already-built pic family descriptor.

## Value

A pic family descriptor.
