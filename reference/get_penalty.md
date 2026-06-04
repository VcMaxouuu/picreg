# Resolve a penalty name into a pic penalty descriptor.

Accepts one of `"lasso"`, `"scad"`, `"mcp"` (case-insensitive). Anything
else raises an error: the build does not support user-defined penalties.

## Usage

``` r
get_penalty(penalty, scad_a = 3.7, mcp_gamma = 3)
```

## Arguments

- penalty:

  A penalty name.

- scad_a:

  SCAD concavity parameter (default 3.7).

- mcp_gamma:

  MCP concavity parameter (default 3.0).

## Value

A `pic.penalty` descriptor.
