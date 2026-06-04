# GLM families for pic — descriptor layer.

Each family is a thin descriptor carrying:

- `name` - used by the dispatchers to route to the C++ implementation.

- `g` - mean-function link (object with `name` and callable `fn`).
  `g$fn(eta)` is applied at predict time to obtain the mean response.

- `phi` - variance-stabilizing transform (object with `name`). Purely
  informational on the R side; the actual stabilization happens inside
  the C++ math.

## Details

Typing `fit$family` gives a one-glance summary of the model's family
(see
[`print.pic.family()`](https://vcmaxouuu.github.io/picreg/reference/print.pic.family.md)).
The actual loss / gradient math lives in C++:

- `src/family_*.cpp` — Gaussian, Binomial, Poisson, Exponential, Gumbel.

- `src/cox.cpp` — Cox.
