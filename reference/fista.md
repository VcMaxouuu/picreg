# Non-monotone Accelerated Proximal Gradient (nmAPG) solver — dispatcher.

Thin R wrapper around the C++ `fista_cpp` implementation. The actual
algorithm (Li & Lin, NeurIPS 2015) lives in `src/fista.cpp`; this
function only extracts the C++-facing arguments and routes the call.

## Usage

``` r
fista(
  X,
  y,
  family,
  penalty,
  lambda_reg,
  fit_intercept = TRUE,
  rel_tol = 1e-04,
  step_size_init = 0.1,
  max_iter = 500L,
  eta_param = 0.8,
  delta_param = 1e-04,
  rho = 0.5,
  bb_growth_cap = 4,
  coef_init = NULL,
  intercept_init = NULL
)
```

## Arguments

- X:

  Standardized design matrix.

- y:

  Response vector (length n).

- family:

  A pic family object (provides `name`).

- penalty:

  A pic penalty object (provides `name`, `params`).

- lambda_reg:

  Regularization parameter.

- fit_intercept:

  Whether to update an unpenalized intercept.

- rel_tol:

  Relative gradient-mapping tolerance.

- step_size_init:

  Initial step size.

- max_iter:

  Hard cap on outer iterations.

- eta_param:

  Memory weight of the non-monotone reference (in \[0,1)).

- delta_param:

  Sufficient-descent constant in the line search.

- rho:

  Step-size reduction factor (\< 1).

- bb_growth_cap:

  Maximum BB step growth over the previous step.

- coef_init, intercept_init:

  Optional warm starts.

## Value

A list with `coef`, `intercept`, `info`.

## Details

Cox is routed to a dedicated C++ FISTA (`fista_cox_cpp`) because its
loss has no intercept and uses a 2-column `(time, event)` response. The
build accepts only the six families and three penalties that have a C++
implementation; the family / penalty registries throw on any other input
upstream.
