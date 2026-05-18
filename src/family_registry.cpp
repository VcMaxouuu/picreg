// String → Family dispatcher.
//
// Single point of truth for the list of available families. Called once
// per fit at the entry of `lambda_pdb` / `fista`; the hot loop then uses
// the returned struct directly without any further dispatch.
//
// To add a new family:
//   1. Create `family_<name>.cpp` defining `picr::<name>::*`.
//   2. Forward-declare the namespace in `family_base.h`.
//   3. Add a matching `if` branch below.
// No other file needs to know about the new family.

#include "family_base.h"

namespace picr {

// Build a Family handle for `name`. Each branch wires the public function
// pointers to the corresponding namespace's symbols. Unknown names abort
// via `Rcpp::stop` so an R-side typo surfaces cleanly.
Family get_family(const std::string& name) {
  if (name == "gaussian") {
    return Family{ "gaussian",
                   &gaussian::evaluate, &gaussian::residual, &gaussian::grad,
                   &gaussian::starting_intercept, &gaussian::generate_y };
  }
  if (name == "binomial") {
    return Family{ "binomial",
                   &binomial::evaluate, &binomial::residual, &binomial::grad,
                   &binomial::starting_intercept, &binomial::generate_y };
  }
  if (name == "poisson") {
    return Family{ "poisson",
                   &poisson::evaluate, &poisson::residual, &poisson::grad,
                   &poisson::starting_intercept, &poisson::generate_y };
  }
  if (name == "exponential") {
    return Family{ "exponential",
                   &exponential::evaluate, &exponential::residual,
                   &exponential::grad,
                   &exponential::starting_intercept, &exponential::generate_y };
  }
  if (name == "gumbel") {
    return Family{ "gumbel",
                   &gumbel::evaluate, &gumbel::residual, &gumbel::grad,
                   &gumbel::starting_intercept, &gumbel::generate_y };
  }
  Rcpp::stop("Unknown family in C++ registry: " + name);
}

}  // namespace picr
