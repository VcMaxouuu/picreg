// String → Penalty dispatcher.
//
// Single point of truth for the list of available penalties. Called once
// per fit at the entry of `fista_cpp`; the hot loop uses the returned
// struct directly without further dispatch.
//
// `param` carries the concavity parameter:
//   - lasso : unused
//   - scad  : a       (also defines max_step = a - 1)
//   - mcp   : gamma   (also defines max_step = gamma)
//
// To add a new penalty:
//   1. Create `penalty_<name>.cpp` defining `picr::<name>::*`.
//   2. Forward-declare the namespace in `penalty_base.h`.
//   3. Add a matching `if` branch below.

#include "penalty_base.h"

namespace picr {

Penalty get_penalty(const std::string& name,
                    double scad_a,
                    double mcp_gamma) {
  if (name == "lasso") {
    return Penalty{
      "lasso",
      std::numeric_limits<double>::infinity(),  // max_step = Inf
      0.0,                                       // param unused
      &lasso::evaluate,
      &lasso::prox
    };
  }
  if (name == "scad") {
    if (scad_a <= 2.0)
      Rcpp::stop("SCAD requires a > 2; got a = " + std::to_string(scad_a));
    return Penalty{
      "scad",
      scad_a - 1.0,                              // max_step = a - 1
      scad_a,
      &scad::evaluate,
      &scad::prox
    };
  }
  if (name == "mcp") {
    if (mcp_gamma <= 1.0)
      Rcpp::stop("MCP requires gamma > 1; got gamma = " +
                 std::to_string(mcp_gamma));
    return Penalty{
      "mcp",
      mcp_gamma,                                 // max_step = gamma
      mcp_gamma,
      &mcp::evaluate,
      &mcp::prox
    };
  }
  Rcpp::stop("Unknown penalty in C++ registry: " + name);
}

}  // namespace picr
