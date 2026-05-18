// picr — Penalty system (C++)
//
// One namespace per penalty. Same design as the family system: a `Penalty`
// struct holds function pointers populated once per fit by
// `get_penalty(name, ...)`. The FISTA solver then calls
// `pen.evaluate(...)` / `pen.prox(...)` in the hot loop without dispatch.
//
// Mathematical contract:
//
//   evaluate(beta, lambda) = lambda * P(beta)        (scalar penalty value)
//   prox(v, lambda, t)     = argmin_u { t * lambda * P(u) + 0.5 * ||u - v||^2 }
//   max_step               = strict upper bound on `t` for which `prox`
//                            is well-defined. Inf for convex penalties (L1),
//                            finite for non-convex (SCAD: a-1, MCP: gamma).
//
// SCAD/MCP have a tuning parameter (`a` / `gamma`) which is stored in the
// struct and passed through to the function pointer at every call.

#ifndef PICR_PENALTY_BASE_H
#define PICR_PENALTY_BASE_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>

namespace picr {

// Function-pointer table for a penalty. Built once per fit by
// `get_penalty(name, scad_a, mcp_gamma)`; the `param` slot carries the
// concavity parameter (a for SCAD, gamma for MCP, ignored for L1).
struct Penalty {
  const char* name;
  double      max_step;
  double      param;
  double      (*evaluate)(const arma::vec& beta, double lambda, double param);
  arma::vec   (*prox)    (const arma::vec& v,    double lambda, double step,
                          double param);
};

// ---------- per-penalty declarations --------------------------------------
//
// Each penalty lives in its own namespace and translation unit
// (penalty_<name>.cpp). The two public methods follow the contract above.

namespace lasso {
  double      evaluate(const arma::vec& beta, double lambda, double param);
  arma::vec   prox    (const arma::vec& v, double lambda, double step,
                       double param);
}

namespace scad {
  double      evaluate(const arma::vec& beta, double lambda, double param);
  arma::vec   prox    (const arma::vec& v, double lambda, double step,
                       double param);
}

namespace mcp {
  double      evaluate(const arma::vec& beta, double lambda, double param);
  arma::vec   prox    (const arma::vec& v, double lambda, double step,
                       double param);
}

// ---------- registry ------------------------------------------------------

// Build a Penalty handle from its canonical name.
// Recognised names: "lasso", "scad", "mcp" (lowercase).
// Throws via `Rcpp::stop` on any other input.
Penalty get_penalty(const std::string& name,
                    double scad_a    = 3.7,
                    double mcp_gamma = 3.0);

}  // namespace picr

#endif
