// MCP (Minimax Concave Penalty) — non-convex. Zhang (2010).
//
// Penalty value:
//   |b| <= gamma*lambda : lambda*|b| - b^2 / (2*gamma)
//   |b| >  gamma*lambda : gamma * lambda^2 / 2
//
// Proximal map (closed-form):
//   |v| <= gamma*lambda : (gamma/(gamma-t)) * soft_threshold(v, t*lambda)
//   |v| >  gamma*lambda : v (no shrinkage)
//
// Step-size domain: t < gamma (otherwise the rescaling blows up).

#include "penalty_base.h"

namespace picr::mcp {

double evaluate(const arma::vec& beta, double lambda, double gamma) {
  arma::vec abs_b = arma::abs(beta);
  arma::vec p(abs_b.n_elem, arma::fill::zeros);
  double const_high = gamma * lambda * lambda / 2.0;
  double thresh     = gamma * lambda;

  for (arma::uword i = 0; i < abs_b.n_elem; ++i) {
    double b = abs_b(i);
    if (b <= thresh) {
      p(i) = lambda * b - (b * b) / (2.0 * gamma);
    } else {
      p(i) = const_high;
    }
  }
  return arma::accu(p);
}

arma::vec prox(const arma::vec& v, double lambda, double step, double gamma) {
  arma::vec abs_v = arma::abs(v);
  arma::vec out(v.n_elem, arma::fill::zeros);
  double thresh = gamma * lambda;
  double denom  = gamma - step;          // caller enforces step < gamma
  double scale  = gamma / denom;

  for (arma::uword i = 0; i < v.n_elem; ++i) {
    double absv = abs_v(i);
    if (absv <= thresh) {
      double s = (v(i) >= 0) ? 1.0 : -1.0;
      double shrunk = absv - lambda * step;
      out(i) = (shrunk > 0) ? scale * s * shrunk : 0.0;
    } else {
      out(i) = v(i);
    }
  }
  return out;
}

}  // namespace picr::mcp
