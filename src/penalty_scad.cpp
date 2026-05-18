// SCAD (Smoothly Clipped Absolute Deviation) — non-convex.
// Fan & Li (2001).
//
// Penalty value:
//   |b| <= lambda       : lambda * |b|
//   lambda < |b| <= a*l : (2*a*lambda*|b| - b^2 - lambda^2) / (2*(a-1))
//   |b| > a*lambda      : (a+1) * lambda^2 / 2
//
// Proximal map (closed-form, piecewise):
//   |v| <= lambda*(1+t)       : soft_threshold(v, t*lambda)
//   lambda*(1+t) < |v| <= a*l : ((a-1)/(a-1-t)) * soft_threshold(v, a*lambda*t/(a-1))
//   |v| > a*lambda            : v (no shrinkage)
//
// Step-size domain: t < a - 1 (otherwise the middle region blows up).

#include "penalty_base.h"

namespace picr::scad {

double evaluate(const arma::vec& beta, double lambda, double a) {
  arma::vec abs_b = arma::abs(beta);
  arma::vec p(abs_b.n_elem, arma::fill::zeros);
  double a_minus_1 = a - 1.0;
  double const_high = (a + 1.0) * lambda * lambda / 2.0;

  for (arma::uword i = 0; i < abs_b.n_elem; ++i) {
    double b = abs_b(i);
    if (b <= lambda) {
      p(i) = lambda * b;
    } else if (b <= a * lambda) {
      p(i) = (2.0 * a * lambda * b - b * b - lambda * lambda) /
             (2.0 * a_minus_1);
    } else {
      p(i) = const_high;
    }
  }
  return arma::accu(p);
}

arma::vec prox(const arma::vec& v, double lambda, double step, double a) {
  arma::vec abs_v = arma::abs(v);
  arma::vec out(v.n_elem, arma::fill::zeros);
  double a_minus_1 = a - 1.0;

  // Caller must enforce step < max_step = a - 1; we still guard against
  // numerical edges.
  double denom_mid = a_minus_1 - step;

  for (arma::uword i = 0; i < v.n_elem; ++i) {
    double absv = abs_v(i);
    double s = (v(i) >= 0) ? 1.0 : -1.0;

    if (absv <= lambda * (1.0 + step)) {
      double shrunk = absv - lambda * step;
      out(i) = (shrunk > 0) ? s * shrunk : 0.0;
    } else if (absv <= a * lambda) {
      double inner_thr = (a * lambda * step) / a_minus_1;
      double shrunk = absv - inner_thr;
      double scaled = (a_minus_1 / denom_mid) *
                      ((shrunk > 0) ? s * shrunk : 0.0);
      out(i) = scaled;
    } else {
      out(i) = v(i);
    }
  }
  return out;
}

}  // namespace picr::scad
