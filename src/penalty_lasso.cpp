// L1 (Lasso) penalty — convex, soft-thresholding prox.
//
//   P(beta)            = lambda * sum(|beta|)
//   prox(v, lambda, t) = sign(v) * max(|v| - t*lambda, 0)

#include "penalty_base.h"

namespace picr::lasso {

double evaluate(const arma::vec& beta, double lambda, double /*param*/) {
  return lambda * arma::accu(arma::abs(beta));
}

arma::vec prox(const arma::vec& v, double lambda, double step,
               double /*param*/) {
  double thr = lambda * step;
  return arma::sign(v) % arma::max(arma::abs(v) - thr, arma::zeros(v.n_elem));
}

}  // namespace picr::lasso
