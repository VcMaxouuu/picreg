// Exponential family — log link, identity stabiliser.
//
//
//   mu           = g(lp) = exp(lp)
//   raw          = mean( log(mu) + y / mu )
//   evaluate     = raw
//   residual     = phi'(raw) * g'(lp) * raw_loss_derivative
//                = mu * (1/n) * (1/mu - y/mu^2)
//                = (1 - y / mu) / n

#include "family_base.h"

namespace picr::exponential {

double evaluate(const arma::vec& y, const arma::vec& lp) {
  arma::vec mu = arma::exp(lp);
  mu = arma::clamp(mu, 1e-10, arma::datum::inf);
  return arma::mean(arma::log(mu) + y / mu);
}

arma::vec residual(const arma::vec& y, const arma::vec& lp) {
  arma::vec mu = arma::exp(lp);
  mu = arma::clamp(mu, 1e-10, arma::datum::inf);
  double n = (double) y.n_elem;
  return (1.0 - y / mu) / n;
}

GradResult grad(const arma::mat& X, const arma::vec& y, const arma::vec& lp) {
  arma::vec r = residual(y, lp);
  return { X.t() * r, arma::accu(r) };
}

double starting_intercept(const arma::vec& y) {
  double m = arma::mean(y);
  m = std::max(m, 1e-10);
  return std::log(m);
}

arma::mat generate_y(int n, int n_simu) {
  arma::mat U = arma::randu<arma::mat>(n, n_simu);
  U = arma::clamp(U, 1e-12, 1.0); 
  return - 5.0 * arma::log(U);
}

}  // namespace picr::exponential
