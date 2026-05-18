// Poisson family — log link, identity stabiliser.
//
//
//   mu           = g(lp) = lp
//   raw          = 2 * mean( y / sqrt(mu) + sqrt(mu) )
//   evaluate     = raw
//   residual     = phi'(raw) * g'(lp) * raw_loss_derivative
//                = mu * (1/n) * (mu - y) / mu^1.5
//                = (lp - y) / (n * sqrt(raw))

#include "family_base.h"

namespace picr::poisson {

double evaluate(const arma::vec& y, const arma::vec& lp) {
  arma::vec mu = arma::exp(lp);
  mu = arma::clamp(mu, 1e-10, arma::datum::inf);
  return 2.0 * arma::mean(y / arma::sqrt(mu) + arma::sqrt(mu));
}

arma::vec residual(const arma::vec& y, const arma::vec& lp) {
  arma::vec mu = arma::exp(lp);
  mu = arma::clamp(mu, 1e-10, arma::datum::inf);
  double n = (double) y.n_elem;
  return (mu - y) / (n * arma::sqrt(mu));
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
  arma::mat Y(n, n_simu);
  Y.imbue([]() {
    return (double) R::rpois(5.0);
  });
  return Y;
}

}  // namespace picr::poisson
