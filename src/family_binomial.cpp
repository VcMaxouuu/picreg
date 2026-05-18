// Binomial family — logistic link, identity stabiliser.
//
//
//   p            = g(lp) = 1 / (1 + exp(-lp))
//   raw          = 2 * mean( y * sqrt((1-p)/p) + (1-y) * sqrt(p/(1-p)) )
//   evaluate     = raw
//   residual     = g'(lp) * raw_loss_derivative
//                = p*(1-p) * (1/n) * (p - y) / (p*(1-p))^1.5
//                = (p - y) / (n * sqrt(p * (1 - p)))

#include "family_base.h"

namespace picr::binomial {

static arma::vec sigmoid(const arma::vec& lp) {
  return 1.0 / (1.0 + arma::exp(-lp));
}

double evaluate(const arma::vec& y, const arma::vec& lp) {
  arma::vec p = sigmoid(lp);
  p = arma::clamp(p, 1e-10, 1.0 - 1e-10);
  arma::vec ratio  = (1.0 - p) / p;                
  arma::vec sqrt_r = arma::sqrt(ratio);
  arma::vec term   = y % sqrt_r + (1.0 - y) / sqrt_r;
  return 2.0 * arma::mean(term);
}

arma::vec residual(const arma::vec& y, const arma::vec& lp) {
  arma::vec p = sigmoid(lp);
  p = arma::clamp(p, 1e-10, 1.0 - 1e-10);
  double n = (double) y.n_elem;
  return (p - y) / (n * arma::sqrt(p % (1.0 - p)));
}

GradResult grad(const arma::mat& X, const arma::vec& y, const arma::vec& lp) {
  arma::vec r = residual(y, lp);
  return { X.t() * r, arma::accu(r) };
}

double starting_intercept(const arma::vec& y) {
  double m = arma::mean(y);
  m = std::clamp(m, 1e-10, 1.0 - 1e-10);
  return std::log(m / (1.0 - m));                   // logit
}

arma::mat generate_y(int n, int n_simu) {
  return arma::conv_to<arma::mat>::from(
    arma::randu<arma::mat>(n, n_simu) < 0.5
  );
}

}  // namespace picr::binomial
