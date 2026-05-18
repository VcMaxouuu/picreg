// Gumbel family — identity link, exponential stabiliser.
//
//
//   g(lp)        = lp
//   z            = (y - lp) / sigma
//   raw          = log(sigma) + mean(z + exp(-z))
//   evaluate     = exp(raw)
//   residual     = phi'(raw) * g'(lp) * raw_loss_derivative
//                = exp(raw) * (exp(-z) - 1) / (sigma * n)

#include "family_base.h"

namespace picr::gumbel {

static constexpr double EULER_GAMMA = 0.57721566490153286060;

static arma::vec safe_exp(const arma::vec& x) {
  return arma::exp(arma::clamp(x, -50.0, 50.0));
}

static double safe_exp_scalar(double x) {
  return std::exp(std::clamp(x, -50.0, 50.0));
}

// Approximation used to estimate sigma from residuals.
// This is close to the method-of-moments estimate:
//
//   Var(Gumbel) = pi^2 sigma^2 / 6
//   sigma       = sqrt(6) / pi * sd(residuals)
//
// It is not exactly the same as scipy.stats.gumbel_r.fit(..., floc=0).
static double estimate_sigma(const arma::vec& y, const arma::vec& lp) {
  arma::vec residuals = y - lp;
  double sigma = std::sqrt(6.0) / arma::datum::pi * arma::stddev(residuals);
  return sigma;
}

static double estimate_sigma_MLE(const arma::vec& y, const arma::vec& lp) {
  arma::vec r = y - lp;
  double sigma = std::sqrt(6.0) / arma::datum::pi * arma::stddev(r);
  if (sigma <= 0) return 1.0;
  
  const int max_iter = 20;       // 4-5 suffisent en pratique
  const double tol = 1e-6;
  const double n = (double) r.n_elem;
  
  for (int iter = 0; iter < max_iter; ++iter) {
    arma::vec e = safe_exp(-r / sigma);
    double h       = -n * sigma + arma::accu(r % (1.0 - e));
    double h_prime = -n - arma::accu(r % r % e) / (sigma * sigma);
    double new_sigma = sigma - h / h_prime;
    if (new_sigma <= 0) new_sigma = sigma * 0.5;
    if (std::abs(new_sigma - sigma) < tol * sigma) { sigma = new_sigma; break; }
    sigma = new_sigma;
  }
  return sigma;
}


double evaluate(const arma::vec& y, const arma::vec& lp) {
  double sigma = estimate_sigma_MLE(y, lp);
  
  arma::vec z = (y - lp) / sigma;
  arma::vec exp_minus_z = safe_exp(-z);
  
  double raw = std::log(sigma) + arma::mean(z + exp_minus_z);
  return safe_exp_scalar(raw);
}

arma::vec residual(const arma::vec& y, const arma::vec& lp) {
  double sigma = estimate_sigma_MLE(y, lp);
  
  arma::vec z = (y - lp) / sigma;
  arma::vec exp_minus_z = safe_exp(-z);
  
  double raw = std::log(sigma) + arma::mean(z + exp_minus_z);
  double phi_prime = safe_exp_scalar(raw);
  
  double n = (double) y.n_elem;
  
  return phi_prime * (exp_minus_z - 1.0) / (sigma * n);
}

GradResult grad(const arma::mat& X, const arma::vec& y, const arma::vec& lp) {
  arma::vec r = residual(y, lp);
  return { X.t() * r, arma::accu(r) };
}

double starting_intercept(const arma::vec& y) {
  double sigma = std::sqrt(6.0) / arma::datum::pi * arma::stddev(y);
  return arma::mean(y) - sigma * EULER_GAMMA;
}

arma::mat generate_y(int n, int n_simu) {
  arma::mat U = arma::randu<arma::mat>(n, n_simu);
  U = arma::clamp(U, 1e-10, 1.0 - 1e-10);
  // If U ~ Uniform(0, 1), then -log(-log(U)) ~ Gumbel(0, 1).
  return -arma::log(-arma::log(U));
}

}  // namespace picr::gumbel