// Gumbel family — identity link, exponential stabiliser.
//
//   z            = (y - lp) / sigma
//   raw          = log(sigma) + mean(z + exp(-z))
//   evaluate     = exp(raw)
//   residual     = exp(raw) * (exp(-z) - 1) / (sigma * n)
//
// The scale is profiled out at its conditional MLE for every lp. Its score is
// strictly decreasing, so a safeguarded Newton solve gives the unique positive
// root. At the null model, location and scale are estimated jointly.

#include "family_base.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace picr::gumbel {
namespace {

constexpr double ROOT_REL_TOL = 1e-12;
constexpr int ROOT_MAX_ITER = 64;
constexpr double LOG_DBL_MAX = 709.0;
constexpr double LOG_DBL_MIN = -745.0;

struct ScoreResult {
  double value;
  double derivative;
};

struct ProfileResult {
  double sigma;
  double raw_loss;
  double value;
  arma::vec residual;
};

struct LocationScaleResult {
  double intercept;
  double centered_intercept;
  double sigma;
};

// Solve a strictly decreasing score on (0, infinity). The positive root is
// first bracketed, then refined with Newton steps guarded by bisection.
template <typename ScoreFunction>
double solve_positive_root(double initial, ScoreFunction score) {
  if (!(initial > 0.0) || !std::isfinite(initial)) {
    Rcpp::stop("Gumbel MLE: could not initialize a positive scale.");
  }

  double lo = 0.0;
  double hi = initial;
  ScoreResult at_hi = score(hi);

  if (at_hi.value > 0.0) {
    lo = hi;
    for (int i = 0; i < ROOT_MAX_ITER && at_hi.value > 0.0; ++i) {
      if (hi > std::numeric_limits<double>::max() / 2.0) {
        Rcpp::stop("Gumbel MLE: scale overflow while bracketing the optimum.");
      }
      hi *= 2.0;
      at_hi = score(hi);
      if (at_hi.value > 0.0) lo = hi;
    }
    if (at_hi.value > 0.0) {
      Rcpp::stop("Gumbel MLE: failed to bracket the scale optimum.");
    }
  }

  double sigma = hi;
  for (int iter = 0; iter < ROOT_MAX_ITER; ++iter) {
    ScoreResult current = score(sigma);
    if (current.value == 0.0) return sigma;

    if (current.value > 0.0) {
      lo = sigma;
    } else {
      hi = sigma;
    }

    if (hi - lo <= ROOT_REL_TOL * hi) {
      return 0.5 * (lo + hi);
    }

    double candidate = sigma - current.value / current.derivative;
    if (!std::isfinite(candidate) || candidate <= lo || candidate >= hi) {
      candidate = 0.5 * (lo + hi);
    }
    sigma = candidate;
  }

  return 0.5 * (lo + hi);
}

double initial_scale(const arma::vec& values) {
  const double max_abs = arma::max(arma::abs(values));
  if (!(max_abs > 0.0) || !std::isfinite(max_abs)) {
    Rcpp::stop("Gumbel MLE is undefined when all residuals are zero.");
  }

  const double moment = std::sqrt(6.0) / arma::datum::pi
    * arma::stddev(values);
  return (moment > 0.0 && std::isfinite(moment)) ? moment : max_abs;
}

// Conditional scale MLE for fixed locations lp. Multiplying the usual scale
// score by -n*sigma^2 gives the strictly decreasing root function below.
double estimate_scale_mle(const arma::vec& r) {
  const double n = static_cast<double>(r.n_elem);
  const double initial = initial_scale(r);

  auto score = [&r, n](double sigma) {
    double h = -n * sigma;
    double h_prime = -n;

    for (arma::uword i = 0; i < r.n_elem; ++i) {
      const double scaled = r[i] / sigma;
      const double exponent = -scaled;

      // Overflow can only come from a negative residual, whose contribution
      // makes the decreasing score positive. Its sign is sufficient while
      // bracketing; no clipped surrogate objective is introduced.
      if (exponent > LOG_DBL_MAX) {
        return ScoreResult{std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity()};
      }

      const double e = exponent < LOG_DBL_MIN ? 0.0 : std::exp(exponent);
      h += r[i] * (1.0 - e);
      h_prime -= scaled * scaled * e;
    }
    return ScoreResult{h, h_prime};
  };

  return solve_positive_root(initial, score);
}

ProfileResult profile_from_residual(arma::vec r, double sigma,
                                    bool need_residual) {
  const double n = static_cast<double>(r.n_elem);
  arma::vec residual = need_residual ? arma::vec(r.n_elem) : arma::vec();
  double loss_sum = 0.0;

  for (arma::uword i = 0; i < r.n_elem; ++i) {
    const double z = r[i] / sigma;
    const double e = std::exp(-z);
    loss_sum += z + e;
    if (need_residual) residual[i] = e - 1.0;
  }

  const double raw_loss = std::log(sigma) + loss_sum / n;
  const double value = std::exp(raw_loss);
  if (need_residual) residual *= value / (sigma * n);

  return {sigma, raw_loss, value, std::move(residual)};
}

ProfileResult profile(const arma::vec& y, const arma::vec& lp,
                      bool need_residual) {
  arma::vec r = y - lp;
  const double sigma = estimate_scale_mle(r);
  return profile_from_residual(std::move(r), sigma, need_residual);
}

// Joint null MLE. For fixed sigma, the location score has the closed form
//   mu(sigma) = min(y) - sigma * log(mean(exp(-(y - min(y))/sigma))).
// Substitution into the scale score leaves one strictly decreasing equation.
LocationScaleResult estimate_null_mle(const arma::vec& y) {
  const double y_min = y.min();
  arma::vec d = y - y_min;
  const double max_d = d.max();
  if (!(max_d > 0.0) || !std::isfinite(max_d)) {
    Rcpp::stop("Gumbel location-scale MLE is undefined for a constant response.");
  }

  const double mean_d = arma::mean(d);
  const double initial = initial_scale(d);

  auto score = [&d, mean_d](double sigma) {
    double sum_w = 0.0;
    double sum_wd = 0.0;
    double sum_wd2 = 0.0;

    for (arma::uword i = 0; i < d.n_elem; ++i) {
      const double w = std::exp(-d[i] / sigma);
      sum_w += w;
      sum_wd += w * d[i];
      sum_wd2 += w * d[i] * d[i];
    }

    const double weighted_mean = sum_wd / sum_w;
    const double weighted_var = std::max(
      sum_wd2 / sum_w - weighted_mean * weighted_mean, 0.0);
    return ScoreResult{
      mean_d - weighted_mean - sigma,
      -1.0 - weighted_var / (sigma * sigma)
    };
  };

  const double sigma = solve_positive_root(initial, score);

  double sum_w = 0.0;
  for (arma::uword i = 0; i < d.n_elem; ++i) {
    sum_w += std::exp(-d[i] / sigma);
  }
  const double centered_intercept =
    -sigma * std::log(sum_w / static_cast<double>(d.n_elem));
  const double intercept = y_min + centered_intercept;

  return {intercept, centered_intercept, sigma};
}

}  // namespace

double evaluate(const arma::vec& y, const arma::vec& lp) {
  return profile(y, lp, false).value;
}

arma::vec residual(const arma::vec& y, const arma::vec& lp) {
  return profile(y, lp, true).residual;
}

GradResult grad(const arma::mat& X, const arma::vec& y, const arma::vec& lp) {
  arma::vec r = residual(y, lp);
  return {X.t() * r, arma::accu(r)};
}

ValueGradResult value_grad(const arma::mat& X, const arma::vec& y,
                           const arma::vec& lp) {
  ProfileResult out = profile(y, lp, true);
  const double intercept_grad = arma::accu(out.residual);
  arma::vec coef_grad = X.t() * out.residual;
  return {out.value, {std::move(coef_grad), intercept_grad}};
}

double starting_intercept(const arma::vec& y) {
  return estimate_null_mle(y).intercept;
}

NullResult null_fit(const arma::vec& y) {
  LocationScaleResult mle = estimate_null_mle(y);
  // Form residuals on the centered scale. Re-subtracting a large absolute
  // location would otherwise discard several digits when sigma is small.
  arma::vec r = (y - y.min()) - mle.centered_intercept;
  ProfileResult out = profile_from_residual(std::move(r), mle.sigma, true);
  return {mle.intercept, mle.sigma, std::move(out.residual)};
}

arma::mat generate_y(int n, int n_simu) {
  arma::mat U = arma::randu<arma::mat>(n, n_simu);
  U = arma::clamp(U, 1e-12, 1.0 - 1e-12);
  return -arma::log(-arma::log(U));
}

Rcpp::List profile_diagnostics(const arma::vec& y, const arma::vec& lp) {
  ProfileResult out = profile(y, lp, true);
  return Rcpp::List::create(
    Rcpp::Named("sigma") = out.sigma,
    Rcpp::Named("raw_loss") = out.raw_loss,
    Rcpp::Named("value") = out.value,
    Rcpp::Named("residual") = out.residual
  );
}

}  // namespace picr::gumbel

// Internal diagnostics shared with tests and R-side assessment.
// [[Rcpp::export]]
Rcpp::List gumbel_profile_cpp(const arma::vec& y, const arma::vec& lp) {
  return picr::gumbel::profile_diagnostics(y, lp);
}

// [[Rcpp::export]]
Rcpp::List gumbel_null_fit_cpp(const arma::vec& y) {
  picr::gumbel::NullResult out = picr::gumbel::null_fit(y);
  return Rcpp::List::create(
    Rcpp::Named("intercept") = out.intercept,
    Rcpp::Named("sigma") = out.sigma,
    Rcpp::Named("residual") = out.residual
  );
}
