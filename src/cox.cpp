// Cox proportional hazards — math implementations.
//
// All formulas use the Breslow approximation for tied event times.

#include "cox.h"

namespace picr::cox {

// --------------------------------------------------------------------------
// Build CoxData from sorted (times, events) in a single O(n) pass.
// --------------------------------------------------------------------------
CoxData make_data(const arma::vec& times, const arma::vec& events) {
  const arma::uword n = times.n_elem;
  if (n == 0)
    Rcpp::stop("Cox: empty (times, events).");
  if (events.n_elem != n)
    Rcpp::stop("Cox: times and events must have the same length.");

  // Single pass to find tie-group boundaries.
  std::vector<arma::uword> counts_v;
  std::vector<arma::uword> starts_v;
  counts_v.reserve(n);
  starts_v.reserve(n);

  arma::uword start = 0;
  for (arma::uword i = 1; i <= n; ++i) {
    if (i == n || times(i) != times(start)) {
      starts_v.push_back(start);
      counts_v.push_back(i - start);
      start = i;
    }
  }

  const arma::uword G = counts_v.size();
  arma::uvec counts(G), starts(G);
  arma::vec  sum_uncens(G);
  for (arma::uword g = 0; g < G; ++g) {
    counts(g) = counts_v[g];
    starts(g) = starts_v[g];
    double s = 0.0;
    for (arma::uword k = 0; k < counts(g); ++k)
      s += events(starts(g) + k);
    sum_uncens(g) = s;
  }

  bool has_ties = arma::any(counts > 1);

  return CoxData{ times, events, counts, starts, sum_uncens, has_ties, n };
}

// --------------------------------------------------------------------------
// Shared workhorse: stabilised log-risk-sum vector, in a single pass.
//   risk_sum_i = sum_{j >= i} exp(lp_j)
//   log_risk_i = log(risk_sum_i)
// Returned in-place via outputs to avoid duplicate work between evaluate
// and residual.
// --------------------------------------------------------------------------
static void compute_risk_sums(const arma::vec& lp,
                              arma::vec& exp_shifted_out,
                              arma::vec& risk_sum_out,
                              arma::vec& log_risk_out,
                              double& m_out) {
  const double m = lp.max();
  m_out = m;
  exp_shifted_out = arma::exp(lp - m);
  // Reverse cumulative sum: risk_sum[i] = sum_{j >= i} exp_shifted[j]
  risk_sum_out = arma::flipud(arma::cumsum(arma::flipud(exp_shifted_out)));
  log_risk_out = arma::log(risk_sum_out) + m;
}

// Raw average partial log-likelihood (positive scalar).
static double raw_partial_loglik(const CoxData& cd,
                                 const arma::vec& lp,
                                 const arma::vec& log_risk) {
  double loss;
  if (!cd.has_ties) {
    loss = -arma::dot(lp, cd.events) + arma::dot(cd.events, log_risk);
  } else {
    arma::vec lr_per_time = log_risk.elem(cd.starts);
    loss = -arma::dot(lp, cd.events)
         + arma::dot(cd.sum_uncensored, lr_per_time);
  }
  return loss / (double) cd.n;
}

// --------------------------------------------------------------------------
// Public math.
// --------------------------------------------------------------------------

double evaluate(const CoxData& cd, const arma::vec& lp) {
  arma::vec exp_shifted, risk_sum, log_risk;
  double m;
  compute_risk_sums(lp, exp_shifted, risk_sum, log_risk, m);
  double raw = raw_partial_loglik(cd, lp, log_risk);
  return std::sqrt(std::max(raw, 0.0));   // phi = sqrt
}

arma::vec residual(const CoxData& cd, const arma::vec& lp) {
  arma::vec exp_shifted, risk_sum, log_risk;
  double m;
  compute_risk_sums(lp, exp_shifted, risk_sum, log_risk, m);

  // raw_loss_derivative
  arma::vec grad_eta(cd.n);
  if (cd.has_ties) {
    arma::vec rs_per_time = risk_sum.elem(cd.starts);
    arma::vec coeff     = cd.sum_uncensored / rs_per_time;
    arma::vec cum_coeff = arma::cumsum(coeff);
    // Repeat each cum_coeff[g] cd.counts[g] times to fill cd.n entries.
    arma::uword idx = 0;
    for (arma::uword g = 0; g < cd.counts.n_elem; ++g) {
      for (arma::uword k = 0; k < cd.counts(g); ++k) {
        grad_eta(idx) = (-cd.events(idx) +
                         exp_shifted(idx) * cum_coeff(g)) / (double) cd.n;
        ++idx;
      }
    }
  } else {
    arma::vec coeff     = cd.events / risk_sum;
    arma::vec cum_coeff = arma::cumsum(coeff);
    grad_eta = (-cd.events + exp_shifted % cum_coeff) / (double) cd.n;
  }

  // phi = sqrt → phi'(raw) = 0.5 / sqrt(raw)
  double raw = raw_partial_loglik(cd, lp, log_risk);
  raw = std::max(raw, 1e-30);
  double phi_prime = 0.5 / std::sqrt(raw);
  return phi_prime * grad_eta;
}

CoxGradResult grad(const arma::mat& X, const CoxData& cd, const arma::vec& lp) {
  arma::vec r = residual(cd, lp);
  return { X.t() * r, 0.0 };
}

}  // namespace picr::cox
