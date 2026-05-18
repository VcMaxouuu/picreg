// Cox proportional hazards (Breslow approximation) — C++ port.
//
// Cox is treated as a special case in picr: its loss involves risk sets and
// tie-group metadata, so it does not fit the generic Family interface. The
// math lives here, in its own namespace, and is consumed by the dedicated
// `fista_cox_cpp` and `cox_null_grad_norms_cpp` entry points.
//
// The pivotal quantity for the PDB selector is family-specific too — a
// permutation-based estimator on the design matrix (see cox_null_grad.cpp).
//
// Design choice (vs the picpy approach): instead of caching tie metadata as
// mutable state inside a stateful family object and resetting it between
// fits, we compute it ONCE up-front in `make_data(times, events)` and pass
// the resulting `CoxData` struct by const reference. No mutable state, no
// `reset_fit_state()`.

#ifndef PICR_COX_H
#define PICR_COX_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

namespace picr::cox {

// Pre-computed survival data ready to be consumed by the loss/grad routines.
//   times          length n, sorted ascendingly (caller's responsibility).
//   events         length n, 0/1 indicator (1 = event, 0 = censored).
//   counts         length G = number of unique times; tie-group sizes.
//   starts         length G; 0-based start index of each tie group in
//                  the (sorted) (times, events) arrays.
//   sum_uncensored length G; sum of events within each tie group.
//   has_ties       true iff any tie group has size > 1.
//   n              total sample size.
struct CoxData {
  arma::vec   times;
  arma::vec   events;
  arma::uvec  counts;
  arma::uvec  starts;
  arma::vec   sum_uncensored;
  bool        has_ties;
  arma::uword n;
};

// Build a CoxData from a (times, events) pair. Assumes `times` is sorted
// ascendingly. Computes tie-group metadata in a single O(n) pass.
CoxData make_data(const arma::vec& times, const arma::vec& events);

// Stabilised partial log-likelihood per observation:
//   l = - (eta . events - sum_t d_t log sum_{j in R_t} exp(eta_j)) / n
//   evaluate = sqrt(l)        (phi = sqrt, like Gaussian)
double evaluate(const CoxData& cd, const arma::vec& lp);

// Residual r such that grad = X^T r. Equal to phi'(raw) * d(raw)/d(eta).
arma::vec residual(const CoxData& cd, const arma::vec& lp);

// Coefficient gradient. Cox has no intercept by construction; the returned
// scalar is 0 and provided only for API symmetry with the generic family
// system.
struct CoxGradResult {
  arma::vec coef;
  double    intercept;  // always 0
};
CoxGradResult grad(const arma::mat& X, const CoxData& cd, const arma::vec& lp);

}  // namespace picr::cox

#endif
