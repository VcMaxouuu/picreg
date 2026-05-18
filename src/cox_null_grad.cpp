// Cox PDB null-gradient — permutation-based pivotal estimator (optimised).
//
//
// Derivation:
//   Original loss gradient :  Xc[i,:] = Xp[i,:] - (1/(i+1)) * sum_{t<=i} Xp[t,:]
//              v       = Xc^T e
//   Expanding the cumulative mean and swapping summation order:
//              v = Xp^T (e - w),   w_t = sum_{i>=t} e_i / (i+1)   (suffix sum)
//   And since Xp = X[perm, :]:
//              Xp^T u = X^T u2,    u2[perm[i]] = u[i]   (scatter the vector)
//   => X is never permuted, never copied; no matrix cumsum, no transpose.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::export]]
Rcpp::List cox_null_grad_norms_cpp(const arma::mat& X,
                                   int n_simu,
                                   double alpha) {
  const uword n = X.n_rows;
  
  // Draw null events ~ Bernoulli(0.5), shape (n x n_simu).
  arma::mat events_null =
    arma::conv_to<arma::mat>::from(
      arma::randu<arma::mat>(n, n_simu) < 0.5
    );
  
  vec counts     = regspace<vec>(1.0, (double) n);
  vec log_counts = log(counts);
  vec inv_counts = 1.0 / counts;

  // Denominator per simulation:
  //   den_j = 2 * sqrt( sum_i e_ij * log_counts_i ) * sqrt(n) + eps
  rowvec sum_e_log = log_counts.t() * events_null;          // (1 x n_simu)
  rowvec den       = 2.0 * sqrt(sum_e_log) * std::sqrt((double) n) + 1e-10;
  
  vec stats_(n_simu);
  vec u2(n);
  
  for (int s = 0; s < n_simu; ++s) {
    vec e = events_null.col(s);
    
    // c_i = e_i / (i+1)
    vec c = e % inv_counts;
    
    // w = suffix sum of c, i.e. w_t = sum_{i>=t} c_i.
    // reverse -> cumsum -> reverse.
    vec w = flipud(cumsum(flipud(c)));
    
    vec u = e - w;                                  // length n, O(n)
    
    // Scatter u into u2 according to a fresh random permutation:
    //   u2[perm[i]] = u[i]   (equivalent to Xp^T u = X^T u2)
    uvec perm = randperm(n);
    u2.zeros();
    u2.elem(perm) = u;                              // O(n)
    
    // Single gemv on the fixed, contiguous matrix X.
    vec v = X.t() * u2;                             // O(n*p), BLAS-friendly
    
    stats_(s) = arma::abs(v).max();
  }
  
  stats_ = stats_ / den.t();
  double q = arma::as_scalar(
    arma::quantile(stats_, arma::vec{1.0 - alpha}));
  
  return Rcpp::List::create(
    Rcpp::Named("value")      = q,
    Rcpp::Named("statistics") = stats_
  );
}