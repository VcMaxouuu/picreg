// lambda_max — closed-form KKT threshold for the warm-start path entry.
//
// For the convex relaxation min_beta L(beta) + lambda * ||beta||_1, the KKT
// conditions give: beta = 0 (with the intercept at its unpenalised optimum
// beta0*) is a minimiser iff
//     ||grad_beta L(beta = 0, beta0*)||_infty <= lambda.
// The smallest such lambda — i.e. the smallest lambda that kills every
// coordinate — is therefore
//     lambda_max = ||grad_beta L(beta = 0, beta0*)||_infty.
//
// SCAD and MCP share this threshold because their derivative at 0 equals
// lambda (lasso-like behaviour near the origin).
//
// This is structurally the same statistic the PDB selector evaluates on
// simulated null draws (see lambda_pdb.cpp / cox_null_grad.cpp), except
// here it is evaluated on the observed y at lp = beta0* (or lp = 0 for
// Cox, which has no intercept).

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "family_base.h"
#include "cox.h"

// [[Rcpp::export]]
double lambda_max_cpp(const arma::mat& X,
                      const arma::vec& y,
                      std::string family_name,
                      bool fit_intercept) {
  picr::Family fam = picr::get_family(family_name);
  const double b0 = fit_intercept ? fam.starting_intercept(y) : 0.0;
  arma::vec lp(X.n_rows, arma::fill::value(b0));
  picr::GradResult g = fam.grad(X, y, lp);
  return arma::max(arma::abs(g.coef));
}

// [[Rcpp::export]]
double lambda_max_cox_cpp(const arma::mat& X,
                          const arma::vec& times,
                          const arma::vec& events) {
  picr::cox::CoxData cd = picr::cox::make_data(times, events);
  arma::vec lp(X.n_rows, arma::fill::zeros);
  picr::cox::CoxGradResult g = picr::cox::grad(X, cd, lp);
  return arma::max(arma::abs(g.coef));
}
