// PDB selector — RcppArmadillo backend.
//
// Two entry points (called from R/lambda-pdb.R):
//   * pdb_mc_gaussian_cpp : pure linear-algebra Monte Carlo, family-agnostic
//   * pdb_mc_exact_cpp    : Monte Carlo under the family's null distribution
//                           — consumes the C++ family system in family_base.h
//
// Custom-pivot families (Gumbel, Cox) stay in R via `family$null_grad_norms`.

// [[Rcpp::depends(RcppArmadillo)]]
#include "family_base.h"

using namespace arma;

// ---------------------------------------------------------------------------
// `mc_gaussian` Monte Carlo.
//
// Draws W ~ N(0, I_n)^{n_simu} and returns the (1 - alpha) empirical
// quantile of (sqrt(c_n) / n) * ||X^T W||_inf along columns.
// Family-agnostic by design: c_n is supplied by the caller.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List pdb_mc_gaussian_cpp(const arma::mat &X, double c_n,
                               int n_simu, double alpha)
{
    const int n = (int)X.n_rows;

    mat W = arma::randn<arma::mat>(n, n_simu);

    mat XtW = X.t() * W; // (p x n_simu)
    vec stats_ = arma::max(arma::abs(XtW), 0).t();

    stats_ *= std::sqrt(c_n) / (double)n;
    double q = arma::as_scalar(arma::quantile(stats_, arma::vec{1.0 - alpha}));

    return Rcpp::List::create(
        Rcpp::Named("value") = q,
        Rcpp::Named("statistics") = stats_);
}

// ---------------------------------------------------------------------------
// `mc_exact` Monte Carlo for the four standard families.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List pdb_mc_exact_cpp(const arma::mat &X, std::string family_name,
                            int n_simu, double alpha)
{
    picr::Family fam = picr::get_family(family_name);

    const int n = (int)X.n_rows;
    mat y_null = fam.generate_y(n, n_simu);

    mat residual_mat(n, n_simu);
    for (int j = 0; j < n_simu; ++j)
    {
        vec yj = y_null.col(j);
        double b0 = fam.starting_intercept(yj);
        vec lp_j(n, fill::value(b0));
        residual_mat.col(j) = fam.residual(yj, lp_j);
    }

    mat grad_mat = X.t() * residual_mat; // (p x n_simu)
    vec stats_ = arma::max(arma::abs(grad_mat), 0).t();

    double q = arma::as_scalar(arma::quantile(stats_, arma::vec{1.0 - alpha}));

    return Rcpp::List::create(
        Rcpp::Named("value") = q,
        Rcpp::Named("statistics") = stats_
    );
}
