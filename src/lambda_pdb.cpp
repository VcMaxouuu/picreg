// PDB selector — RcppArmadillo backend.
//
// Two entry points (called from R/lambda-pdb.R):
//   * pdb_mc_gaussian_cpp : pure linear-algebra Monte Carlo, family-agnostic
//   * pdb_mc_exact_cpp    : Monte Carlo under the family's null distribution
//                           — consumes the C++ family system in family_base.h
//
// Cox uses its dedicated permutation pivot; Gumbel is handled here at its
// joint null location/scale MLE.
//
// Both selectors process the n_simu draws in blocks of columns. Draws are
// generated in column-major order, i.e. in the same RNG order as a single
// n x n_simu matrix, so the statistics do not depend on the block size; only
// the peak memory does.

// [[Rcpp::depends(RcppArmadillo)]]
#include "family_base.h"

#include <algorithm>

using namespace arma;

// Number of draws per block: keeps each n x block matrix around 32 MB.
static int mc_block_size(int n, int n_simu)
{
    const double max_doubles = 4e6;
    const int block = (int)(max_doubles / std::max(n, 1));
    return std::max(1, std::min(n_simu, block));
}

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
    const int block = mc_block_size(n, n_simu);

    vec stats_(n_simu);
    for (int start = 0; start < n_simu; start += block)
    {
        const int b = std::min(block, n_simu - start);
        mat W = arma::randn<arma::mat>(n, b);
        mat XtW = X.t() * W; // (p x b)
        stats_.subvec(start, start + b - 1) = arma::max(arma::abs(XtW), 0).t();
    }

    stats_ *= std::sqrt(c_n) / (double)n;
    double q = arma::as_scalar(arma::quantile(stats_, arma::vec{1.0 - alpha}));

    return Rcpp::List::create(
        Rcpp::Named("value") = q,
        Rcpp::Named("statistics") = stats_);
}

// ---------------------------------------------------------------------------
// `mc_exact` Monte Carlo for the five non-Cox families.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List pdb_mc_exact_cpp(const arma::mat &X, std::string family_name,
                            int n_simu, double alpha)
{
    picr::Family fam = picr::get_family(family_name);

    const int n = (int)X.n_rows;
    const int block = mc_block_size(n, n_simu);
    const bool is_gumbel = (family_name == "gumbel");

    vec stats_(n_simu);
    for (int start = 0; start < n_simu; start += block)
    {
        const int b = std::min(block, n_simu - start);
        mat y_null = fam.generate_y(n, b);

        mat residual_mat(n, b);
        for (int j = 0; j < b; ++j)
        {
            if (is_gumbel)
            {
                // The exact null statistic is evaluated at the joint
                // location/scale MLE. null_fit() also reuses that solve to
                // form the residual.
                residual_mat.col(j) = picr::gumbel::null_fit(y_null.col(j)).residual;
            }
            else
            {
                vec yj = y_null.col(j);
                double b0 = fam.starting_intercept(yj);
                vec lp_j(n, fill::value(b0));
                residual_mat.col(j) = fam.residual(yj, lp_j);
            }
        }

        mat grad_mat = X.t() * residual_mat; // (p x b)
        stats_.subvec(start, start + b - 1) = arma::max(arma::abs(grad_mat), 0).t();
    }

    double q = arma::as_scalar(arma::quantile(stats_, arma::vec{1.0 - alpha}));

    return Rcpp::List::create(
        Rcpp::Named("value") = q,
        Rcpp::Named("statistics") = stats_
    );
}
