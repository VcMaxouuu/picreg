// picr — Family system (C++)
//
// One namespace per family. The `Family` struct is a 
// function-pointer table populated once by `get_family(name)` and then
// consumed without further dispatch by `lambda_pdb` and `fista`.
//
// Mathematical contract:
//
//   y_pred              = g(lp)                     [mean-function link]
//   raw                 = raw_loss(y, y_pred)       [non-transformed loss]
//   evaluate(y, lp)     = phi(raw)                  [transformed loss]
//   residual(y, lp)     = phi'(raw) * g'(lp) * d/dy_pred raw_loss(y, y_pred)
//   grad(X, y, lp).coef = X^T residual
//   grad(X, y, lp).int  = sum(residual)

#ifndef PICR_FAMILY_BASE_H
#define PICR_FAMILY_BASE_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>

namespace picr {

// Joint gradient of the stabilised loss at a given linear predictor.
//   coef:      gradient w.r.t. the coefficient vector beta, length p.
//   intercept: gradient w.r.t. the unpenalised intercept, scalar.
struct GradResult {
  arma::vec coef;
  double    intercept;
};

// Family handle: a table of function pointers, populated once per fit by
// `get_family(name)` and then consumed by `lambda_pdb` / `fista` without
// further dispatch.
//
// Slots:
//   name              human-readable identifier ("gaussian", ...).
//   evaluate          stabilised loss   phi(raw_loss(y, g(lp))).
//   residual          element-wise residual r used internally for `grad`;
//                     X^T r yields the coefficient gradient.
//   grad              full (coef, intercept) gradient at (X, y, lp).
//   starting_intercept g^{-1}(mean(y)); the canonical intercept warm start.
//   generate_y        null draws (n x n_simu) for Monte Carlo selectors.
struct Family {
  const char* name;
  double      (*evaluate)         (const arma::vec&, const arma::vec&);
  arma::vec   (*residual)         (const arma::vec&, const arma::vec&);
  GradResult  (*grad)             (const arma::mat&, const arma::vec&,
                                   const arma::vec&);
  double      (*starting_intercept)(const arma::vec&);
  arma::mat   (*generate_y)       (int n, int n_simu);
};

// ---------- per-family declarations ---------------------------------------
//
// Each family lives in its own namespace and its own translation unit
// (family_<name>.cpp). The five public methods follow the contract above;
// see the corresponding .cpp for the closed-form derivations.

namespace gaussian {
  double      evaluate         (const arma::vec& y, const arma::vec& lp);
  arma::vec   residual         (const arma::vec& y, const arma::vec& lp);
  GradResult  grad             (const arma::mat& X, const arma::vec& y,
                                const arma::vec& lp);
  double      starting_intercept(const arma::vec& y);
  arma::mat   generate_y       (int n, int n_simu);
}

namespace binomial {
  double      evaluate         (const arma::vec& y, const arma::vec& lp);
  arma::vec   residual         (const arma::vec& y, const arma::vec& lp);
  GradResult  grad             (const arma::mat& X, const arma::vec& y,
                                const arma::vec& lp);
  double      starting_intercept(const arma::vec& y);
  arma::mat   generate_y       (int n, int n_simu);
}

namespace poisson {
  double      evaluate         (const arma::vec& y, const arma::vec& lp);
  arma::vec   residual         (const arma::vec& y, const arma::vec& lp);
  GradResult  grad             (const arma::mat& X, const arma::vec& y,
                                const arma::vec& lp);
  double      starting_intercept(const arma::vec& y);
  arma::mat   generate_y       (int n, int n_simu);
}

namespace exponential {
  double      evaluate         (const arma::vec& y, const arma::vec& lp);
  arma::vec   residual         (const arma::vec& y, const arma::vec& lp);
  GradResult  grad             (const arma::mat& X, const arma::vec& y,
                                const arma::vec& lp);
  double      starting_intercept(const arma::vec& y);
  arma::mat   generate_y       (int n, int n_simu);
}

namespace gumbel {
  double      evaluate         (const arma::vec& y, const arma::vec& lp);
  arma::vec   residual         (const arma::vec& y, const arma::vec& lp);
  GradResult  grad             (const arma::mat& X, const arma::vec& y,
                                const arma::vec& lp);
  double      starting_intercept(const arma::vec& y);
  arma::mat   generate_y       (int n, int n_simu);
}

// ---------- registry ------------------------------------------------------

// Build a Family handle from its canonical name. The returned struct holds
// pointers to the corresponding namespace's functions and can be used
// freely afterwards without further dispatch.
//
// Recognised names: "gaussian", "binomial", "poisson", "exponential",
// "gumbel". Cox is handled separately (see src/cox.h).
// Throws via `Rcpp::stop` on any other input.
Family get_family(const std::string& name);

}  // namespace picr

#endif
