// Gaussian family — identity link, sqrt stabiliser.
//
//
//   g(lp)        = lp
//   raw          = mean((y - lp)^2)
//   evaluate     = sqrt(raw)
//   residual     = phi'(raw) * g'(lp) * raw_loss_derivative
//                = (0.5 / sqrt(raw)) * 1 * (2/n) * (lp - y)
//                = (lp - y) / (n * sqrt(raw))

#include "family_base.h"

namespace picr::gaussian
{

    double evaluate(const arma::vec &y, const arma::vec &lp)
    {
        arma::vec d = y - lp;
        double raw = arma::mean(arma::square(d));
        return std::sqrt(raw);
    }

    arma::vec residual(const arma::vec &y, const arma::vec &lp)
    {
        arma::vec d = y - lp;
        double raw = arma::mean(arma::square(d));
        raw = std::max(raw, 1e-30);
        double n = (double)y.n_elem;
        return (-d) / (n * std::sqrt(raw));
    }

    GradResult grad(const arma::mat &X, const arma::vec &y, const arma::vec &lp)
    {
        arma::vec r = residual(y, lp);
        return {X.t() * r, arma::accu(r)};
    }

    double starting_intercept(const arma::vec &y)
    {
        return arma::mean(y);
    }

    arma::mat generate_y(int n, int n_simu)
    {
        return arma::randn<arma::mat>(n, n_simu);
    }

} // namespace picr::gaussian
