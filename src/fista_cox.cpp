// FISTA dedicated to Cox proportional hazards.
//
// Mirrors the algorithm of `fista.cpp` (nmAPG, Li & Lin 2015) but
// specialised for Cox:
//   - y is replaced by a CoxData struct (pre-computed tie metadata)
//   - no intercept (Cox is identifiable only up to additive shifts in lp)
//   - calls picr::cox::{evaluate, residual, grad} for the Cox math
//
// The penalty layer is shared with the generic FISTA.

// [[Rcpp::depends(RcppArmadillo)]]
#include "cox.h"
#include "penalty_base.h"

using namespace arma;

// ---------- helpers --------------------------------------------------------
// Coefficient-only BB step (Cox has no intercept).
static double bb_step_coef(const vec& coef_a, const vec& coef_b,
                           const vec& grad_a, const vec& grad_b,
                           double fallback,
                           double growth_cap, double alpha_pen_cap,
                           double alpha_min = 1e-10, double alpha_max = 1e10,
                           double eps = 1e-12) {
  vec s = coef_a - coef_b;
  vec r = grad_a - grad_b;
  double ss = dot(s, s);
  double sr = dot(s, r);
  if (ss <= eps || sr <= eps) {
    return std::min(fallback, alpha_pen_cap);
  }
  double upper = std::min({alpha_max, growth_cap * fallback, alpha_pen_cap});
  return std::max(alpha_min, std::min(ss / sr, upper));
}

// Line search result from y_k.
namespace {
struct LsResult {
  vec    z_coef;
  double F_z;
  double alpha;
  bool   non_monotone;
};
}  // namespace

static LsResult ls_from_y(const mat& X, const picr::cox::CoxData& cd,
                          const picr::Penalty& pen, double lambda_reg,
                          const vec& y_coef, double f_y,
                          const vec& grad_y,
                          double alpha, double c_k,
                          double delta_param, double rho,
                          int max_inner = 30) {
  vec z_coef = y_coef;
  double f_z = f_y;

  for (int i = 0; i < max_inner; ++i) {
    z_coef = pen.prox(y_coef - alpha * grad_y, lambda_reg, alpha, pen.param);

    vec lp = X * z_coef;
    f_z = picr::cox::evaluate(cd, lp);

    vec    d        = z_coef - y_coef;
    double inner_p  = dot(d, grad_y);
    double diff_sq  = dot(d, d);

    if (f_z <= f_y + inner_p + diff_sq / (2.0 * alpha)) {
      double F_z_val = f_z + pen.evaluate(z_coef, lambda_reg, pen.param);
      bool nm = (F_z_val <= c_k - delta_param * diff_sq);
      return { z_coef, F_z_val, alpha, nm };
    }
    alpha *= rho;
  }
  double F_z_val = f_z + pen.evaluate(z_coef, lambda_reg, pen.param);
  return { z_coef, F_z_val, alpha, false };
}

// v-branch fallback from x_k.
namespace {
struct VResult {
  vec    coef;
  double F_val;
  double alpha;
  char   picked;        // 'z' or 'v'
};
}  // namespace

static VResult v_branch(const mat& X, const picr::cox::CoxData& cd,
                        const picr::Penalty& pen, double lambda_reg,
                        const vec& x_coef,
                        const vec& y_prev_coef,
                        const vec& grad_y_prev,
                        const vec& z_next_coef, double F_z,
                        double alpha_x, double rho,
                        double bb_growth_cap, double alpha_pen_cap,
                        int max_inner = 30) {
  vec lp_x = X * x_coef;
  double f_x = picr::cox::evaluate(cd, lp_x);
  picr::cox::CoxGradResult gx = picr::cox::grad(X, cd, lp_x);

  alpha_x = bb_step_coef(x_coef, y_prev_coef, gx.coef, grad_y_prev,
                         alpha_x, bb_growth_cap, alpha_pen_cap);
  alpha_x = std::min(alpha_x, alpha_pen_cap);

  vec    v_coef = x_coef;
  double f_v    = f_x;
  double F_v    = std::numeric_limits<double>::infinity();

  for (int i = 0; i < max_inner; ++i) {
    v_coef = pen.prox(x_coef - alpha_x * gx.coef, lambda_reg, alpha_x,
                      pen.param);
    vec lp = X * v_coef;
    f_v = picr::cox::evaluate(cd, lp);

    vec    d       = v_coef - x_coef;
    double inner_p = dot(d, gx.coef);
    double diff_sq = dot(d, d);

    if (f_v <= f_x + inner_p + diff_sq / (2.0 * alpha_x)) {
      F_v = f_v + pen.evaluate(v_coef, lambda_reg, pen.param);
      break;
    }
    alpha_x *= rho;
  }

  if (!std::isfinite(F_v)) {
    F_v = f_v + pen.evaluate(v_coef, lambda_reg, pen.param);
  }

  if (F_z <= F_v) {
    return { z_next_coef, F_z, alpha_x, 'z' };
  }
  return { v_coef, F_v, alpha_x, 'v' };
}

// ---------- entry point ---------------------------------------------------

// [[Rcpp::export]]
Rcpp::List fista_cox_cpp(const arma::mat& X,
                          const arma::vec& times,
                          const arma::vec& events,
                          std::string penalty_name,
                          double lambda_reg,
                          double scad_a       = 3.7,
                          double mcp_gamma    = 3.0,
                          double rel_tol      = 1e-4,
                          double step_size_init = 1e-2,
                          int    max_iter     = 500,
                          double eta_param    = 0.8,
                          double delta_param  = 1e-4,
                          double rho          = 0.5,
                          double bb_growth_cap = 2.0,
                          Rcpp::Nullable<Rcpp::NumericVector> coef_init = R_NilValue) {

  picr::cox::CoxData cd = picr::cox::make_data(times, events);
  picr::Penalty     pen = picr::get_penalty(penalty_name, scad_a, mcp_gamma);

  const double step_margin = 0.95;
  const double alpha_pen_cap = step_margin * pen.max_step;
  step_size_init = std::min(step_size_init, alpha_pen_cap);

  const int p = X.n_cols;

  vec x_coef = (coef_init.isNull())
               ? vec(p, fill::zeros)
               : Rcpp::as<vec>(coef_init);

  vec x_prev_coef = x_coef;
  vec z_coef      = x_coef;
  vec y_prev_coef = x_coef;

  vec lp = X * x_coef;
  picr::cox::CoxGradResult g0 = picr::cox::grad(X, cd, lp);
  vec grad_y_prev = g0.coef;

  double t_prev = 0.0, t_curr = 1.0;
  double F_x = picr::cox::evaluate(cd, lp) + pen.evaluate(x_coef, lambda_reg, pen.param);
  double c_k = F_x, q_k = 1.0;

  double alpha_y = step_size_init;
  double alpha_x = step_size_init;

  std::vector<double>           objective_path;
  std::vector<std::vector<int>> features_path;
  std::vector<double>           gradient_mapping_path;
  int    n_iter = 0;
  std::string status = "max_iter_reached";

  for (int k = 1; k <= max_iter; ++k) {
    double a = t_prev / t_curr;
    double b = (t_prev - 1.0) / t_curr;

    vec y_coef = x_coef + a * (z_coef - x_coef) + b * (x_coef - x_prev_coef);

    vec lp_y = X * y_coef;
    double f_y = picr::cox::evaluate(cd, lp_y);
    picr::cox::CoxGradResult gy = picr::cox::grad(X, cd, lp_y);

    alpha_y = bb_step_coef(y_coef, y_prev_coef, gy.coef, grad_y_prev,
                           alpha_y, bb_growth_cap, alpha_pen_cap);

    LsResult ls_y = ls_from_y(X, cd, pen, lambda_reg,
                              y_coef, f_y, gy.coef,
                              std::min(alpha_y, alpha_pen_cap),
                              c_k, delta_param, rho);
    alpha_y = ls_y.alpha;

    vec    x_next_coef;
    double F_next;
    vec    gm_anchor_coef;
    double gm_alpha;

    if (ls_y.non_monotone) {
      x_next_coef    = ls_y.z_coef;
      F_next         = ls_y.F_z;
      gm_anchor_coef = y_coef;
      gm_alpha       = alpha_y;
    } else {
      VResult vb = v_branch(X, cd, pen, lambda_reg,
                            x_coef, y_prev_coef, grad_y_prev,
                            ls_y.z_coef, ls_y.F_z,
                            alpha_x, rho, bb_growth_cap, alpha_pen_cap);
      x_next_coef = vb.coef;
      F_next      = vb.F_val;
      alpha_x     = vb.alpha;
      if (vb.picked == 'z') {
        gm_anchor_coef = y_coef;
        gm_alpha       = alpha_y;
      } else {
        gm_anchor_coef = x_coef;
        gm_alpha       = alpha_x;
      }
    }

    objective_path.push_back(F_next);
    {
      std::vector<int> nz;
      nz.reserve(x_next_coef.n_elem);
      for (uword i = 0; i < x_next_coef.n_elem; ++i)
        if (x_next_coef(i) != 0.0) nz.push_back((int) i + 1);
      features_path.push_back(nz);
    }
    n_iter = k;

    vec    G_coef    = (gm_anchor_coef - x_next_coef) / gm_alpha;
    double G_norm    = std::sqrt(dot(G_coef, G_coef));
    double beta_norm = std::sqrt(dot(gm_anchor_coef, gm_anchor_coef));
    double rel_gm    = G_norm / std::max(beta_norm, 1.0);
    gradient_mapping_path.push_back(rel_gm);

    if (rel_gm < rel_tol) {
      x_coef = x_next_coef;
      status = "converged";
      break;
    }

    x_prev_coef = x_coef;
    z_coef      = ls_y.z_coef;
    y_prev_coef = y_coef;
    grad_y_prev = gy.coef;
    x_coef      = x_next_coef;

    double t_next = (std::sqrt(4 * t_curr * t_curr + 1.0) + 1.0) / 2.0;
    t_prev = t_curr;
    t_curr = t_next;

    double q_next = eta_param * q_k + 1.0;
    c_k = (eta_param * q_k * c_k + F_next) / q_next;
    q_k = q_next;
  }

  Rcpp::List info = Rcpp::List::create(
    Rcpp::Named("n_iter")                = n_iter,
    Rcpp::Named("lambda")                = lambda_reg,
    Rcpp::Named("objective_path")        = objective_path,
    Rcpp::Named("features_path")         = features_path,
    Rcpp::Named("gradient_mapping_path") = gradient_mapping_path,
    Rcpp::Named("final_gradient_mapping") =
        gradient_mapping_path.empty() ? NA_REAL : gradient_mapping_path.back(),
    Rcpp::Named("status")                = status
  );

  return Rcpp::List::create(
    Rcpp::Named("coef")      = x_coef,
    Rcpp::Named("intercept") = R_NilValue,   // Cox: no intercept
    Rcpp::Named("info")      = info
  );
}
