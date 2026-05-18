// Non-monotone Accelerated Proximal Gradient (nmAPG) solver — C++.
//
// Faithful port of R/fista.R (Li & Lin, NeurIPS 2015 Algorithm 4).
// Consumes the picr Family and Penalty systems for everything beyond
// pure FISTA mechanics.

// [[Rcpp::depends(RcppArmadillo)]]
#include "family_base.h"
#include "penalty_base.h"

using namespace arma;

// ---------- helpers --------------------------------------------------------

// Joint Barzilai–Borwein step over (coef, intercept). Clipped to
// [alpha_min, min(alpha_max, growth_cap * fallback, alpha_pen_cap)] so the
// step never grows uncontrollably nor leaves the penalty's safe domain.
static double bb_step(const vec& coef_a, const vec& coef_b,
                      double int_a, double int_b,
                      const vec& grad_a_coef, const vec& grad_b_coef,
                      double grad_a_int, double grad_b_int,
                      double fallback, bool fit_intercept,
                      double growth_cap, double alpha_pen_cap,
                      double alpha_min = 1e-10, double alpha_max = 1e10,
                      double eps = 1e-12) {
  vec s_coef = coef_a - coef_b;
  vec r_coef = grad_a_coef - grad_b_coef;
  double ss = dot(s_coef, s_coef);
  double sr = dot(s_coef, r_coef);
  if (fit_intercept) {
    double s_int = int_a - int_b;
    double r_int = grad_a_int - grad_b_int;
    ss += s_int * s_int;
    sr += s_int * r_int;
  }
  if (ss <= eps || sr <= eps) {
    return std::min(fallback, alpha_pen_cap);
  }
  double upper = std::min({alpha_max, growth_cap * fallback, alpha_pen_cap});
  return std::max(alpha_min, std::min(ss / sr, upper));
}

// Result of the line search from the extrapolated point y_k.
struct LsResult {
  vec    z_coef;
  double z_int;
  double F_z;
  double alpha;
  bool   non_monotone;
};

// Backtracking proximal line search from y_k, with descent-lemma check
// and a non-monotone test against the running reference c_k.
static LsResult ls_from_y(const mat& X, const vec& y,
                          const picr::Family& fam, const picr::Penalty& pen,
                          double lambda_reg,
                          const vec& y_coef, double y_int,
                          double f_y,
                          const vec& grad_y_coef, double grad_y_int,
                          double alpha, double c_k,
                          double delta_param, double rho,
                          bool fit_intercept, int max_inner = 30) {
  vec    z_coef = y_coef;
  double z_int  = y_int;
  double f_z    = f_y;
  
  for (int i = 0; i < max_inner; ++i) {
    z_coef = pen.prox(y_coef - alpha * grad_y_coef, lambda_reg, alpha,
                      pen.param);
    z_int  = fit_intercept ? (y_int - alpha * grad_y_int) : 0.0;
    
    vec lp = X * z_coef + z_int;
    f_z = fam.evaluate(y, lp);
    
    vec d_coef = z_coef - y_coef;
    double inner_p = dot(d_coef, grad_y_coef);
    double diff_sq = dot(d_coef, d_coef);
    if (fit_intercept) {
      double d_int = z_int - y_int;
      inner_p += d_int * grad_y_int;
      diff_sq += d_int * d_int;
    }
    
    if (f_z <= f_y + inner_p + diff_sq / (2.0 * alpha)) {
      double F_z_val = f_z + pen.evaluate(z_coef, lambda_reg, pen.param);
      bool nm = (F_z_val <= c_k - delta_param * diff_sq);
      return { z_coef, z_int, F_z_val, alpha, nm };
    }
    alpha *= rho;
  }
  
  // Max inner reached: accept last trial, treat as not non-monotone so the
  // v-branch fallback runs.
  double F_z_val = f_z + pen.evaluate(z_coef, lambda_reg, pen.param);
  return { z_coef, z_int, F_z_val, alpha, false };
}

// Result of the v-branch fallback line search from x_k.
struct VResult {
  vec    coef;
  double intercept;
  double F_val;
  double alpha;
  char   picked;        // 'z' or 'v'
};

static VResult v_branch(const mat& X, const vec& y,
                        const picr::Family& fam, const picr::Penalty& pen,
                        double lambda_reg,
                        const vec& x_coef, double x_int,
                        const vec& y_prev_coef, double y_prev_int,
                        const vec& grad_y_prev_coef, double grad_y_prev_int,
                        const vec& z_next_coef, double z_next_int, double F_z,
                        double alpha_x, double rho,
                        bool fit_intercept,
                        double bb_growth_cap, double alpha_pen_cap,
                        int max_inner = 30) {
  vec lp_x = X * x_coef + x_int;
  double f_x = fam.evaluate(y, lp_x);
  picr::GradResult gx = fam.grad(X, y, lp_x);
  
  alpha_x = bb_step(x_coef, y_prev_coef, x_int, y_prev_int,
                    gx.coef, grad_y_prev_coef, gx.intercept, grad_y_prev_int,
                    alpha_x, fit_intercept, bb_growth_cap, alpha_pen_cap);
  alpha_x = std::min(alpha_x, alpha_pen_cap);
  
  vec    v_coef = x_coef;
  double v_int  = x_int;
  double f_v    = f_x;
  double F_v    = std::numeric_limits<double>::infinity();
  
  for (int i = 0; i < max_inner; ++i) {
    v_coef = pen.prox(x_coef - alpha_x * gx.coef, lambda_reg, alpha_x,
                      pen.param);
    v_int  = fit_intercept ? (x_int - alpha_x * gx.intercept) : 0.0;
    
    vec lp = X * v_coef + v_int;
    f_v = fam.evaluate(y, lp);
    
    vec d_coef = v_coef - x_coef;
    double inner_p = dot(d_coef, gx.coef);
    double diff_sq = dot(d_coef, d_coef);
    if (fit_intercept) {
      double d_int = v_int - x_int;
      inner_p += d_int * gx.intercept;
      diff_sq += d_int * d_int;
    }
    
    if (f_v <= f_x + inner_p + diff_sq / (2.0 * alpha_x)) {
      F_v = f_v + pen.evaluate(v_coef, lambda_reg, pen.param);
      break;
    }
    alpha_x *= rho;
  }
  
  // Whether ls converged or hit max_inner, F_v holds the best trial we have.
  if (!std::isfinite(F_v)) {
    F_v = f_v + pen.evaluate(v_coef, lambda_reg, pen.param);
  }
  
  if (F_z <= F_v) {
    return { z_next_coef, z_next_int, F_z, alpha_x, 'z' };
  }
  return { v_coef, v_int, F_v, alpha_x, 'v' };
}

// ---------- main entry point ----------------------------------------------

// [[Rcpp::export]]
Rcpp::List fista_cpp(const arma::mat& X,
                     const arma::vec& y,
                     std::string family_name,
                     std::string penalty_name,
                     double lambda_reg,
                     double scad_a       = 3.7,
                     double mcp_gamma    = 3.0,
                     bool   fit_intercept = true,
                     double rel_tol      = 1e-4,
                     double step_size_init = 1e-2,
                     int    max_iter     = 500,
                     double eta_param    = 0.8,
                     double delta_param  = 1e-4,
                     double rho          = 0.5,
                     double bb_growth_cap = 2.0,
                     Rcpp::Nullable<Rcpp::NumericVector> coef_init     = R_NilValue,
                     Rcpp::Nullable<double>              intercept_init = R_NilValue) {
  
  picr::Family  fam = picr::get_family(family_name);
  picr::Penalty pen = picr::get_penalty(penalty_name, scad_a, mcp_gamma);
  
  const double step_margin = 0.95;
  const double alpha_pen_cap = step_margin * pen.max_step;
  step_size_init = std::min(step_size_init, alpha_pen_cap);
  
  const int n_features = X.n_cols;
  
  // ---- init coef / intercept ----
  vec x_coef = (coef_init.isNull())
    ? vec(n_features, fill::zeros)
      : Rcpp::as<vec>(coef_init);
  
  double x_int = 0.0;
  if (fit_intercept) {
    x_int = intercept_init.isNull()
    ? fam.starting_intercept(y)
      : Rcpp::as<double>(intercept_init);
  }
  
  vec    x_prev_coef = x_coef;
  double x_prev_int  = x_int;
  vec    z_coef      = x_coef;
  double z_int       = x_int;
  
  vec    y_prev_coef = x_coef;
  double y_prev_int  = x_int;
  
  vec lp = X * x_coef + x_int;
  picr::GradResult g0 = fam.grad(X, y, lp);
  vec    grad_y_prev_coef = g0.coef;
  double grad_y_prev_int  = g0.intercept;
  
  double t_prev = 0.0, t_curr = 1.0;
  double F_x = fam.evaluate(y, lp) + pen.evaluate(x_coef, lambda_reg, pen.param);
  double c_k = F_x, q_k = 1.0;
  
  double alpha_y = step_size_init;
  double alpha_x = step_size_init;
  
  for (int k = 1; k <= max_iter; ++k) {
    double a = t_prev / t_curr;
    double b = (t_prev - 1.0) / t_curr;
    
    vec    y_coef = x_coef + a * (z_coef - x_coef) + b * (x_coef - x_prev_coef);
    double y_int  = fit_intercept
    ? (x_int + a * (z_int - x_int) + b * (x_int - x_prev_int))
      : 0.0;
    
    vec lp_y = X * y_coef + y_int;
    double f_y = fam.evaluate(y, lp_y);
    picr::GradResult gy = fam.grad(X, y, lp_y);
    
    alpha_y = bb_step(y_coef, y_prev_coef, y_int, y_prev_int,
                      gy.coef, grad_y_prev_coef, gy.intercept, grad_y_prev_int,
                      alpha_y, fit_intercept, bb_growth_cap, alpha_pen_cap);
    
    LsResult ls_y = ls_from_y(X, y, fam, pen, lambda_reg,
                              y_coef, y_int, f_y,
                              gy.coef, gy.intercept,
                              std::min(alpha_y, alpha_pen_cap),
                              c_k, delta_param, rho, fit_intercept);
    alpha_y = ls_y.alpha;
    
    vec    x_next_coef;
    double x_next_int;
    double F_next;
    vec    gm_anchor_coef;
    double gm_alpha;
    
    if (ls_y.non_monotone) {
      x_next_coef    = ls_y.z_coef;
      x_next_int     = ls_y.z_int;
      F_next         = ls_y.F_z;
      gm_anchor_coef = y_coef;
      gm_alpha       = alpha_y;
    } else {
      VResult vb = v_branch(X, y, fam, pen, lambda_reg,
                            x_coef, x_int,
                            y_prev_coef, y_prev_int,
                            grad_y_prev_coef, grad_y_prev_int,
                            ls_y.z_coef, ls_y.z_int, ls_y.F_z,
                            alpha_x, rho, fit_intercept,
                            bb_growth_cap, alpha_pen_cap);
      x_next_coef = vb.coef;
      x_next_int  = vb.intercept;
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
    
    // ---- gradient-mapping convergence ----
    vec    G_coef    = (gm_anchor_coef - x_next_coef) / gm_alpha;
    double G_norm    = std::sqrt(dot(G_coef, G_coef));
    double beta_norm = std::sqrt(dot(gm_anchor_coef, gm_anchor_coef));
    double rel_gm    = G_norm / std::max(beta_norm, 1.0);
    
    if (rel_gm < rel_tol) {
      x_coef = x_next_coef;
      x_int  = x_next_int;
      break;
    }
    
    // ---- state roll-forward ----
    x_prev_coef = x_coef;
    x_prev_int  = x_int;
    z_coef      = ls_y.z_coef;
    z_int       = ls_y.z_int;
    y_prev_coef = y_coef;
    y_prev_int  = y_int;
    grad_y_prev_coef = gy.coef;
    grad_y_prev_int  = gy.intercept;
    x_coef      = x_next_coef;
    x_int       = x_next_int;
    
    double t_next = (std::sqrt(4 * t_curr * t_curr + 1.0) + 1.0) / 2.0;
    t_prev = t_curr;
    t_curr = t_next;
    
    double q_next = eta_param * q_k + 1.0;
    c_k = (eta_param * q_k * c_k + F_next) / q_next;
    q_k = q_next;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("coef")      = x_coef,
    Rcpp::Named("intercept") = fit_intercept ? Rcpp::wrap(x_int)
      : R_NilValue
  );
}