#' Non-monotone Accelerated Proximal Gradient (nmAPG) solver — dispatcher.
#'
#' Thin R wrapper around the C++ `fista_cpp` implementation. The actual
#' algorithm (Li & Lin, NeurIPS 2015) lives in `src/fista.cpp`; this
#' function only extracts the C++-facing arguments and routes the call.
#'
#' Cox is routed to a dedicated C++ FISTA (`fista_cox_cpp`) because its
#' loss has no intercept and uses a 2-column `(time, event)` response.
#' The build accepts only the six families and three penalties that have
#' a C++ implementation; the family / penalty registries throw on any
#' other input upstream.
#'
#' @param X Standardized design matrix.
#' @param y Response vector (length n).
#' @param family A pic family object (provides `name`).
#' @param penalty A pic penalty object (provides `name`, `params`).
#' @param lambda_reg Regularization parameter.
#' @param fit_intercept Whether to update an unpenalized intercept.
#' @param rel_tol Relative gradient-mapping tolerance.
#' @param step_size_init Initial step size.
#' @param max_iter Hard cap on outer iterations.
#' @param eta_param Memory weight of the non-monotone reference (in [0,1)).
#' @param delta_param Sufficient-descent constant in the line search.
#' @param rho Step-size reduction factor (< 1).
#' @param bb_growth_cap Maximum BB step growth over the previous step.
#' @param coef_init,intercept_init Optional warm starts.
#' @return A list with `coef`, `intercept`, `info`.
#' @keywords internal
fista <- function(X, y, family, penalty, lambda_reg,
                  fit_intercept = TRUE,
                  rel_tol = 1e-4,
                  step_size_init = 1e-1,
                  max_iter = 500L,
                  eta_param = 0.8,
                  delta_param = 1e-4,
                  rho = 0.5,
                  bb_growth_cap = 4.0,
                  coef_init = NULL,
                  intercept_init = NULL) {
  scad_a    <- if (!is.null(penalty$params$a))     penalty$params$a     else 3.7
  mcp_gamma <- if (!is.null(penalty$params$gamma)) penalty$params$gamma else 3.0

  # Cox has its own dedicated C++ FISTA (separate y interface, no intercept).
  if (family$name == "cox") {
    if (!is.matrix(y) || ncol(y) != 2L)
      stop("Cox: y must be a 2-column matrix (time, event).", call. = FALSE)
    return(fista_cox_cpp(
      X              = X,
      times          = y[, 1L],
      events         = y[, 2L],
      penalty_name   = penalty$name,
      lambda_reg     = lambda_reg,
      scad_a         = scad_a,
      mcp_gamma      = mcp_gamma,
      rel_tol        = rel_tol,
      step_size_init = step_size_init,
      max_iter       = as.integer(max_iter),
      eta_param      = eta_param,
      delta_param    = delta_param,
      rho            = rho,
      bb_growth_cap  = bb_growth_cap,
      coef_init      = coef_init
    ))
  }

  y_vec <- if (is.matrix(y)) y[, 1L] else as.numeric(y)

  fista_cpp(
    X              = X,
    y              = y_vec,
    family_name    = family$name,
    penalty_name   = penalty$name,
    lambda_reg     = lambda_reg,
    scad_a         = scad_a,
    mcp_gamma      = mcp_gamma,
    fit_intercept  = fit_intercept,
    rel_tol        = rel_tol,
    step_size_init = step_size_init,
    max_iter       = as.integer(max_iter),
    eta_param      = eta_param,
    delta_param    = delta_param,
    rho            = rho,
    bb_growth_cap  = bb_growth_cap,
    coef_init      = coef_init,
    intercept_init = intercept_init
  )
}
