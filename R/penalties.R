#' Sparsity-inducing penalties for pic.
#'
#' Three penalties are supported, identified by lowercase name to match
#' the C++ registry. Each penalty enters the [pic()] objective as
#' \deqn{\mathrm{pen}(\beta) = \sum_{j=1}^p p_\lambda(|\beta_j|),}
#' where \eqn{p_\lambda(\cdot)} depends on the penalty.
#'
#' \describe{
#'   \item{`"lasso"`}{
#'     L1 (soft-thresholding) penalty:
#'     \deqn{p_\lambda(|t|) = \lambda |t|.}
#'     Convex, gives the strongest shrinkage on large coefficients 
#'     bias does not vanish as \eqn{|t| \to \infty}.
#'   }
#'   \item{`"scad"` (Smoothly Clipped Absolute Deviation, Fan & Li 2001)}{
#'     Non-convex penalty with concavity parameter `scad_a > 2`
#'     (default 3.7):
#'     \deqn{p_\lambda'(|t|) = \lambda\!\left\{
#'        \mathbf{1}\{|t| \le \lambda\}
#'        + \frac{(a\lambda - |t|)_+}{(a - 1)\lambda}
#'          \mathbf{1}\{|t| > \lambda\}\right\}.}
#'     Behaves like the lasso for small \eqn{|t|}, then tapers off so
#'     large coefficients are barely penalized - yields nearly unbiased
#'     estimates on strong signals.
#'   }
#'   \item{`"mcp"` (Minimax Concave Penalty, Zhang 2010)}{
#'     Non-convex penalty with concavity parameter `mcp_gamma > 1`
#'     (default 3.0):
#'     \deqn{p_\lambda'(|t|) = \left(\lambda - \frac{|t|}{\gamma}\right)_+.}
#'     Similar motivation as SCAD but a smoother transition: starts at
#'     the lasso derivative for small \eqn{|t|} and tapers linearly to
#'     zero at \eqn{|t| = \gamma\lambda}.
#'   }
#' }
#'
#' The actual evaluation and proximal operators live in C++
#' (`src/penalty_*.cpp`). Larger `scad_a` / `mcp_gamma` make the
#' penalty closer to the lasso; smaller values amplify the
#' non-convexity (and the bias reduction on strong signals).
#'
#' @name pic_penalties
NULL


# Internal constructors — return a thin descriptor (no R-side prox/evaluate).
.lasso <- function() {
  structure(list(name = "lasso", params = list()),
            class = c("pic.penalty.lasso", "pic.penalty"))
}

.scad <- function(a = 3.7) {
  if (a <= 2) stop("SCAD requires a > 2; got a = ", a, call. = FALSE)
  structure(list(name = "scad", params = list(a = a)),
            class = c("pic.penalty.scad", "pic.penalty"))
}

.mcp <- function(gamma = 3.0) {
  if (gamma <= 1) stop("MCP requires gamma > 1; got gamma = ", gamma,
                       call. = FALSE)
  structure(list(name = "mcp", params = list(gamma = gamma)),
            class = c("pic.penalty.mcp", "pic.penalty"))
}

#' @export
print.pic.penalty <- function(x, ...) {
  if (length(x$params) == 0L) {
    cat(x$name, "()\n", sep = "")
  } else {
    args <- paste(names(x$params), unlist(x$params),
                  sep = " = ", collapse = ", ")
    cat(x$name, "(", args, ")\n", sep = "")
  }
  invisible(x)
}

#' Resolve a penalty name into a pic penalty descriptor.
#'
#' Accepts one of `"lasso"`, `"scad"`, `"mcp"` (case-insensitive). Anything
#' else raises an error: the build does not support user-defined penalties.
#'
#' @param penalty A penalty name.
#' @param scad_a SCAD concavity parameter (default 3.7).
#' @param mcp_gamma MCP concavity parameter (default 3.0).
#' @return A `pic.penalty` descriptor.
#' @keywords internal
get_penalty <- function(penalty, scad_a = 3.7, mcp_gamma = 3.0) {
  if (!is.character(penalty) || length(penalty) != 1L)
    stop("`penalty` must be one of \"lasso\", \"scad\", \"mcp\".",
         call. = FALSE)
  key <- match.arg(tolower(penalty), c("lasso", "scad", "mcp"))
  switch(
    key,
    lasso = .lasso(),
    scad  = .scad(a = scad_a),
    mcp   = .mcp(gamma = mcp_gamma)
  )
}
