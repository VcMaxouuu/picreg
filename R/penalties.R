#' Sparsity-inducing penalties for pic.
#'
#' Three penalties are supported and identified by lowercase name to match
#' the C++ registry: `"lasso"`, `"scad"`, `"mcp"`. Each penalty object is a
#' lightweight list carrying its name and concavity parameters. The actual
#' evaluation and proximal operators live in C++ (`src/penalty_*.cpp`).
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
