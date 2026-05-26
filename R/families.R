#' GLM families for pic — descriptor layer.
#'
#' Each family is a thin descriptor carrying:
#'   * `name` - used by the dispatchers to route to the C++ implementation.
#'   * `g`    - mean-function link (object with `name` and callable `fn`).
#'              `g$fn(eta)` is applied at predict time to obtain the mean
#'              response.
#'   * `phi`  - variance-stabilising transform (object with `name`).
#'              Purely informational on the R side; the actual stabilisation
#'              happens inside the C++ math.
#'
#' Typing `fit$family` gives a one-glance summary of the model's family
#' (see [print.pic.family()]). The actual loss / gradient math lives in C++:
#'   * `src/family_*.cpp`     — Gaussian, Binomial, Poisson, Exponential, Gumbel.
#'   * `src/cox.cpp`          — Cox.
#'
#' @name pic_families
NULL

# ---- internal link descriptors ------------------------------------------
# Each carries a display name (`name`) and a callable (`fn`). Only `g$fn` is
# actually called from R (in `predict.pic`); `phi$fn` is informational.

.identity_link <- list(name = "identity", fn = function(eta) eta)
.logistic_link <- list(name = "logistic", fn = function(eta) 1 / (1 + exp(-eta)))
.exp_link      <- list(name = "exp",      fn = function(eta) exp(eta))
.sqrt_link     <- list(name = "sqrt",     fn = function(x)   sqrt(x))

# ---- family descriptors ------------------------------------------------

.make_family <- function(name, g, phi) {
  structure(
    list(name = name, g = g, phi = phi),
    class = c(paste0("pic.family.", name), "pic.family")
  )
}

gaussian    <- function() .make_family("gaussian",    g = .identity_link, phi = .sqrt_link)
binomial    <- function() .make_family("binomial",    g = .logistic_link, phi = .identity_link)
poisson     <- function() .make_family("poisson",     g = .exp_link,      phi = .identity_link)
exponential <- function() .make_family("exponential", g = .exp_link,      phi = .identity_link)
gumbel      <- function() .make_family("gumbel",      g = .identity_link, phi = .exp_link)
cox         <- function() .make_family("cox",         g = .identity_link, phi = .sqrt_link)

# ---- registry / dispatcher ---------------------------------------------

#' Resolve a family name into a pic family descriptor.
#'
#' Accepts one of the six supported names: `"gaussian"`, `"binomial"`,
#' `"poisson"`, `"exponential"`, `"gumbel"`, `"cox"`. Anything else raises
#' an error — the build does not support user-defined families.
#' Already-built descriptors are returned unchanged (internal re-entry).
#'
#' @param family A family name or already-built pic family descriptor.
#' @return A pic family descriptor.
#' @keywords internal
get_family <- function(family) {
  if (inherits(family, "pic.family")) return(family)
  if (!is.character(family) || length(family) != 1L)
    stop("`family` must be one of \"gaussian\", \"binomial\", \"poisson\", ",
         "\"exponential\", \"gumbel\", \"cox\".", call. = FALSE)
  family <- match.arg(
    family,
    c("gaussian", "binomial", "poisson", "exponential", "gumbel", "cox")
  )
  switch(
    family,
    gaussian    = gaussian(),
    binomial    = binomial(),
    poisson     = poisson(),
    exponential = exponential(),
    gumbel      = gumbel(),
    cox         = cox()
  )
}

#' Pretty-print a pic family descriptor.
#'
#' Three-line summary showing the family name and its (g, phi) link pair.
#'
#' @param x A `pic.family` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.pic.family <- function(x, ...) {
  cat("family: ", x$name, " (link g = ", x$g$name, ", phi = ", x$phi$name, ")",
      sep = "")
  invisible(x)
}
