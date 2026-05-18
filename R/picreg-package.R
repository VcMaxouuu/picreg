#' picreg: Variable Selection using the Pivotal Information Criterion
#'
#' Sparse regression and classification with automatic regularisation
#' parameter selection via the Pivotal Detection Boundary (PDB) method.
#' Implements FISTA optimisation for L1, SCAD and MCP penalties across
#' Gaussian, Binomial, Poisson, Exponential, Gumbel and Cox families.
#'
#' @keywords internal
#' @useDynLib picreg, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rnorm rbinom rpois rexp runif sd qnorm quantile aggregate stepfun
#' @importFrom graphics abline axis box hist legend lines par points segments
#' @importFrom grDevices adjustcolor grey.colors hcl.colors
"_PACKAGE"
