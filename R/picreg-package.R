#' @keywords internal
#' @useDynLib picreg, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats predict setNames aggregate quantile rnorm rbinom rpois rexp runif median sd qnorm pnorm dnorm var coef
#' @importFrom graphics par axis box abline segments points hist plot.new plot.window mtext text lines polygon legend strwidth strheight
#' @importFrom grDevices adjustcolor
#' @importFrom utils head tail
"_PACKAGE"
