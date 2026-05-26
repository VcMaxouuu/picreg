#' Small Gaussian dataset for the introductory vignette.
#'
#' A synthetic Gaussian regression dataset used to illustrate `pic()`
#' throughout the introductory vignette. It contains \eqn{n = 100}
#' observations of \eqn{p = 30} predictors, of which only \eqn{s = 5}
#' carry signal; the remaining 25 are pure noise.
#'
#' Column names are chosen to make the underlying support obvious at a
#' glance:
#'
#' \itemize{
#'   \item \code{gene_1, ..., gene_5}: the five active variables, whose
#'         non-zero coefficients are drawn uniformly in
#'         \eqn{[0.5,\, 1.5]} with random sign.
#'   \item \code{noise_1, ..., noise_25}: the remaining inactive
#'         variables, with true coefficient \eqn{0}.
#' }
#'
#' The columns are interleaved in random order; column names are the
#' only indicator of which features are part of the true support.
#'
#' The response is generated as \eqn{y = X\beta + \varepsilon} with
#' \eqn{\varepsilon \sim \mathcal{N}(0, 1)}.
#'
#' @format A list with two components:
#' \describe{
#'   \item{\code{X}}{Numeric matrix of dimension \eqn{100 \times 30}
#'         with column names \code{gene_*} and \code{noise_*}.}
#'   \item{\code{y}}{Numeric vector of length \eqn{100}.}
#' }
#'
#' @examples
#' data(QuickStartExample)
#' fit <- pic(QuickStartExample$X, QuickStartExample$y)
#' fit$selected
"QuickStartExample"

#' Small survival dataset for the Cox section of the vignette.
#'
#' A synthetic survival dataset used to illustrate `pic()` with the Cox
#' family. It contains \eqn{n = 250} subjects observed on \eqn{p = 50}
#' covariates, of which \eqn{s = 5} carry signal; the remaining 45 are
#' noise.
#'
#' Column names follow the same convention as
#' \code{\link{QuickStartExample}}: active variables are labelled
#' \code{gene_1, ..., gene_5} and inactive ones \code{noise_1, ...,
#' noise_45}, interleaved in random order.
#'
#' Event times are drawn from an exponential proportional-hazards
#' model
#' \deqn{T_i \sim \mathrm{Exp}\!\bigl(e^{X_i\beta}\bigr),}
#' and independent censoring times from \eqn{C_i \sim
#' \mathrm{Exp}(0.5)}. The observed response is the standard two-column
#' \eqn{(\min(T_i, C_i),\, \mathbf{1}\{T_i \le C_i\})}. The censoring
#' rate is roughly \eqn{40\%}.
#'
#' @format A list with two components:
#' \describe{
#'   \item{\code{X}}{Numeric matrix of dimension \eqn{250 \times 50}
#'         with column names \code{gene_*} and \code{noise_*}.}
#'   \item{\code{y}}{Numeric matrix of dimension \eqn{250 \times 2}
#'         with columns \code{time} and \code{event}.}
#' }
#'
#' @examples
#' data(CoxExample)
#' fit <- pic(CoxExample$X, CoxExample$y, family = "cox")
#' fit$selected
"CoxExample"

#' Small binary-classification dataset for the Binomial section of the vignette.
#'
#' A synthetic logistic-regression dataset used to illustrate `pic()` with
#' the Binomial family. It contains \eqn{n = 300} observations of
#' \eqn{p = 50} predictors, of which \eqn{s = 5} carry signal; the
#' remaining 45 are noise.
#'
#' Column names follow the same convention as
#' \code{\link{QuickStartExample}}: active variables are labelled
#' \code{gene_1, ..., gene_5} and inactive ones \code{noise_1, ...,
#' noise_45}, interleaved in random order.
#'
#' The binary response is generated as
#' \deqn{Y_i \sim \mathrm{Bernoulli}\!\left(\frac{1}{1 + e^{-X_i\beta}}\right),}
#' with non-zero coefficients drawn uniformly in \eqn{[1.5,\, 3]} with
#' random sign. The class balance is roughly \eqn{45\%} of positives.
#'
#' @format A list with two components:
#' \describe{
#'   \item{\code{X}}{Numeric matrix of dimension \eqn{300 \times 50}
#'         with column names \code{gene_*} and \code{noise_*}.}
#'   \item{\code{y}}{Integer vector of length \eqn{300} containing
#'         \eqn{0/1} class labels.}
#' }
#'
#' @examples
#' data(BinomialExample)
#' fit <- pic(BinomialExample$X, BinomialExample$y, family = "binomial")
#' fit$selected
"BinomialExample"
