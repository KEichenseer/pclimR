#' The Main MCMC function for the ALPHA BETA model
#'
#' \code{gradient} runs the MCMC algorithm to estimate ALPHA (shift)
#' and BETA (stretch) of a correlation section to match the reference section.
#'
#' @param x dataset as a \code{CSMproxydata} object.
#' @param coeff \code{CSMautoinits} object: if supplied many values can
#'                    be inherited from this (marked ** below). Not currently
#'                    being used.

#' @return An object of class CSM_MCMC consisting of a list containing:
#'   \item{call}{the parameters of the call.}
#'   \item{ALPHA}{2D numeric array of ALPHA at each iteration. Contains one row per chain.

#' @import fda msm conquer tidyquant zoo splines MCMCpack mvnfast
#'
#' @export
#'
#'


gradient <- function(x, coeff, sdy) { # sigma is labelled "sdy"
  A = coeff[1]
  K = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  nu = coeff[5]
  return(A + max(c(K-A,0))/((1+(nu*exp(Q*(x-M))))^(1/nu)) + rnorm(length(x),0,sdy))
}
