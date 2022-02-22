#' Title
#'
#' @param x
#' @param coeff
#' @param sdy
#'
#' @return
#' @export
#'
#' @examples
#'


gradient <- function(x, coeff, sdy) { # sigma is labelled "sdy"
  A = coeff[1]
  K = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  nu = coeff[5]
  return(A + max(c(K-A,0))/((1+(nu*exp(Q*(x-M))))^(1/nu)) + rnorm(length(x),0,sdy))
}

