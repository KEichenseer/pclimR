#' gradient generates a latitudinal temperature gradient
#' using a generalised logistic function
#'
#' @param x absolute latitudes
#' @param coeff coefficients of the generalised logistic function
#' @param sdy standard deviation of the points from the curve
#'
#' @return one temperature values for each latitude x, based on coeff and sdy
#' @export
#'
#' @examples
#'


gradient <- function(x, coeff, sdy) { # sigma is labelled "sdy"
  A = coeff[1]
  K = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  return(A + max(c(K-A,0))/((1+(exp(Q*(x-M))))) + rnorm(length(x),0,sdy))
}
