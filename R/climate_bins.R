#' Climatic bins
#'
#' A function which return climatic bins for a given time interval.
#'
#' @param interval A character string matching an interval available in NAME of \code{GTS2020}.
#' @return A \code{data.frame} of climatic bins for a specified interval.
#' @examples
#' lat_bins(size = 20)
#' lat_bins(size = 13, fit = TRUE)
climate_bins <- function(interval){
  df <- GTS2020
  w <- which(df$NAME == interval)
  df <- subset(df, MIDPOINT_Ma >= df$TOP_Ma[w] & MIDPOINT_Ma <= df$BASE_Ma)
  df <- subset(df, RANK == "Stage")
  df <- df[,c("INDEX")]
}
