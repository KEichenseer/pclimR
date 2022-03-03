#' Latitudinal bins
#'
#' A function to generate latitudinal bins for a given bin size.
#'
#' @param size A single numeric value of > 0 and <= 90.
#' @return A \code{data.frame} of latitudinal bins for a given size.
#' @examples
#' lat_bins(size = 10)
#' lat_bins(size = 30)

lat_bins <- function(size = 10){
  bins <- 180/size
  df <- seq(from = -90, to = 90, by = size)
  min <- df[1:bins]
  max <- df[1:bins]+size
  mid <- (max+min)/2
  bin <- 1:bins
  df <- cbind(max, mid, min)
  df <- df[order(-max),]
  df <- cbind.data.frame(bin, df)
  return(df)
}
