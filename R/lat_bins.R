#' Latitudinal bins
#'
#' A function to generate latitudinal bins for a given bin size.
#'
#' @param size A single numeric value of > 0 and <= 90.
#' @param fit A logical value indicating whether a bin size check should be performed to ensure that the
#' entire latitudinal range is covered (-90 to 90).
#' @return A \code{data.frame} of latitudinal bins for a given size.
#' @examples
#' lat_bins(size = 20)
#' lat_bins(size = 13, fit = TRUE)

lat_bins <- function(size = 10, fit = FALSE){
  bins <- 180/size
  if(fit == TRUE){
    if(is.integer(bins) == FALSE){
      int <- 180/seq(from = 1, to = 90, by = 1)
      int <- which(int%%1==0)
      size <- int[which.min(abs(int - size))]
      bins <- 180/size
      message(paste0("Bin size set to ",size," degrees to cover latitudinal range."))
    }
  }
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
