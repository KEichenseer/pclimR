#' Climate proxy data
#'
#' This function downloads a specified dataset from those available. Available datasets can also be called via this function.
#'
#' @param x \code{character}. A specified ID from those available. Available datasets can be viewed via the \code{available} parameter.
#' @param x \code{logical}. Should a dataframe detailing available datasets be returned? Defaults to \code{FALSE}.
#' @return A \code{data.frame} of climate proxy data for a specified proxy.
#' @examples
#' climate_proxy(available = TRUE)
#' climate_proxy(x = "D18O")
climate_proxy <- function(x, available = FALSE){
  #available datasets
  datasets <- c("D18O")
  #check available datasets
  if(available == TRUE){
    #datasets <- read.csv() available datasets to be added here
    return(datasets)
  }
  #check if specified dataset is available
  if(!(x %in% datasets)){
    stop(paste0("x must equal: ", datasets[1]))
  }
  #retrieve specified dataset
  if(x == "D18O"){
    df <- read.csv(file = "https://drive.google.com/uc?export=download&id=1WZu0PUitdpInYxxfK__lcsW8KA0jn5vS")
  }
  return(df)
}

