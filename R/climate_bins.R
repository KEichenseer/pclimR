#' Climatic bins
#'
#' A function which return climatic bins for a given time interval.
#'
#' @param interval \code{character or numeric}. Interval name or age available in \code{GTS2020}. If a single interval name is provided, this interval is returned.
#' If two interval names are provided, these intervals and those existing between are returned. If a single interval age is provided, the age matching this interval is returned.
#' If two interval ages are provided, the intervals occurring in the range of these ages are returned.
#' @return A \code{data.frame} of climatic bins for a specified interval.
#' @examples
#' Using a single interval
#' climate_bins(interval = "Maastrichtian")
#' Using multiple intervals
#' climate_bins(interval = c("Albian", "Maastrichtian"))
#' Using interval midpoint age
#' climate_bins(interval = 100)
#' Using interval age range
#' climate_bins(interval = c(50, 100))
climate_bins <- function(interval){
  #download data
  df <- read.csv("https://drive.google.com/uc?export=download&id=1_bhI10hJRp_21dGLLyoX-C_jJQHUyjcz")
  #subset to stages
  df <- subset(df, RANK == "Stage")
  #character string entered
  if(is.character(interval)){
    if(length(interval) == 1){
    w <- which(df$NAME %in% interval)
      if(length(w) != length(interval)){stop(paste0("Checking spelling of specified intervals. \n Available intervals are: ",
                                                    toString(df$NAME)))}
    df <- df[w,c("STAGE_NUMBER", "NAME", "BASE_Ma", "MIDPOINT_Ma", "DURATION_Myr",
                 "CLIMATE_BIN", "CLIMATE_STATE", "FONT_COLOUR", "HTML_RGB_HASH")]
    }
    if(length(interval) == 2){
      w <- which(df$NAME %in% interval)
      df <- df[w[1]:w[2], c("STAGE_NUMBER", "NAME", "BASE_Ma", "MIDPOINT_Ma", "DURATION_Myr",
                              "CLIMATE_BIN", "CLIMATE_STATE", "FONT_COLOUR", "HTML_RGB_HASH")]
      }
  }
  #numeric ages entered
  if(is.numeric(interval)){
    if(length(interval) == 1){
      if(interval > max(df$BASE_Ma) | interval < min(df$TOP_Ma)){
        stop("Value does not appear in range of available intervals: 0 to 541")
      }
     df <- df[which.min(abs(df$MIDPOINT_Ma - interval)),
              c("STAGE_NUMBER", "NAME", "BASE_Ma", "MIDPOINT_Ma", "DURATION_Myr",
                "CLIMATE_BIN", "CLIMATE_STATE", "FONT_COLOUR", "HTML_RGB_HASH"
              )]
    }

    if(length(interval) == 2){
      max_int <- max(interval)
      min_int <- min(interval)

      if(max_int > max(df$BASE_Ma) | min_int < min(df$TOP_Ma)){
        stop("Values do not appear in range of available intervals: 0 to 541")
      }

      df <- df[which.min(abs(df$TOP_Ma - min_int)):which.min(abs(df$BASE_Ma - max_int)),
         c("STAGE_NUMBER", "NAME", "BASE_Ma", "MIDPOINT_Ma", "DURATION_Myr",
           "CLIMATE_BIN", "CLIMATE_STATE", "FONT_COLOUR", "HTML_RGB_HASH"
         )]
    }

    if(length(interval) > 2){
      stop("interval must be a character vector or a numeric vector of length 1 or 2")
    }

  }
    return(df)
  }
