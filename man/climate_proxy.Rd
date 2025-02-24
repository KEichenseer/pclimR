% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/climate_proxy.R
\name{climate_proxy}
\alias{climate_proxy}
\title{Climate proxy data}
\usage{
climate_proxy(x, available = FALSE)
}
\arguments{
\item{x}{\code{logical}. Should a dataframe detailing available datasets be returned? Defaults to \code{FALSE}.}
}
\value{
A \code{data.frame} of climate proxy data for a specified proxy.
}
\description{
This function downloads a specified dataset from those available. Available datasets can also be called via this function.
}
\examples{
climate_proxy(available = TRUE)
climate_proxy(x = "D18O")
}
