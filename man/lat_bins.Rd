% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lat_bins.R
\name{lat_bins}
\alias{lat_bins}
\title{Latitudinal bins}
\usage{
lat_bins(size = 10, fit = FALSE)
}
\arguments{
\item{size}{A single numeric value of > 0 and <= 90.}

\item{fit}{A logical value indicating whether a bin size check should be performed to ensure that the
entire latitudinal range is covered (-90 to 90).}
}
\value{
A \code{data.frame} of latitudinal bins for a given size.
}
\description{
A function to generate latitudinal bins for a given bin size.
}
\examples{
lat_bins(size = 20)
lat_bins(size = 13, fit = TRUE)
}
