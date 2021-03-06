% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gcImp.R
\name{gc.as.mids}
\alias{gc.as.mids}
\title{Convert Gaussian copula to mids}
\usage{
gc.as.mids(out)
}
\arguments{
\item{out}{A gcImp calss object, usually the output from \code{gcImp}}
}
\value{
An object of class "mira," see ?mice::mira for help.
}
\description{
\code{gc.as.mids} converts a gcImp object into a "mira" format. 
The mira format can then be used by the mice package to analyze imputed data.
}
