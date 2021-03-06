% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ArchitectureModel.R
\name{lossfunctionp}
\alias{lossfunctionp}
\title{Loss function evaluation from architecture and data}
\usage{
lossfunctionp(s, ref)
}
\arguments{
\item{s}{Genomic Architecture parameter for estimation}

\item{ref}{A list containing: ref$dat: A data frame of (f,beta) pairs; ref$bins: the histogram bins to use; and ref$thresholds: the tail areas to compare (in s.d units).}
}
\value{
The loss as returned by \code{\link{lossfunction}}, run of the summary statistics returned by \code{\link{lossfunctionss}}.
}
\description{
A convenience function to calculate summary statistics and evaluate their loss. Not recommended for regular use, for which direct inference via \code{\link{lossfunctionss}} and \code{\link{lossfunction}} is recommended.
}
\seealso{
\code{\link{fullsearch}} for inference.
}
