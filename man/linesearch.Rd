% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ArchitectureModel.R
\name{linesearch}
\alias{linesearch}
\title{Perform a line search of the Genomic Architecture parameter s}
\usage{
linesearch(
  nout,
  dat,
  bins = defaultbins(),
  thresholds = c(2.5, 3),
  range = c(-2, 2),
  sigma = NULL
)
}
\arguments{
\item{nout}{The number of values to evaluate.}

\item{dat}{A data frame of (f,beta) pairs}

\item{bins}{(Default: \code{\link{defaultbins}}) the histogram bins to use}

\item{thresholds}{(Default: \code{\link{defaultthresholds}}) the tail areas to compare (in s.d units).}

\item{range}{(Default: c(-2,2)) The range of s values to evaluate between}

\item{sigma}{(Default: NULL, meaning use the best value) Genomic architecture sigma}
}
\value{
a dataframe contining s and the score for that s.
}
\description{
Evaluate the loss function \code{\link{lossfunctionp}} for each value of s in \code{seq(range[1],range[2],length.out=nout)}.
}
\seealso{
\code{\link{fullsearch}} to infer a single p, or \code{\link{bootstrap}} to generate a distribution of estimates of s.
}
