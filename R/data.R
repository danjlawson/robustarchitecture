#' @title Simulated summary data with Genomic Architecture
#'
#' @description
#' Simulated data using s=-0.5 and only the top 10\% (10k) SNPs
#' retained by variance explained. This masking leads to a complex
#' inference problem typical of real data, but is not a problem for
#' the \code{link{genomicarchitecture}} approach.
#'
#' @docType data
#' @name topdata
#' @format An data frame containing 4 columns, "f","beta","varexplained","propexplained". Only the first two are required for inference.
#'
#' @keywords datasets
#'
#' @seealso \code{\link{make_test_data}} to generate your own simulated data, and \code{\link{top_data}} to threshold it. \code{\link{fullsearch}} to perform inference and \code{\link{architectureplot}} to plot the results.
#' @examples
#' \dontrun{
#' data(topdata)
#' topres=fullsearch(dat=topdata)
#' architectureplot(topdata,topres)
#' }

NULL
