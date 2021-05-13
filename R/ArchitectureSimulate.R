###############################
#' @title Simulate from a Genomic Architecture
#'
#' @description
#' Simulate a genomic architecture from:
#'
#' sigma^2 (f (1-f))^(s)
#'
#' Where:
#' \code{f=unif(frange[1],frange[2])} and \code{beta=N(0,sigma*(f*(1-f))^(s/2))}.
#' 
#' @param N Number of SNPs to sample
#' @param sigma (default: 1) Baseline variance
#' @param s (default: -1) Genomic archictecture shape
#' @param frange (default: 0.01,0.5) range of frequencies to simulate for f
#' @return A data frame containing N rows and four columns: f, beta, varexplained (variance explained) and propexplained (proportion of variance explained).
#' @export
make_test_data=function(N,sigma=1,s=-1,frange=c(0.01,0.5)){
    f=stats::runif(N,frange[1],frange[2])
    beta=stats::rnorm(N,0,sigma*(f*(1-f))^(s/2))
    varexplained=beta^2*f*(1-f)
    data.frame(f=f,
               beta=beta,
               varexplained=varexplained,
               propexplained=varexplained/sum(varexplained))
}

#' @title Threshold a dataset based on frequency
#'
#' @description Assuming the first column is frequency, return only rows where f>=t.
#' 
#' @param testdata The data frame containing frequency in the first column.
#' @param t The desired threshold
#' @return A dataframe the same shape as provided with only common SNPs retained.
#' @export
thresh_data=function(testdata,t){
    testdata=testdata[testdata[,1]>=t,,drop=FALSE]
}

#' @title Retain only top SNPs by variance explained
#'
#' @description Filter data to keep only the top k by proportion of variance explained, which must be provided in a column called 'propexplained'.
#' 
#' @param testdata The data frame containing frequency in the first column.
#' @param k The number of SNPs to retain
#' @return A dataframe the same shape as provided with only the top k SNPs.
#' @export
top_data=function(testdata,k=100){
    to=order(testdata[,"propexplained"])
    testdata[utils::tail(to,n=k),,drop=FALSE]
}

#' @title Choose k random SNPs
#'
#' @description Keep at most k random SNPs from a dataframe.
#' 
#' @param testdata The data frame containing frequency in the first column.
#' @param k The number of SNPs to retain
#' @return A dataframe the same shape as provided with only min(k,dim(testdata)[1]) SNPs.
#' @export
random_data=function(testdata,k=100){
    ts=sample(1:dim(testdata)[1],k)
    testdata[ts,,drop=FALSE]
}
