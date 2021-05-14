#' @title robustarchitecture : Robustly Learn Genomic Architecture from Summary Statistics
#'
#' @description
#'
#' One key aspect of Genomic Architecture is the relationship between SNP frequency and the effect size each SNP has on a trait. This has implications for whether the trait is has been subject to natural selection. It also controls whether it will be easy or complex to predict from the genome when learned via Genome Wide Association Studies.
#'
#' This package uses a robust statistical procedure to learn this from Summary Statistic data, i.e. (MAF,beta) pairs where MAF is minor allele frequency and beta is the effect size. Only the SNPs with the largest effect are required to learn this.
#' 
#' The code for robustarchitecture is at \url{https://github.com/danjlawson/robustarchitecture}.
#'
#' @references To follow
#' @importFrom methods is
"_PACKAGE"
#> [1] "_PACKAGE"
