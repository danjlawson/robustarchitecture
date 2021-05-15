###############################
#' @title Effect size of a SNP
#'
#' @description
#' Get the mean effect of a genomic architecture model, i.e.
#' sigma (f (1-f))^(s/2)
#' 
#' @param f The freqeuncy of the SNP
#' @param sigma The standard deviation 
#' @param s The Selection parameter
#' @return The variance explained
#' @export
getbeta=function(f,sigma,s){
  sigma*(f*(1-f))^(s/2)
}

###############################
#' @title Default bin locations for general purpose usage
#'
#' @description
#' Histogram bins that should be ok for most use; these are:
#' 
#' seq(0.05,0.5,by=0.025)
#' 
#' @return a vector of bin separators
#' @export
defaultbins=function(){
    seq(0.05,0.5,by=0.025)
}

###############################
#' @title Default thresholds
#'
#' @description
#' Thresholds that should be accurate if they can be used. These are:
#' 
#' c(2.5,3)
#'
#' in standard deviation units.
#'
#' If you don't have data that has this large a variance, you may choose a smaller
#' threshold pair; c(1,2) is about as small as you should go.
#' 
#' @return A pair of threshold values, the first for data to be considered, the second to be predicted.
#' @export
defaultthresholds=function(){
    c(2.5,3)
}

###############################
#' @title Get the best estimate of sigma conditional on s 
#'
#' @description
#' Estimate sigma using the provided SNPs and a given s. This is done by ensuring that the whole distribution is scaled to obtain the theoretical proportion of SNPs above the second threshold, given the first, thereby ensuring that if s is correct then there should be no relationship between MAF and whether beta is over each threshold.
#' 
#' This uses a binary search of the possible thresholds that impact the choice, and is therefore logarithmic in the number of SNPs.
#' 
#' @param ratio The ratio of beta/betahat(f,s,sigma=1)
#' @param thresholds The architecture thresholds to be used in s.d. units.
#' @return The estimate of sigma
#' @export
getsigma=function(ratio,thresholds){
    qs=(1-(stats::pnorm(thresholds)))*2
    qrat=qs[2]/qs[1]
    psigma=function(tr,x,thresholds){
       sum(tr>= x* thresholds[2])/sum(tr>= x* thresholds[1])
    }
    tr=sort(ratio)
    ## Binary search for the best value
    li=1; ui=length(tr)-1; mi=floor(length(tr)/2) ## Indices
    ls=(tr[li])/thresholds[2];
    us=(tr[ui])/thresholds[2];
    ms=(tr[mi])/thresholds[2]; ## Sigmas
    lv=psigma(tr,ls,thresholds);uv=psigma(tr,us,thresholds);mv=psigma(tr,ms,thresholds); ## values (proportions)
    while((ui-li)>1){
        if(mv<qrat){
            uv=uv;us=ms;ui=mi
            mi=li+floor((ui-li)/2)
            ms=(tr[mi])/thresholds[2];
            mv=psigma(tr,ms,thresholds);
        }else{
            lv=mv;ls=ms;li=mi
            mi=li+floor((ui-li)/2)
            ms=(tr[mi])/thresholds[2];
            mv=psigma(tr,ms,thresholds);
        }
    }
    ms
}

#' @title Annotate SNPs in a genetic architecture according to a model
#'
#' @description
#' Compute the effect size prediction for a specific genomic architecture model, thereby annotating SNPs into classes according to how they fit with respect to the given thresholds.
#'
#' Recall that we are not expecting individual SNPs to match their prediction, only the variance conditional on MAF to do so.
#' 
#' @param dat data frame with the first two columns being "f" and "beta"
#' @param s The architecture model value for s
#' @param thresholds (Default: \code{\link{defaultthresholds}}) The architecture thresholds to be used in s.d. units.
#' @param sigma (Default: NULL) The sigma to use for the data. Default: use that estimated from the data.
#' @return The class of each datapoint
#' @export

annotatearchitecture=function(dat,s,
                              thresholds=c(2.5,3),sigma=NULL){
    dat<-as.data.frame(dat[,1:2])
    colnames(dat)=c("f","beta")
    if(is.null(sigma)) {
        betapred=getbeta(dat[,"f"],1,s)
        sigma=getsigma(ratio=dat[,"beta"]/betapred,thresholds)
        dat[,"betapred"] = betapred*sigma
    }else{
        dat[,"betapred"] = getbeta(dat[,"f"],sigma,s)
    }
    dat[,"betapred1"]=dat[,"betapred"]*thresholds[1]
    dat[,"betapred2"]=dat[,"betapred"]*thresholds[2]
    dat[,"class"] =  1+ (abs(dat$beta)>dat[,"betapred1"])+(abs(dat$beta)>dat[,"betapred2"])
    dat
}
###############################
#' @title Calculation of histogram-based summary statistics of a distribution
#'
#' @description
#'
#' Compute summary statistics from a dataset of MAF,beta pairs, by binning data into histograms and comparing whether the theoretical proportion of SNPs above a threshold matches the observed proportion.
#'
#' Uses \code{\link{annotatearchitecture}} to annotate each SNP by with its prediction, and \code{\link{getsigma}} to estimate sigma given s and the dataset.
#' 
#' @param dat A dataframe containing summary statistics from a GWAS/BayesS/LDAK analysis. In its first two columns we need: MAF (f) and beta.
#' @param s The selection coefficient assumed in the architecture model
#' @param bins (Default: bins=\code{\link{defaultbins}} The histogram bin boundaries assumed for creating summary statistcs, which need to be large enough to contain enough samples and small enough to retain resolution.
#' @param thresholds (Default: \code{\link{defaultthresholds}}) Standard deviation units for the summary statistic generation. Only SNPs above the first are used for model learning. The proportion above the second is predicted from the number above the first.
#' @param sigma (Default: NULL meaning estimate from data) The scale parameter of the Genomic Architecture. Do not provide this with any search procedure as the results can be poor.
#' @return A data frame containing: bins, count0 (The number of SNPs in the bin), count (the number of SNPs above the first threshold), count2 (the number of SNPs above the second threshold)
#' @export
#' @seealso \code{\link{lossfunction}} to turn this summary into a loss, and \code{\link{fullsearch}} to search over s.
lossfunctionss=function(dat,s,
                        bins=defaultbins(),
                        thresholds=defaultthresholds(),
                        sigma=NULL){
    
    if(is.null(sigma)) {
        sigma=getsigma(dat[,2]/getbeta(dat[,1],1,s),
                       thresholds)
    }
    dat=annotatearchitecture(dat,s,thresholds,sigma)
    
    ## Remove SNPs that are too rare to be included in the base class
    if(any(dat$class==1)){
        fmin=min(dat[dat$class==1,1])
        dat=dat[dat$f>=fmin,,drop=FALSE]
    }
    
    ## Iterate through SNPs in order of frequency, and assign them to bins
    dat=dat[order(dat[,1]),]
    if(dim(dat)[1]==0) return(ret)
    tbins=sapply(dat[,1],function(x)which.min(x>bins))
    ubins=sort(unique(tbins))
    ret=t(sapply(1:length(bins),function(x){
        tmp=dat[tbins==x,]
        c(f=bins[x],
          count0=dim(tmp)[1]-sum(abs(tmp[,"beta"])>tmp[,"betapred1"]),
          count1=sum(abs(tmp[,"beta"])>tmp[,"betapred1"]),
          count2=sum(abs(tmp[,"beta"])>tmp[,"betapred2"]))
    }))
    ret=ret[-1,,drop=FALSE]
    ## Expand the results with empirical frequencies
    ret=cbind(ret,"p"=ret[,"count2"]/ret[,"count1"])
    ret[is.na(ret[,"p"]),"p"]=0
    ret=cbind(ret,"sigma"=sigma)

    ## Expand the results with per-bin scores
    ret=as.data.frame(ret)
    ret[,"loss"]=stats::dbinom(ret[,"count2"],ret[,"count1"],
                        stats::weighted.mean(ret[,"count2"]/ret[,"count1"],ret[,"count1"]),
                        log=T)/(ret[,"count1"])

    ret
}

## ###############################
#' @title Compute the loss from summary statistics
#'
#' @description
#'
#' Takes the loss function summary statistics table as returned from \code{\link{lossfunctionss}} and evaluates the loss, defined by the mean of the loss function for all SNPs across all MAF bins.
#' 
#' @param ss The ss as returned by \code{\link{lossfunctionss}}.
#' @return The score of the summary statistics as a loss (higher is better)
#' @export
#' @seealso \code{\link{lossfunctionss}} to generate summary statistics, and \code{\link{fullsearch}} to search over s.
lossfunction=function(ss){

    ss[ss[,"count0"]==0,c("count1","count2","loss")]=0
    
    ## if((any(ss[,"count2"]<=1))||(any(ss[,"count1"]<=1))) {
    ##     warning("Some bins are empty! This typically means that the architecture model is poorly fitting. Consider increasing thresholds or providing a smaller range of bins.")
    ##     return(-Inf) # terrible score: no data
    ## }
    if(sum(ss[,"count2"])<=1) return(-Inf) # terrible score: no data
    if(sum(ss[,"count1"])<=1) return(-Inf) # terrible score: no data
    if(sum(ss[,"count0"]>0)<=4) return(-Inf) # terrible score: no data

    ss[,"wt"]=ss[,"count1"] #+ss[,"count0"] 
    if(sum(ss[,"wt"])==0) return(-Inf)
    ##    stats::weighted.mean(ss[,"loss"],ss[,"wt"])
    sum(ss[,"loss"])
}

###############################
#' @title Loss function evaluation from architecture and data
#'
#' @description
#'
#' A convenience function to calculate summary statistics and evaluate their loss. Not recommended for regular use, for which direct inference via \code{\link{lossfunctionss}} and \code{\link{lossfunction}} is recommended.
#' 
#' @param s Genomic Architecture parameter for estimation
#' @param ref A list containing: ref$dat: A data frame of (f,beta) pairs; ref$bins: the histogram bins to use; and ref$thresholds: the tail areas to compare (in s.d units).
#' @return The loss as returned by \code{\link{lossfunction}}, run of the summary statistics returned by \code{\link{lossfunctionss}}.
#' @export
#' @seealso \code{\link{fullsearch}} for inference.
lossfunctionp=function(s,ref){
    tss=lossfunctionss(ref$dat,s,ref$bins,ref$thresholds,ref$sigma)
    lossfunction(tss)
}

###############################
#' @title Details of an architecture
#'
#' @description
#'
#' Details of an architecture, which may be inferred or defined manually.
#' 
#' @param dat A data frame of (f,beta) pairs
#' @param s The architecture parameter to use
#' @param thresholds (Default :\code{\link{defaultthresholds}}) the tail areas to compare (in s.d units).
#' @param bins (Default: \code{\link{defaultbins}}) the histogram bins to use
#' @param range (Default: c(-2,2) The range of s that could be considered with finite values (if known).
#' @param range0 (Default: c(-2,2) The range of s that was chosen by the user.
#' @param sigma (Default: NULL, meaning use the best value) Genomic architecture sigma. Don't change unless you know what you are doing...
#' @export
#' @seealso \code{\link{architectureplot}} for visualisation.
detail=function(dat,
                s,
                thresholds=defaultthresholds(),
                bins=defaultbins(),
                range=c(-2,2),
                range0=c(-2,2),
                sigma=NULL) {
    ret=list()
    ret$dat=dat
    if(!all(is.null(sigma))){
        ret$sigma=getsigma(dat[,2]/getbeta(dat[,1],1,s),
                           thresholds)
    }else{
        ret$sigma=sigma
    }
    ret$ss=lossfunctionss(dat,s,bins,thresholds,sigma)
    ret$thresholds=thresholds
    ret$bins=bins
    ret$range0=range0
    ret$range=range
    ret$s=s
    ret
}


#' @title Search s for the best parameter
#'
#' @description
#'
#' Optimize \code{\link{lossfunctionp}}. For basic usage, just give \code{dat}
#' and it should do everything else for you. You sometimes have to
#' change \code{thresholds} to fit the shape of your data; check this with
#' \code{\link{architectureplot}}.
#'
#' This is a 1-D non-convex optimization problem. The first challenge is that many
#' parameters don't have a finite loss, so we first do a \code{\link{linesearch}}
#' to identify a valid range.
#'
#' We then either do a sequence of linesearches to arrive at a region, or go straight for a
#' convex optimizer ("optimize" or "Brent" from \code{\link[stats]{optim}}). Line searches
#' are slower but cope better with small datasets.
#' 
#' @param dat A data frame of (f,beta) pairs
#' @param bins (Default: \code{\link{defaultbins}}) the histogram bins to use
#' @param thresholds (Default :\code{\link{defaultthresholds}}) the tail areas to compare (in s.d units).
#' @param sigma (Default: NULL, meaning use the best value) Genomic architecture sigma. Don't change unless you know what you are doing...
#' @param method (Default: c("linesearch","optimize")) The set of all methods to use (out of "linesearch", "optimize" and "Brent"). Good values are c("linesearch","optimize"), which tries to find global maxima, and "Brent" which tries to find a local maxima only, but is slightly faster. Only specify up to one of "Brent" and "optimize".
#' @param nline (Defaul: 21) Number of linesearch parameters to evaluate in .
#' @param nlinesearch (Defaul: 1) Number of linesearch parameters to evaluate, controlling the resolution near the edges of the valid range.
#' @param range (Default: c(-2,2)) The lower and upper bounds provided to \code{\link{linesearch}} to find a region that returns a finite score.
#' @return A list as returned by the optimiser, with additionally \code{dat}, \code{thresholds}, \code{bins}, \code{range} and \code{sigma} as provided, and \code{s} as inferred.
#' @export
#' @seealso \code{\link{linesearch}} to explore a range of s, or \code{\link{bootstrap}} to generate a distribution of estimates of s.
fullsearch=function(dat,
                    bins=defaultbins(),
                    thresholds=c(2.5,3),sigma=NULL,
                    method=c("linesearch","optimize"),
                    nlinesearch=2,
                    nline=21,
                    range=c(-2,2)
                    ){
    ## First we identify a valid region with linesearch
    range0=range
    testline=suppressWarnings(
        linesearch(nline,dat,bins=bins,
                   thresholds=thresholds,range=range,sigma=sigma)
        )
    testline[!is.finite(testline[,2]),2]<-NA
    range=range(stats::na.omit(testline)[,1])

    ## Second we apply this iteratively
    if("linesearch"%in% method){
        rangetest=range
        for(i in 1:nlinesearch){
            test=linesearch(nline,dat,bins=bins,
                            thresholds=thresholds,range=rangetest,sigma=sigma)
            rangetest=test[order(test[,2],decreasing=T)[1:2],1]
        }
        range=rangetest
        ret=list(
            par=test[order(test[,2],decreasing=T)[1],1],
            value=test[order(test[,2],decreasing=T)[1],2]
        )
    }
    ## Third we optimize locally
    if("optimize" %in% method){
        ret=stats::optimize(lossfunctionp,range,maximum=TRUE,
                     ref=list(dat=dat,bins=bins,thresholds=thresholds,sigma=sigma))
        ret$par=ret$maximum
        ret$value=ret$objective
    }else if("Brent" %in% method){
        s=sum(range)/2 # Initial guess
        ret=stats::optim(s,
          lossfunctionp,
          ref=list(dat=dat,bins=bins,thresholds=thresholds,sigma=sigma),
          method="Brent",
##          method="L-BFGS-B",
          lower=range[1],upper=range[2],
          control=list(fnscale=-1))
    }else{
        stop("Invalid method")
    }
    if(!is.finite(ret$value)){
        ret$s<-ret$par<-NA
        warning("No values of s result in valid data! You may need to reduce your bin range or increase your thresholds.")
    }else{
        ret=detail(dat,ret$par,thresholds,bins,range,range0,sigma)
    }
    ret
}
#' @title Obtain a bootstrap sample of Genomic Architecture s
#'
#' @description
#'
#' Resample the data \code{dat} with replacement a given number \code{nbs} times, for each running \code{\link{fullsearch}} and reporting the inferred parameter.
#'
#' @param nbs The number of bootstrap samples to obtain. We recommend 20 for exploration, 100 for estimating standard deviations and more for estimating the shape of the distribution.
#' @param dat A data frame of (f,beta) pairs
#' @param bins (default: \code{\link{defaultbins}}) the histogram bins to use
#' @param thresholds (Default: \code{\link{defaultthresholds}}) the tail areas to compare (in s.d units).
#' @param sigma (Default: NULL, meaning use the best value) Genomic architecture sigma
#' @param verbose (Default: TRUE) Whether to print iteration progress.
#' @param ... Additional parameters to \code{\link{fullsearch}}.
#' @return A vector of length \code{nbs} of estimates of \code{s}.
#' @export
#' @seealso \code{\link{linesearch}} to explore a range of s, or \code{\link{fullsearch}} to search a single dataset.
bootstrap=function(nbs,dat,bins=defaultbins(),
                   thresholds=c(2.5,3),sigma=NULL,verbose=TRUE,...){
    ret=sapply(1:nbs,function(i){
        if(verbose)cat(paste("Bootstrap",i,"\n"))
        tdat=sample(1:dim(dat)[1],dim(dat)[1],replace=TRUE)
        tres=suppressWarnings(
            fullsearch(dat=dat[tdat,,drop=FALSE],
                        bins=bins,
                        thresholds=thresholds,
                        sigma=sigma,
                        ...))
        tres$par
    })
    if(any(is.na(ret))) warning("Some values of s result in invalid data. You may need to reduce your bin range or increase your thresholds.")
    ret
}
#' @title Perform a line search of the Genomic Architecture parameter s
#'
#' @description
#'
#' Evaluate the loss function \code{\link{lossfunctionp}} for each value of s in \code{seq(range[1],range[2],length.out=nout)}.
#'
#' @param nout The number of values to evaluate.
#' @param dat A data frame of (f,beta) pairs
#' @param bins (Default: \code{\link{defaultbins}}) the histogram bins to use
#' @param thresholds (Default: \code{\link{defaultthresholds}}) the tail areas to compare (in s.d units).
#' @param range (Default: c(-2,2)) The range of s values to evaluate between 
#' @param sigma (Default: NULL, meaning use the best value) Genomic architecture sigma
#' @return a dataframe contining s and the score for that s.
#' @export
#' @seealso \code{\link{fullsearch}} to infer a single p, or \code{\link{bootstrap}} to generate a distribution of estimates of s.
linesearch=function(nout,dat,bins=defaultbins(),thresholds=c(2.5,3),range=c(-2,2),sigma=NULL){
    tx=seq(range[1],range[2],length.out=nout)
    ty=sapply(tx,lossfunctionp,
              ref=list(dat=dat,
                       bins=bins,
                       thresholds=thresholds,
                       sigma=sigma))
    if(any(!is.finite(ty))) warning("Some values of s result in invalid data. Check whether the plausible range matches. You may need to reduce your bin range or increase your thresholds.")
    data.frame(s=tx,score=ty)
}

