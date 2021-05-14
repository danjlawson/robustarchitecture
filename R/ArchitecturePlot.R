##################
## plotting functions
#' @title Plot data annotated 
#'
#' @description Filter data to keep only the top k by proportion of variance explained, which must be provided in a column called 'propexplained'.
#' 
#' @param dat Data frame containing MAF and beta in first two columns.
#' @param res (Default: NULL)
#' @param mycols (Default: c(grey,blue,red)) three colours, one for each class of SNP
#' @param mar (Default: c(6,5,3,5)) Margins, passed to par(mar=mar) if not NA.
#' @param xlim (Default: c(0.01,0.5)) As provided to \code{plot}.
#' @param ylim (Default: NA, meaning use c(0,max(dat[,2])) As provided to \code{plot}.
#' @param xlab (Default: "MAF") As provided to \code{plot}.
#' @param ylab (Default: "Effect magnitude") As provided to \code{plot}.
#' @param log (Default: "") As provided to \code{plot}.
#' @param tx (Default: seq(0.005,0.5,length.out=100)) x values to evaluate the architecture curves at
#' @param main (Default: "Genomic Architecture structure") as provided to plot.
#' @param show.data (Default: TRUE) whether to plot the data.
#' @param show.bins (Default: TRUE) whether to plot the bins, if res provided.
#' @param show.res (Default: TRUE) whether to show the result proportions, if res provided.
#' @param data.pch (Default: 19) pch for data.
#' @param data.cex (Default: 0.5) cex for data.
#' @param bins.lty (Default: 2) line type of the bin boundaries (if res provided).
#' @param res.scale (Default: 1) vertical scaling of the results. Set to values < 1 to restrict the class proportion data into the bottom segment, which is often excluded from architecture plots.
#' @param res.lwd (Default: 1) Line width of the class proportion results.
#' @param res.col (Default: "#008800FF" which is dark green) colour for all results to be displayed.
#' @param res.shade (Default: "#88FF8844" which is pale green) colour for the theoretically expected variation in proportions per frequency bin.
#' @param res.lab (Default: "Extreme class proportion") axis label for the results. Written in res.col.
#' @param res.lab.line (Default: 3) line to write res.lab with mtext.
#' @param res.adj (Default:0.5) vertical location to write res.lab. Try setting to 0 if you have used res.scale<1.
#' @param ... Additional parameters to \code{plot}.
#' @return A dataframe the same shape as provided with only the top k SNPs.
#' @export
architectureplot=function(dat,
                          res=NULL,
                          mycols=c("#AAAAAA44","#0000FF44","#FF000044"),
                          mar=c(6,5,3,5),
                          xlim=c(0.01,0.5),
                          ylim=NULL,
                          xlab="MAF",
                          ylab="Effect magnitude",
                          log="",
                          tx=seq(0.005,0.5,length.out=100),
                          main="Genomic Architecture structure",
                          show.data=TRUE,
                          show.bins=TRUE,
                          show.res=TRUE,
                          data.pch=19,
                          data.cex=0.5,
                          bins.lty=2,
                          res.scale=1,
                          res.lwd=1,
                          res.col="#008800FF",
                          res.shade="#88FF8844",
                          res.lab="Extreme class proportion",
                          res.lab.line=3,
                          res.adj=0.5,
                          ...){
    if(any(!is.null(res))){
        dat=annotatearchitecture(dat[,1:2],res$s,res$thresholds,res$sigma)
    }else if (!("class"%in% colnames(dat))){
        dat[,"class"]=rep(1,dim(dat)[1])
        show.bins=FALSE
        warning("architectureplot works best if you call annotatearchitecture or provide results from which this information is available.")
    }
    if(all(is.null(ylim))) ylim=c(0,max(dat[,2]))
    if(!all(is.na(mar))) graphics::par(mar=mar)
    plot(abs(dat[,c(1,2)]),type="n",
         col=mycols[dat[,"class"]],
         xlab=xlab,
         ylab=ylab,
         main=main,
         log=log,xlim=xlim,ylim=ylim,...)
    if(show.bins) graphics::abline(v=res$bins,lty=bins.lty)
    if(show.data){
        graphics::points(abs(dat[,c(1,2)]),
                         col=mycols[dat[,"class"]],pch=data.pch,cex=data.cex)
    }
    if(any(!is.null(res))){
        if(show.res){
            fvals=res$ss$f-diff(res$bins)/2
            scalefn=function(x)x*diff(ylim)*res.scale+ylim[1]
            pscaled=scalefn(res$ss$p)
            qs=(1-(stats::pnorm(res$thresholds)))*2
            qrat=qs[2]/qs[1]
            qvals=sapply(1:dim(res$ss)[1],function(i){
                c(stats::qbinom((0.05/dim(res$ss)[1]),
                                res$ss[i,"count1"],qrat)/ res$ss[i,"count1"],
                  stats::qbinom((1-0.05/dim(res$ss)[1]),
                                res$ss[i,"count1"],qrat)/ res$ss[i,"count1"])
            })
            graphics::polygon(c(fvals,rev(fvals)),
                              scalefn(c(qvals[1,],rev(qvals[2,]))),
                              density=NA,border=NA,
                              col=res.shade)
            graphics::lines(fvals,pscaled,
                            col=res.col,
                            lwd=res.lwd)
            graphics::axis(4,at=scalefn(seq(0,1,by=0.2)),
                           labels=seq(0,1,by=0.2),
                           col=res.col,
                           col.axis=res.col,
                           las=2)
            graphics::mtext(res.lab,side=4,
                            line=res.lab.line,
                            col=res.col,
                            adj=res.adj)
        }
    }
}
