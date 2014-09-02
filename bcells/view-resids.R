#!/usr/bin/Rscript
source("../context-inference-fns.R")

resultsfile <- commandArgs(TRUE)[1]
if (length(commandArgs(TRUE))<2) {
    residsfile <- gsub("-results.RData","-resids.tsv",resultsfile)
} else {
    residsfile <- commandArgs(TRUE)[2]
}

load(resultsfile)
load(opt$gmfile)

projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin )
pcounts <- function (params) { predictcounts(shortwin, leftwin, rightwin, initcounts=rowSums(counts), mutrates=params[1:nmuts], selcoef=numeric(0), scale=params[length(params)], genmatrix=genmatrix, projmatrix=projmatrix, time="fixed" ) }
expected <- pcounts(point.estimate)

resids <- listresids(counts,expected,file=residsfile)

## exploratory code
if (FALSE) {
    resultsfile <- "02-C-M_in_frame/win-2-2-2-1-results.RData"
    load(resultsfile)

    resultsfile <- "02-C-M_in_frame/02-C-M_in_frame.tuples.6.2.counts-genmatrix-6-model-02-112777-results.RData"
    load(resultsfile)
    load(opt$gmfile)

    projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin )
    pcounts <- function (params) { predictcounts(shortwin, leftwin, rightwin, initcounts=rowSums(counts), mutrates=params[1:nmuts], selcoef=numeric(0), scale=params[length(params)], genmatrix=genmatrix, projmatrix=projmatrix, time="fixed" ) }
    expected <- pcounts(point.estimate)

    resids <- listresids(counts,expected,trim=0)
    old.mutpats <- names(point.estimate)[grepl("->",names(point.estimate))]

    load("genmatrices/genmatrix-6-dualbases.RData")

    mutpats <- genmatrix@mutpats

    stopifnot( all( old.mutpats %in% mutnames(mutpats) ) )

    new.params <- which( !( mutnames(mutpats) %in% old.mutpats ) )
    adhoc <- countmuts(counts=counts,mutpats=mutpats,leftwin=leftwin)
    adhoc <- adhoc[1,]/adhoc[2,]
    base.mutrates <- numeric(length(mutpats))
    names(base.mutrates) <- mutnames(mutpats)
    base.mutrates[ match( old.mutpats, mutnames(mutpats) ) ] <- point.estimate[grepl("->",names(point.estimate))]
    require(parallel)
    marginal.params <- mclapply( new.params, function (k) {
                # (quasi)-likelihood function using all counts -- multinomial
                likfun <- function (params) {
                    # params are: mutrates*tlen, shape
                    base.mutrates[k] <- params
                    genmatrix@x <- update(genmatrix,mutrates=base.mutrates,selcoef=numeric(0))
                    # this is collapsed transition matrix
                    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, tlen=1, time="fixed") # shape=params[length(params)], time="gamma" )
                    # return NEGATIVE log-likelihood 
                    ans <- (-1) * sum( counts * log(subtransmatrix) ) + (params[length(params)]-1)^2
                    if (!is.finite(ans)) { print(params) }
                    return(ans)
                }
                # mrates <- seq(adhoc[k]/5,adhoc[k]*20,length.out=10)
                # likvals <- sapply(mrates, likfun)
                optimize( likfun, interval=c(0,10*adhoc[k]) )
        }, mc.cores=8 )
    # save( marginal.params, file="temp.RData")
    names(marginal.params) <- mutnames(mutpats)[new.params]


    # not so good:
    all.mutpats <- getmutpats(2)
    mutpatlen <- sapply( all.mutpats, function (x) nchar(unlist(x)[1]) )
    residmat <- 0 * counts 
    stopifnot( all(toupper(resids$inpat) %in% rownames(residmat)) & all(resids$outpat %in% colnames(residmat)) )
    residmat[cbind( match(toupper(resids$inpat),rownames(residmat)), match(resids$outpat,colnames(residmat)) )] <- resids$resid
    resid.counts <- countmuts( residmat, all.mutpats, leftwin=leftwin )
    plot(t(resid.counts[2:1,]),col=mutpatlen)
    plot(resid.counts[2,]/resid.counts[1,], resid.counts[1,], col=mutpatlen)

}
