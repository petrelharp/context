#!/usr/bin/Rscript
source("../context-inference-fns.R")

resultsfile <- commandArgs(TRUE)[1]
if (length(commandArgs(TRUE))<2) {
    residsfile <- gsub("-results.RData","-resids.tsv",resultsfile)
} else {
    residsfile <- commandArgs(TRUE)[2]
}

load(resultsfile)

resids <- listresids(counts,expected,file=residsfile)

## exploratory code
if (FALSE) {
    resultsfile <- "02-C-M_in_frame/win-2-2-2-1-results.RData"
    load(resultsfile)

    resids <- listresids(counts,expected,trim=0)
    old.mutpats <- names(point.estimate)[grepl("->",names(point.estimate))]

    load("genmatrices/genmatrix-6-none-0-2.RData")
    stopifnot( all( old.mutpats %in% mutnames(mutpats) ) )

    new.params <- which( !( mutnames(mutpats) %in% old.mutpats ) )
    adhoc <- countmuts(counts=counts,mutpats=mutpats,lwin=lwin)
    adhoc <- adhoc[1,]/adhoc[2,]
    base.mutrates <- numeric(length(mutpats))
    names(base.mutrates) <- mutnames(mutpats)
    base.mutrates[ match( old.mutpats, mutnames(mutpats) ) ] <- point.estimate[grepl("->",names(point.estimate))]
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
        } )


    # # not so good:
    # all.mutpats <- getmutpats(2)
    # mutpatlen <- sapply( all.mutpats, function (x) nchar(unlist(x)[1]) )
    # residmat <- 0 * counts 
    # stopifnot( all(toupper(resids$inpat) %in% rownames(residmat)) & all(resids$outpat %in% colnames(residmat)) )
    # residmat[cbind( match(toupper(resids$inpat),rownames(residmat)), match(resids$outpat,colnames(residmat)) )] <- resids$resid
    # resid.counts <- countmuts( residmat, all.mutpats, lwin=lwin )
    # plot(t(resid.counts[2:1,]),col=mutpatlen)
    # plot(resid.counts[2,]/resid.counts[1,], resid.counts[1,], col=mutpatlen)

}