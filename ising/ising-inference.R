#!/usr/bin/Rscript

usage <- "\
Infer parameters from output of sim-ising\
\
Usage:\
   Rscript ising-inference.R simfile [win] [maxwin]\
"

args <- commandArgs(TRUE)
if (length(args)<1) {
    stop(usage)
} else {
    thissim <- args[1]
    if (length(args)>=2) {
        win <- as.numeric(args[2])
    } else {
        win <- 2
    }
    if (length(args)>=3) {
        maxwin <- as.numeric(args[3])
    } else {
        maxwin <- 4
    }
}


scriptdir <- "../"
source(paste(scriptdir,"codon-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))

if (interactive()) {

    # available simulated sequences
    simdir <- "ising-sims"
    simfiles <- list.files(simdir,"*.RData",full.names=TRUE)
    siminfo <- do.call(rbind, lapply(simfiles, function (x) {
            load(x)
            y <- do.call( data.frame, c(list(date=now,seqlen=seqlen, tlen=tlen, file=x), as.list(mutrates), as.list(selcoef), stringsAsFactors=FALSE) )
            colnames(y) <- c( colnames(y)[1:4], paste("mutrate",seq_along(mutrates),sep=''), paste("selcoef",seq_along(selcoef),sep='') )
            y
        } ) )
    siminfo$meandist <- siminfo$tlen * colSums(siminfo[,grep("mutrate",colnames(siminfo)),drop=FALSE])
    siminfo <- siminfo[ order(siminfo$meandist), ]

    win <- 2

    # pick one
    thissim <- siminfo[1,"file"]
}

load(thissim)
basedir <- gsub(".RData","",thissim,fixed=TRUE)
if (!file.exists(basedir)) {
    dir.create(basedir)
}

# set-up
bases <- c("X","O")

lapply( 2:maxwin, function (thiswin) {

    basename <- paste(basedir,"/win-",thiswin,sep='')
    datafile <- paste( basename ,"-results.RData",sep='')
    resultsfile <- paste( basename ,"-results.tsv",sep='')
    plotfile <- paste( basename ,"-plot",sep='')

    # Inference.
    lwin <- rwin <- thiswin
    winlen <- lwin+win+rwin

    genmatrix <- meangenmatrix( lwin=1, rwin=1, patlen=winlen, mutpats=mutpats, selpats=selpats, mutrates=mutrates*tlen, selcoef=selcoef, boundary="none" )
    projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin )
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

    counts <- list(
                counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[1]]$initseq, simseqs[[1]]$finalseq, lwin=lwin )
            )

    nmuts <- length(mutpats)
    nsel <- length(selpats)
    # move from base frequencies (what we estimate) to pattern frequencies
    likfun <- function (params) {
            # params are: mutrates, selcoef
            # this is collapsed transition matrix
            mutrates <- params[1:nmuts]
            selcoef <- params[nmuts+1:nsel]
            genmatrix@x <- update(genmatrix,mutrates=mutrates,selcoef=selcoef)
            subtransmatrix <- computetransmatrix( genmatrix, projmatrix )
            # return negative log-likelihood 
            (-1) * sum( counts[[1]] * log(subtransmatrix) )
    }

    initpar <- c( 2 * runif( length(mutpats) ) * mean(mutrates) * tlen, 2 * runif( length(selpats) ) * mean(selcoef) ) # random init
    truth <- c( mutrates * tlen, selcoef )  # truth
    lbs <- c( rep(1e-6,nmuts), rep(-Inf,nsel) )
    cheating.ans <- optim( par=truth, fn=likfun, method="L-BFGS-B", lower=lbs, control=list(trace=3) )
    random.ans <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, control=list(trace=3) )

    estimates <- data.frame( rbind(init=initpar, ans=random.ans$par, cheating=cheating.ans$par, truth=truth ) )
    colnames(estimates) <- c( paste("muttime",seq_along(mutrates),sep=''), paste("selcoef",seq_along(selcoef),sep='') )
    estimates$likfun <- apply( estimates, 1, likfun )
    write.table( estimates, file=resultsfile, quote=FALSE, sep="\t" )

    # look at observed/expected counts
    all.expected <- lapply( 1:nrow(estimates), function (k) {
                x <- unlist(estimates[k,])
                predictcounts( win, lwin, rwin, initcounts=rowSums(counts[[1]]), mutrates=x[1:nmuts], selcoef=x[nmuts+(1:nsel)], genmatrix=genmatrix, projmatrix=projmatrix )
        } )
    names(all.expected) <- rownames(estimates)

    # look at observed/expected counts in smaller windows
    cwin <- 2
    subcounts <- projectcounts( lwin=lwin, countwin=cwin, lcountwin=0, rcountwin=0, counts=counts[[1]] )
    all.subexpected <- lapply( all.expected, function (x) { projectcounts( lwin=lwin, countwin=cwin, lcountwin=0, rcountwin=0, counts=x ) } )

    # # look at observed and expected.
    # subtransmatrix <- computetransmatrix( genmatrix, projmatrix )
    # true.expected <- lapply( counts, function (cnt) { rowSums(cnt) * subtransmatrix } )
    # est.genmatrix <- genmatrix; est.genmatrix@x <- update(est.genmatrix,estimates["ans",1],selcoef=numeric(0),Ne=1)
    # est.subtransmatrix <- computetransmatrix( est.genmatrix, projmatrix )
    # est.expected <- lapply( counts, function (cnt) { rowSums(cnt) * est.subtransmatrix } )
    # cheating.genmatrix <- genmatrix; cheating.genmatrix@x <- update(cheating.genmatrix,estimates["cheating",1],selcoef=numeric(0),Ne=1)
    # cheating.subtransmatrix <- computetransmatrix( cheating.genmatrix, projmatrix )
    # cheating.expected <- lapply( counts, function (cnt) { rowSums(cnt) * cheating.subtransmatrix } )
    # init.counts <- rowSums(counttrans(rownames(genmatrix),colnames(projmatrix),simseqs= simseqs[[1]]))
    # save( counts, genmatrix, subtransmatrix, truth, cheating.ans, random.ans, true.expected, est.expected, cheating.expected, file=datafile )

    save( counts, genmatrix, subtransmatrix, truth, cheating.ans, random.ans, all.expected, cwin, subcounts, all.subexpected, file=datafile )

    pdf(file=paste(plotfile,"-1.pdf",sep=''),width=6, height=4, pointsize=10)
    lord <- order( all.expected[["truth"]][[1]][,1] )
    layout(1)
    plot( counts[[1]][lord,1], type='n', xaxt='n', xlab='', ylim=range(c(unlist(all.expected[["truth"]]),unlist(lapply(counts,as.matrix)),unlist(all.expected[["ans"]]))) )
    axis(1,at=1:nrow(counts[[1]]),labels=rownames(counts[[1]])[lord],las=3)
    for (k in 1:ncol(counts[[1]])) {
        for (j in 1:length(counts)) {
            points( counts[[j]][lord,k], pch=j )
            lines(all.expected[["truth"]][[j]][lord,k],col='red', lty=j)
            lines(all.expected[["ans"]][[j]][lord,k],col='green', lty=j, lwd=2)
            lines(all.expected[["cheating"]][[j]][lord,k],col='blue',lty=j)
        }
        legend("topleft",legend=c("expected","estimated","cheating"),lty=1,col=c("red","green","blue"))
    }
    dev.off()

    pdf(file=paste(plotfile,"-2.pdf",sep=''),width=6, height=4, pointsize=10)
    layout(matrix(1:4,nrow=2))
    cols <- rainbow(2+length(all.expected))[1:length(all.expected)]
    for (k in 1:4) {
        lord <- order( all.subexpected[["truth"]][,k] )
        plot( subcounts[lord,k], xaxt='n', xlab='', main=colnames(subcounts)[k] )
        axis(1,at=1:nrow(subcounts),labels=rownames(subcounts)[lord],las=3)
        invisible( lapply(seq_along(all.subexpected),function(j) { lines(all.subexpected[[j]][lord,k],col=cols[j]) } ) )
        legend("topleft",legend=names(all.subexpected),lty=1,col=cols)
    }
    dev.off()

    pdf(file=paste(plotfile,"-3.pdf",sep=''),width=6, height=4, pointsize=10)
    layout(1)
    plot( as.vector(true.expected[[1]]), as.vector(counts[[1]]), log='xy' )
    abline(0,1)
    points(as.vector(true.expected[[1]]), as.vector(est.expected[[1]]), col='red', pch=20 )
    points(as.vector(true.expected[[1]]), as.vector(cheating.expected[[1]]), col='green', pch=20, cex=.5)
    dev.off()


} )

print(format(Sys.time(),"%Y-%m-%d-%H-%M"))
