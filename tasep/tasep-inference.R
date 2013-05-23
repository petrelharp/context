#!/usr/bin/Rscript

usage <- "\
Infer parameters from output of sim-tasep\
\
Usage:\
   Rscript tasep-inference.R simfile [win] [maxwin]\
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
        maxwin <- 6
    }
}


scriptdir <- "../"
source(paste(scriptdir,"codon-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))

if (interactive()) {
    # available simulated sequences
    simdir <- "tasep-sims"
    simfiles <- list.files(simdir,"*.RData",full.names=TRUE)
    siminfo <- do.call(rbind, lapply(simfiles, function (x) {
            load(x)
            data.frame(date=now,seqlen=seqlen, tlen=tlen, mutrate=mutrates[1],file=x,stringsAsFactors=FALSE)
        } ) )
    siminfo$meandist <- siminfo$tlen * siminfo$mutrate
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

lapply( 1:maxwin, function (thiswin) {

    basename <- paste(basedir,"/win-",thiswin,sep='')
    datafile <- paste( basename ,"-results.RData",sep='')
    resultsfile <- paste( basename ,"-results.tsv",sep='')
    plotfile <- paste( basename ,"plot",sep='')

    # Inference.
    lwin <- rwin <- thiswin
    winlen <- lwin+win+rwin

    genmatrix <- meangenmatrix( lwin=1, rwin=1, patlen=winlen, mutpats=mutpats, selpats=c(), boundary="wrap" )
    genmatrix@x <- update(genmatrix,mutrates,selcoef=numeric(0),Ne=1)
    projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin )
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

    counts <- list(
                counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[1]]$initseq, simseqs[[1]]$finalseq, lwin=lwin )
            )
    # allow nonmatches with small prob
    eps <- 1e-6

    # move from base frequencies (what we estimate) to pattern frequencies
    likfun <- function (mutrates) {
            # params are: Ne, branchlens[-1], mutrates, selcoef, basefreqs
            # this is collapsed transition matrix
            genmatrix@x <- update(genmatrix,mutrates,selcoef=numeric(0),Ne=1)
            subtransmatrix <- computetransmatrix( genmatrix, projmatrix )
            # return negative log-likelihood 
            (-1) * sum( counts[[1]] * log(eps+subtransmatrix) )
    }

    initpar <- c( runif( length(mutpats) ) * mean(mutrates) * tlen ) # random init
    truth <- c( mutrates * tlen )  # truth
    lbs <- c( 1e-3 )
    scalevec <- c( 1 )
    cheating.ans <- optim( par=truth, fn=likfun, method="L-BFGS-B", lower=lbs, control=list(trace=3, parscale=scalevec) )
    random.ans <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, control=list(trace=3, parscale=scalevec) )

    estimates <- data.frame( rbind(init=initpar, ans=random.ans$par, cheating=cheating.ans$par, truth=truth ) )
    colnames(estimates) <- "muttime"
    estimates$likfun <- apply( estimates, 1, likfun )
    write.table( estimates, file=resultsfile, quote=FALSE, sep="\t" )

    # look at observed and expected.
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix )
    true.expected <- lapply( counts, function (cnt) { rowSums(cnt) * subtransmatrix } )
    est.genmatrix <- genmatrix; est.genmatrix@x <- update(est.genmatrix,estimates["ans",1],selcoef=numeric(0),Ne=1)
    est.subtransmatrix <- computetransmatrix( est.genmatrix, projmatrix )
    est.expected <- lapply( counts, function (cnt) { rowSums(cnt) * est.subtransmatrix } )
    cheating.genmatrix <- genmatrix; cheating.genmatrix@x <- update(cheating.genmatrix,estimates["cheating",1],selcoef=numeric(0),Ne=1)
    cheating.subtransmatrix <- computetransmatrix( cheating.genmatrix, projmatrix )
    cheating.expected <- lapply( counts, function (cnt) { rowSums(cnt) * cheating.subtransmatrix } )
    init.counts <- rowSums(counttrans(rownames(genmatrix),colnames(projmatrix),simseqs= simseqs[[1]]))

    save( counts, genmatrix, subtransmatrix, truth, cheating.ans, random.ans, true.expected, est.expected, cheating.expected, file=datafile )

    pdf(file=paste(plotfile,"-1.pdf",sep=''),width=6, height=4, pointsize=10)
    lord <- order( true.expected[[1]][,1] )
    layout(1)
    plot( counts[[1]][lord,1], type='n', xaxt='n', xlab='', ylim=range(c(unlist(true.expected),unlist(lapply(counts,as.matrix)),unlist(est.expected))) )
    axis(1,at=1:nrow(counts[[1]]),labels=rownames(counts[[1]])[lord],las=3)
    for (k in 1:ncol(counts[[1]])) {
        for (j in 1:length(counts)) {
            points( counts[[j]][lord,k], pch=j )
            lines(true.expected[[j]][lord,k],col='red', lty=j)
            lines(est.expected[[j]][lord,k],col='green', lty=j, lwd=2)
            lines(cheating.expected[[j]][lord,k],col='blue',lty=j)
        }
        legend("topleft",legend=c("expected","estimated","cheating"),lty=1,col=c("red","green","blue"))
    }
    dev.off()

    pdf(file=paste(plotfile,"-2.pdf",sep=''),width=6, height=4, pointsize=10)
    layout(matrix(1:4,nrow=2))
    for (k in 1:4) {
        nonz <- pmax(true.expected[[1]][,k],counts[[1]][,k],est.expected[[1]][,k])>.5
        lord <- order( true.expected[[1]][,k] )[nonz]
        plot( counts[[1]][lord,k], type='n', xaxt='n', xlab='', main=paste(paste(rep(".",lwin),collapse=""),colnames(counts[[1]])[k],paste(rep(".",rwin),collapse=''),sep=''), ylim=range(c(unlist(true.expected),unlist(lapply(counts,as.matrix)),unlist(est.expected))) )
        axis(1,at=1:length(lord),labels=rownames(counts[[1]])[lord],las=3)
        for (j in 1:1) {
            points( counts[[j]][lord,k], pch=j )
            lines(true.expected[[j]][lord,k],col='red', lty=j)
            lines(est.expected[[j]][lord,k],col='green', lty=j, lwd=2)
            lines(cheating.expected[[j]][lord,k],col='blue',lty=j)
        }
        legend("topleft",legend=c("expected","estimated","cheating"),lty=1,col=c("red","green","blue"))
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
