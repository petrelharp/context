#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Infer parameters from output of sim-ising.R .\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file containing simulation." ),
        make_option( c("-w","--win"), type="integer", default=2, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--lwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-n","--nbatches"), type="integer", default=1000, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-m","--mmean"), type="double", default=1, help="Prior mean on mutation rates. [default \"%default\"]" ),
        make_option( c("-v","--svar"), type="double", default=1, help="Prior variance on selection coefficents. [default \"%default\"]" ),
        make_option( c("-c","--continue"), action="store_true", default=FALSE, help="Continue with previous MCMC run?" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)

if (is.null(infile) | is.null(nbatches)) { cat("Run\n  ising-inference.R -h\n for help.") }

scriptdir <- "../"
source(paste(scriptdir,"codon-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))

require(mcmc)

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
    infile <- siminfo[1,"file"]
}

load(infile)
basedir <- gsub(".RData","",infile,fixed=TRUE)
if (!file.exists(basedir)) {
    dir.create(basedir)
}

# set-up
bases <- c("X","O")


basename <- paste(basedir,"/win-",lwin,"-",win,"-",rwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
resultsfile <- paste( basename ,"-results.tsv",sep='')
plotfile <- paste( basename ,"-plot",sep='')

# Inference.
winlen <- lwin+win+rwin

genmatrix <- meangenmatrix( lwin=1, rwin=1, patlen=winlen, mutpats=mutpats, selpats=selpats, mutrates=mutrates*tlen, selcoef=selcoef, boundary="none" )
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

counts <- list(
            counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[1]]$initseq, simseqs[[1]]$finalseq, lwin=lwin )
        )
# want only patterns with leftmost possible position changed
nonoverlapping <- leftchanged(rownames(counts[[1]]),colnames(counts[[1]]),lwin=lwin,win=win)
nov.counts <- counts[[1]][nonoverlapping]
initcounts <- rowSums(counts[[1]])

nmuts <- length(mutpats)
nsel <- length(selpats)
# (quasi)-likelihood function using all counts -- binomial
likfun <- function (params) {
        # params are: mutrates, selcoef
        mutrates <- params[1:nmuts]
        selcoef <- params[nmuts+1:nsel]
        genmatrix@x <- update(genmatrix,mutrates=mutrates,selcoef=selcoef)
        # this is collapsed transition matrix
        subtransmatrix <- computetransmatrix( genmatrix, projmatrix )
        # return negative log-likelihood 
        (-1) * sum( counts[[1]] * log(subtransmatrix) )
}
# using only nonoverlapping counts, plus priors -- indep't poisson.
lud <- function (params) {
    # params are: mutrates*tlen, selcoef
    mutrates <- params[1:nmuts]
    selcoef <- params[nmuts+1:nsel]
    if (any(mutrates<0)) {
        return( -Inf )
    } else {
        genmatrix@x <- update(genmatrix,mutrates=mutrates,selcoef=selcoef)
        # this is collapsed transition matrix
        meancounts <- initcounts * computetransmatrix( genmatrix, projmatrix )
        # return (positive) log-posterior
        return( (-1)*sum(meancounts[nonoverlapping]) + sum( nov.counts * log(meancounts[nonoverlapping]) ) - sum(mmean*mutrates) - sum(selcoef^2)/svar )
    }
}

# point estimates
initpar <- c( 2 * runif( length(mutpats) ) * mean(mutrates) * tlen, 2 * runif( length(selpats) ) * mean(selcoef) ) # random init
truth <- c( mutrates * tlen, selcoef )  # truth
lbs <- c( rep(1e-6,nmuts), rep(-Inf,nsel) )
cheating.ans <- optim( par=truth, fn=likfun, method="L-BFGS-B", lower=lbs, control=list(trace=3) )
random.ans <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, control=list(trace=3) )

estimates <- data.frame( rbind(init=initpar, ans=random.ans$par, cheating=cheating.ans$par, truth=truth ) )
colnames(estimates) <- c( paste("muttime",seq_along(mutrates),sep=''), paste("selcoef",seq_along(selcoef),sep='') )
estimates$likfun <- apply( estimates, 1, likfun )
write.table( estimates, file=resultsfile, quote=FALSE, sep="\t" )

# bayesian
mrun <- metrop( lud, initial=random.ans$par, nbatch=nbatches, blen=blen, scale=1e-2 )

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

save( counts, genmatrix, subtransmatrix, lud, likfun, truth, cheating.ans, random.ans, nonoverlapping, nov.counts, mmean, svar, all.expected, cwin, subcounts, all.subexpected, mrun, file=datafile )

pdf(file=paste(plotfile,"-mcmc.pdf",sep=''),width=6, height=4, pointsize=10)
matplot( mrun$batch, type='l', lty=c(rep(1,nmuts),rep(2,nsel)), col=1:length(truth) )
abline(h=truth, col=adjustcolor(1:length(truth),.5), lty=c(rep(1,nmuts),rep(2,nsel)), lwd=2)
abline(h=estimates["ans",], col=1:length(truth), lty=c(rep(1,nmuts),rep(2,nsel)) )
legend("topright",lty=c(rep(1,nmuts),rep(2,nsel),1,1), col=c(1:length(truth),1,adjustcolor(1,.5)), lwd=c(rep(1,length(truth)),1,2),legend=c(colnames(estimates)[1:length(truth)],"point estimate","truth"))
dev.off()

pdf(file=paste(plotfile,"-1.pdf",sep=''),width=6, height=4, pointsize=10)
lord <- order( all.expected[["truth"]][,1] )
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



print(format(Sys.time(),"%Y-%m-%d-%H-%M"))
