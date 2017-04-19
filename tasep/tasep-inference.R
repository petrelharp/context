#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Infer parameters from output of sim-tasep.R .\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file containing simulation." ),
        make_option( c("-w","--shortwin"), type="integer", default=3, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--leftwin"), type="integer", default=3, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rightwin"), type="integer", default=3, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-n","--nbatches"), type="integer", default=20, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-m","--mmean"), type="double", default=1, help="Prior mean on jump rate. [default \"%default\"]" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)
options(error=traceback)

if (interactive()) { shortwin <- leftwin <- rightwin <- 3; nbatches <- 10; blen <- 10; mmean <- 1 }

if (is.null(infile) | is.null(nbatches)) { cat("Run\n  tasep-inference.R -h\n for help.") }

if (logfile!="" & !interactive()) { 
    logfile <- gsub(".RData",".Rout",infile,fixed=TRUE)
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message") 
    sink(file=logcon, type="output")   # send both to log file
}

library(contextual)
library(contextutils)
library(simcontext)

require(mcmc)

load(infile)
basedir <- gsub(".RData","",infile,fixed=TRUE)
if (!file.exists(basedir)) { dir.create(basedir) }
basename <- paste(basedir,"/win-",leftwin,"-",shortwin,"-",rightwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
resultsfile <- paste( basename ,"-results.tsv",sep='')
plotfile <- paste( basename ,"-plot",sep='')

# set-up
bases <- c("X","O")

# Inference.
longwin <- leftwin+shortwin+rightwin

genmatrix <- meangenmatrix( leftwin=1, rightwin=1, patlen=longwin, mutpats=mutpats, selpats=selpats, mutrates=mutrates*tlen, selcoef=numeric(0), boundary="none" )
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

counts <- list(
            counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[1]]$initseq, simseqs[[1]]$finalseq, leftwin=leftwin )
        )
# want only patterns that overlap little
initcounts <- rowSums(counts[[1]])
nonoverlapping <- ( leftchanged(rownames(counts[[1]]),colnames(counts[[1]]),leftwin=leftwin,shortwin=shortwin) & (initcounts>0) )
nov.counts <- counts[[1]][nonoverlapping]

# (quasi)-likelihood function using all counts -- binomial
likfun <- function (mutrates) {
    # params are: mutrates*tlen
    genmatrix@x <- update(genmatrix,mutrates=mutrates,selcoef=numeric(0))
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix )
    # return negative log-likelihood 
    (-1) * sum( counts[[1]] * log(subtransmatrix) )
}
# using only nonoverlapping counts, plus priors -- indep't poisson.
lud <- function (mutrates) {
    # params are: mutrates*tlen
    if (any(mutrates<0)) {
        return( -Inf )
    } else {
        genmatrix@x <- update(genmatrix,mutrates=mutrates,selcoef=numeric(0))
        # this is collapsed transition matrix
        meancounts <- initcounts * computetransmatrix( genmatrix, projmatrix )
        # return (positive) log-posterior
        return( (-1)*sum(meancounts[nonoverlapping]) + sum( nov.counts * log(meancounts[nonoverlapping]) ) - sum(mmean*mutrates) )
    }
}


truth <- c( mutrates * tlen )  # truth
ans <- optimize(likfun,lower=1e-8,upper=10)
if (interactive()) {
    parrange <- sort(unique(c(seq(0,5,length.out=100),seq(max(0,truth-.1),truth+.1,length.out=40))))
    likrange <- sapply(parrange,likfun)
    plot(parrange,likrange)
}

estimates <- data.frame( rbind(ans=ans$minimum,truth=truth ) )
colnames(estimates) <- "muttime"
estimates$likfun <- apply( estimates, 1, likfun )
write.table( estimates, file=resultsfile, quote=FALSE, sep="\t" )

# bayesian
mrun <- metrop( lud, initial=estimates["ans","muttime"], nbatch=nbatches, blen=blen, scale=1e-2 )

# look at observed/expected counts
all.expected <- lapply( 1:nrow(estimates), function (k) {
            x <- unlist(estimates[k,])
            list( predictcounts( shortwin, leftwin, rightwin, initcounts=rowSums(counts[[1]]), mutrates=x[1], selcoef=numeric(0), genmatrix=genmatrix, projmatrix=projmatrix ) )
    } )
names(all.expected) <- rownames(estimates)

# look at observed/expected counts in smaller windows
cwin <- 2
subcounts <- projectcounts( counts[[1]], new.shortwin=cwin, new.leftwin=0, new.longwin=cwin )
all.subexpected <- lapply( all.expected, function (x) { list( projectcounts( x[[1]], new.shortwin=cwin, new.leftwin=0, new.longwin=cwin ) ) } )

save( counts, genmatrix, projmatrix, subtransmatrix, lud, likfun, truth, ans, estimates, nonoverlapping, nov.counts, mmean, all.expected, cwin, subcounts, all.subexpected, mrun, shortwin, leftwin, rightwin, file=datafile )

pdf(file=paste(plotfile,"-mcmc.pdf",sep=''),width=6, height=4, pointsize=10)
plot( mrun$batch, type='l', xlab="mcmc gens", ylab="tlen*jumprate" )
abline(h=truth, lwd=2)
abline(h=estimates["ans",], lty=2 )
legend("topright",lty=1:2, legend=c("truth","point estimate"))
dev.off()

pdf(file=paste(plotfile,"-1.pdf",sep=''),width=6, height=4, pointsize=10)
layout(matrix(1:4,nrow=2))
for (k in 1:ncol(counts[[1]])) {
    lord <- order( all.expected[["truth"]][[1]][,k] )
    plot( counts[[1]][lord,k], type='n', xaxt='n', xlab='', ylim=range(c(unlist(all.expected[["truth"]]),unlist(lapply(counts,as.matrix)),unlist(all.expected[["ans"]]))), ylab='counts' )
    axis(1,at=1:nrow(counts[[1]]),labels=rownames(counts[[1]])[lord],las=3)
    for (j in 1:length(counts)) {
        points( counts[[j]][lord,k], pch=j )
        lines(all.expected[["truth"]][[j]][lord,k],col='red', lty=j)
        lines(all.expected[["ans"]][[j]][lord,k],col='green', lty=j, lwd=2)
    }
    legend("topleft",legend=c("expected","estimated"),lty=1,col=c("red","green"))
}
dev.off()

pdf(file=paste(plotfile,"-2.pdf",sep=''),width=6, height=4, pointsize=10)
layout(matrix(1:4,nrow=2))
cols <- rainbow(2+length(all.expected))[1:length(all.expected)]
for (k in 1:4) {
    lord <- order( all.subexpected[["truth"]][[1]][,k] )
    plot( subcounts[lord,k], xaxt='n', xlab='', main=colnames(subcounts)[k] )
    axis(1,at=1:nrow(subcounts),labels=rownames(subcounts)[lord],las=3)
    invisible( lapply(seq_along(all.subexpected),function(j) { lines(all.subexpected[[j]][[1]][lord,k],col=cols[j]) } ) )
    legend("topleft",legend=names(all.subexpected),lty=1,col=cols)
}
dev.off()

pdf(file=paste(plotfile,"-3.pdf",sep=''),width=6, height=4, pointsize=10)
layout(1)
plot( as.vector(all.expected[["truth"]][[1]]), as.vector(counts[[1]]), log='xy', xlab="true expected counts", ylab="counts" )
abline(0,1)
points(as.vector(all.expected[["truth"]][[1]]), as.vector(all.expected[["ans"]][[1]]), col='green', pch=20, cex=.5)
legend("topleft",pch=c(1,20),col=c('black','green'),legend=c('observed','estimated'))
dev.off()



print(format(Sys.time(),"%Y-%m-%d-%H-%M"))
