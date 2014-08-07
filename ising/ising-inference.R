#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Infer parameters from output of sim-ising.R .\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file containing simulation." ),
        make_option( c("-w","--shortwin"), type="integer", default=3, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--leftwin"), type="integer", default=3, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rightwin"), type="integer", default=3, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-n","--nbatches"), type="integer", default=20, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-m","--mmean"), type="double", default=1, help="Prior mean on mutation rates. [default \"%default\"]" ),
        make_option( c("-v","--svar"), type="double", default=1, help="Prior variance on selection coefficents. [default \"%default\"]" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)
options(error=traceback)

if (interactive()) { shortwin <- leftwin <- rightwin <- 3; nbatches <- 10; blen <- 10; mmean <- 1; svar <- 1 }

if (is.null(infile)) { cat("Run\n  ising-inference.R -h\n for help.") }

if (logfile!="" & !interactive()) { 
    logfile <- gsub(".RData",".Rout",infile,fixed=TRUE)
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message") 
    sink(file=logcon, type="output")   # send both to log file
}

scriptdir <- "../"
source(paste(scriptdir,"context-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))

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

genmatrix <- meangenmatrix( leftwin=1, rightwin=1, patlen=longwin, mutpats=mutpats, selpats=selpats, mutrates=mutrates*tlen, selcoef=selcoef, boundary="none" )
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

counts <- list(
            counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[1]]$initseq, simseqs[[1]]$finalseq, leftwin=leftwin )
        )
# want only patterns that overlap little
initcounts <- rowSums(counts[[1]])
nonoverlapping <- ( leftchanged(rownames(counts[[1]]),colnames(counts[[1]]),leftwin=leftwin,shortwin=shortwin) & (initcounts>0) )
nov.counts <- counts[[1]][nonoverlapping]

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
lbs <- c( rep(1e-6,nmuts), rep(-20,nsel) )
ubs <- c( rep(20,nmuts), rep(20,nsel) )
cheating.ans <- optim( par=truth, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(trace=3,fnscale=likfun(truth)) )
random.ans <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(trace=3,fnscale=likfun(truth)) )

estimates <- data.frame( rbind(init=initpar, ans=random.ans$par, cheating=cheating.ans$par, truth=truth ) )
colnames(estimates) <- c( paste("muttime",seq_along(mutrates),sep=''), paste("selcoef",seq_along(selcoef),sep='') )
estimates$likfun <- apply( estimates, 1, likfun )
write.table( estimates, file=resultsfile, quote=FALSE, sep="\t" )

# bayesian
mrun <- metrop( lud, initial=random.ans$par, nbatch=nbatches, blen=blen, scale=1e-2 )

# look at observed/expected counts
all.expected <- lapply( 1:nrow(estimates), function (k) {
            x <- unlist(estimates[k,])
            list( predictcounts( shortwin, leftwin, rightwin, initcounts=rowSums(counts[[1]]), mutrates=x[1:nmuts], selcoef=x[nmuts+(1:nsel)], genmatrix=genmatrix, projmatrix=projmatrix ) )
    } )
names(all.expected) <- rownames(estimates)

# look at observed/expected counts in smaller windows
cwin <- 3; lrcwin <- 1
if ( shortwin >= cwin && longwin >= cwin + 2*lrcwin ) {
    subcounts <- projectcounts( leftwin=leftwin, countwin=cwin, lcountwin=lrcwin, rcountwin=lrcwin, counts=counts[[1]] )
    all.subexpected <- lapply( all.expected, function (x) { list( projectcounts( leftwin=leftwin, countwin=cwin, lcountwin=lrcwin, rcountwin=lrcwin, counts=x[[1]] ) ) } )
} else {
    subcounts <- all.subexpected <- NULL
}

save( counts, genmatrix, projmatrix, subtransmatrix, lud, likfun, truth, cheating.ans, random.ans, estimates, initpar, nonoverlapping, nov.counts, mmean, svar, all.expected, cwin, subcounts, all.subexpected, mrun, shortwin, leftwin, rightwin, nmuts, nsel, file=datafile )

pdf(file=paste(plotfile,"-mcmc.pdf",sep=''),width=6, height=4, pointsize=10)
layout(matrix(c(1,4,2,3),nrow=2))
par(mar=c(4,4,0,0)+.1)
for (k in 1:(ncol(mrun$batch)-1)) { 
    for (j in (k+1):ncol(mrun$batch)) {
        plot(mrun$batch[,j],mrun$batch[,k],xlab=colnames(estimates)[j],ylab=colnames(estimates)[k],pch=20,cex=.5,col=adjustcolor('black',.5), xlim=range(c(mrun$batch[,j],estimates[c("truth","ans"),j])), ylim=range(c(mrun$batch[,k],estimates[c("truth","ans"),k])))
        points( estimates["truth",j], estimates["truth",k], pch=20, col='green' )
        points( estimates["ans",j], estimates["ans",k], pch=20, col='red' )
    }
}
legend("bottomright",pch=20,col=c('green','red'),legend=c("truth","estimated"))
matplot( mrun$batch, type='l', lty=c(rep(1,nmuts),rep(2,nsel)), col=1:length(truth), xlab="mcmc gens" )
abline(h=truth, col=adjustcolor(1:length(truth),.5), lty=c(rep(1,nmuts),rep(2,nsel)), lwd=2)
abline(h=estimates["ans",], col=1:length(truth), lty=c(rep(1,nmuts),rep(2,nsel)) )
legend("topright",lty=c(rep(1,nmuts),rep(2,nsel),1,1), col=c(1:length(truth),1,adjustcolor(1,.5)), lwd=c(rep(1,length(truth)),1,2),legend=c(colnames(estimates)[1:length(truth)],"point estimate","truth"))
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
        lines(all.expected[["cheating"]][[j]][lord,k],col='blue',lty=j)
    }
    legend("topleft",legend=c("expected","estimated","cheating"),lty=1,col=c("red","green","blue"))
}
dev.off()

pdf(file=paste(plotfile,"-2.pdf",sep=''),width=6, height=4, pointsize=10)
layout(matrix(1:ncol(all.subexpected[["truth"]][[1]]),nrow=2))
cols <- rainbow(2+length(all.expected))[1:length(all.expected)]
for (k in 1:ncol(all.subexpected[["truth"]][[1]])) {
    lord <- order( all.subexpected[["truth"]][[1]][,k] )
    plot( subcounts[lord,k], xaxt='n', xlab='', main=colnames(subcounts)[k], log='y' )
    axis(1,at=1:nrow(subcounts),labels=rownames(subcounts)[lord],las=3)
    invisible( lapply(seq_along(all.subexpected),function(j) { lines(all.subexpected[[j]][[1]][lord,k],col=cols[j]) } ) )
    legend("topleft",legend=names(all.subexpected),lty=1,col=cols)
}
dev.off()

pdf(file=paste(plotfile,"-3.pdf",sep=''),width=6, height=4, pointsize=10)
layout(1)
plot( as.vector(all.expected[["truth"]][[1]]), as.vector(counts[[1]]), log='xy', xlab="true expected counts", ylab="counts" )
abline(0,1)
points(as.vector(all.expected[["truth"]][[1]]), as.vector(all.expected[["cheating"]][[1]]), col='red', pch=20 )
points(as.vector(all.expected[["truth"]][[1]]), as.vector(all.expected[["ans"]][[1]]), col='green', pch=20, cex=.5)
legend("topleft",pch=c(1,20,20),col=c('black','red','green'),legend=c('observed','cheating','estimated'))
dev.off()



print(format(Sys.time(),"%Y-%m-%d-%H-%M"))
