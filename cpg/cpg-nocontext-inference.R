#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Infer parameters from output of sim-cpg.R ,\
in the case where there is NO CONTEXT-DEPENDENCE.
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file containing simulation." ),
        make_option( c("-n","--nbatches"), type="integer", default=20, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-m","--mmean"), type="double", default=1, help="Prior mean on single base mutation rates. [default \"%default\"]" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)
options(error=traceback)

leftwin <- rightwin <- 0; shortwin <- 1
longwin <- leftwin + rightwin + shortwin

if (is.null(infile)) { cat("Run\n  cpg-inference.R -h\n for help.") }

if (logfile=="" & !interactive()) { 
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
basename <- paste(basedir,"/nocontext",sep='')
datafile <- paste( basename ,"-results.RData",sep='')
resultsfile <- paste( basename ,"-results.tsv",sep='')
plotfile <- paste( basename ,"-plot",sep='')

# single-base subs only
singlebase <- ( sapply( lapply( lapply( mutpats, lapply, nchar ), unlist ), max ) == 1 )
orig.mutpats <- mutpats
orig.mutrates <- mutrates
c.to.t <- ( sapply( lapply(mutpats, "[[", 1), "[", 1 ) == "C" ) & ( sapply( lapply(mutpats, "[[", 1), "[", 2 ) == "T" )
g.to.a <- ( sapply( lapply(mutpats, "[[", 1), "[", 1 ) == "G" ) & ( sapply( lapply(mutpats, "[[", 1), "[", 2 ) == "A" )
mutrates[ c.to.t ] <- mutrates[ c.to.t ] + mutrates[ !singlebase ] * initfreqs[match("C",bases)]
mutrates[ g.to.a ] <- mutrates[ g.to.a ] + mutrates[ !singlebase ] * initfreqs[match("G",bases)]
mutpats <- mutpats[ singlebase ]
mutrates <- mutrates[ singlebase ]

genmatrix <- makegenmatrix( patlen=1, mutpats=mutpats, selpats=list(), mutrates=mutrates*tlen, selcoef=numeric(0), boundary='none' )
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

counts <- list(
            counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[1]]$initseq, simseqs[[1]]$finalseq, leftwin=leftwin )
        )
# want only patterns that overlap little
initcounts <- rowSums(counts[[1]])
nonoverlapping <- ( leftchanged(rownames(counts[[1]]),colnames(counts[[1]]),leftwin=leftwin,shortwin=shortwin) & (initcounts>0) )
nov.counts <- counts[[1]][nonoverlapping]

nmuts <- length(mutrates)
# (quasi)-likelihood function using all counts -- binomial
likfun <- function (mutrates) {
    # params are: mutrates*tlen
    genmatrix@x <- update(genmatrix,mutrates=mutrates,selcoef=numeric(0))
    # this is collapsed transition matrix
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix )
    # return negative log-likelihood 
    (-1) * sum( counts[[1]] * log(subtransmatrix) )
}
# using only nonoverlapping counts, plus priors -- indep't poisson.
mmeans <- c( rep(mmean,nmuts) )
lud <- function (mutrates) {
    # params are: mutrates*tlen
    if (any(mutrates<0)) {
        return( -Inf )
    } else {
        genmatrix@x <- update(genmatrix,mutrates=mutrates,selcoef=numeric(0))
        # this is collapsed transition matrix
        meancounts <- initcounts * computetransmatrix( genmatrix, projmatrix )
        # return (positive) log-posterior
        return( (-1)*sum(meancounts[nonoverlapping]) + sum( nov.counts * log(meancounts[nonoverlapping]) ) - sum(mmeans*mutrates) )
    }
}

# point estimates
initpar <- c( 2 * runif( nmuts ) * mean(mutrates) * tlen ) # random init
truth <- c( mutrates * tlen )  # truth
lbs <- c( rep(1e-6,nmuts) )
ubs <- c( rep(20,nmuts) )
cheating.ans <- optim( par=truth, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(trace=3,fnscale=abs(likfun(truth))) )
random.ans <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(trace=3,fnscale=abs(likfun(truth))) )

estimates <- data.frame( rbind(init=initpar, ans=random.ans$par, cheating=cheating.ans$par, truth=truth ) )
colnames(estimates) <- paste("muttime",seq_along(mutrates),sep='')
estimates$likfun <- apply( estimates, 1, likfun )
write.table( estimates, file=resultsfile, quote=FALSE, sep="\t" )

# bayesian
mrun <- metrop( lud, initial=random.ans$par, nbatch=nbatches, blen=blen, scale=1e-3 )

# look at observed/expected counts
all.expected <- lapply( 1:nrow(estimates), function (k) {
            x <- unlist(estimates[k,])
            list( predictcounts( shortwin, leftwin, rightwin, initcounts=rowSums(counts[[1]]), mutrates=x[1:nmuts], selcoef=numeric(0), genmatrix=genmatrix, projmatrix=projmatrix ) )
    } )
names(all.expected) <- rownames(estimates)

# look at observed/expected counts in smaller windows
cwin <- min(2,shortwin); lrcwin <- min(1,leftwin,rightwin)
subcounts <- projectcounts( leftwin=leftwin, countwin=cwin, lcountwin=lrcwin, rcountwin=lrcwin, counts=counts[[1]] )
all.subexpected <- lapply( all.expected, function (x) { list( projectcounts( leftwin=leftwin, countwin=cwin, lcountwin=lrcwin, rcountwin=lrcwin, counts=x[[1]] ) ) } )

save( opt, counts, genmatrix, projmatrix, subtransmatrix, lud, likfun, truth, cheating.ans, random.ans, estimates, initpar, nonoverlapping, nov.counts, mmeans, all.expected, cwin, subcounts, all.subexpected, mrun, shortwin, leftwin, rightwin, nmuts, file=datafile )

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
matplot( mrun$batch, type='l', col=1:length(truth), xlab="mcmc gens" )
abline(h=truth, col=adjustcolor(1:length(truth),.5), lwd=2)
abline(h=estimates["ans",], col=1:length(truth) )
legend("topright",col=c(1:length(truth),1,adjustcolor(1,.5)), lwd=c(rep(1,length(truth)),1,2),legend=c(colnames(estimates)[1:length(truth)],"point estimate","truth"))
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
