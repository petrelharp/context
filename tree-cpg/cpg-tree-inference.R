#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Infer parameters from output of sim-tree-cpg.R .\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file containing simulation." ),
        make_option( c("-w","--win"), type="integer", default=1, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--lwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-n","--nbatches"), type="integer", default=20, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-s","--stepscale"), type="numeric", default=1e-4, help="Scale of proposal steps for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-m","--mmean"), type="double", default=1, help="Prior mean on single base mutation rates. [default \"%default\"]" ),
        make_option( c("-c","--cpgmean"), type="double", default=1, help="Prior variance on CpG rate. [default \"%default\"]" ),
        make_option( c("-p","--pmean"), type="double", default=1, help="Prior Dirichlet parameter for base frequencies. [default \"%default\"]" ),
        make_option( c("-v","--tvar"), type="double", default=.5, help="Prior variance for Gaussian prior on branch length. [default \"%default\"]" ),
        make_option( c("-d","--boundary"), type="character", default="none", help="Boundary conditions for generator matrix. [default \"%default\"]"),
        make_option( c("-y","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases in computing generator matrix. [default \"%default\"]" ),
        make_option( c("-g","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix, or TRUE [default] to look for one. (otherwise, will compute)"),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)
options(error=traceback)

winlen <- lwin+win+rwin

if (gmfile=="TRUE") { gmfile <- paste(paste("genmatrices/genmatrix",winlen,boundary,meanboundary,sep="-"),".RData",sep='') }

if (is.null(infile)) { cat("Run\n  cpg-tree-inference.R -h\n for help.") }

if (logfile!="" & !interactive()) { 
    logfile <- gsub(".RData",".Rout",infile,fixed=TRUE)
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message") 
    sink(file=logcon, type="output")   # send both to log file
}

scriptdir <- "../"
source(paste(scriptdir,"codon-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))

require(mcmc)

load(infile)
basedir <- gsub(".RData","",infile,fixed=TRUE)
if (!file.exists(basedir)) { dir.create(basedir) }
basename <- paste(basedir,"/win-",lwin,"-",win,"-",rwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
resultsfile <- paste( basename ,"-results.tsv",sep='')
plotfile <- paste( basename ,"-plot",sep='')

if (file.exists(gmfile)) {
    load(gmfile)
} else {
    if (meanboundary>0) {
        genmatrix <- meangenmatrix( lwin=1, rwin=1, patlen=winlen, mutpats=mutpats, selpats=list(), mutrates=mutrates*tlen[1], selcoef=numeric(0), boundary=boundary )
    } else {
        genmatrix <- makegenmatrix( patlen=winlen, mutpats=mutpats, selpats=list(), mutrates=mutrates*tlen[1], selcoef=numeric(0), boundary=boundary )
    }
}
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

counts <- list(
            "1.2"=counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[1]]$finalseq, simseqs[[2]]$finalseq, lwin=lwin ),
            "2.1"=counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[2]]$finalseq, simseqs[[1]]$finalseq, lwin=lwin )
        )
# want only patterns that overlap little
initcounts <- lapply( counts, rowSums )
nonoverlapping <- lapply( seq_along(counts), function (k) ( leftchanged(rownames(counts[[k]]),colnames(counts[[k]]),lwin=lwin,win=win) & (initcounts[[k]]>0) ) )
nov.counts <- lapply(seq_along(counts), function (k) counts[[k]][nonoverlapping[[k]]] )

# move from base frequencies (what we estimate) to pattern frequencies
nmuts <- length(mutpats)
nfreqs <- length(initfreqs)
npats <- nrow(genmatrix)
patcomp <- apply( do.call(rbind, strsplit(rownames(genmatrix),'') ), 2, match, bases )  # which base is at each position in each pattern
likfun <- function (params) {
    # params are: tlen[2]/tlen[1], tlen[1]*mutrates, basefreqs
    branchlens <- c(1,params[1])
    mutrates <- params[1+(1:nmuts)]
    initfreqs <- params[1+nmuts+(1:(nfreqs-1))]
    initfreqs <- c(1-sum(initfreqs),initfreqs)
    patfreqs <- initfreqs[patcomp]
    dim(patfreqs) <- dim(patcomp)
    patfreqs <- apply( patfreqs, 1, prod )
    # these are collapsed transition matrix
    updownbranch <- list(
            getupdowntrans( genmatrix, projmatrix, mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=branchlens ),
            getupdowntrans( genmatrix, projmatrix, mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=rev(branchlens) )
        )
    # return negative log-likelihood 
    (-1) * ( sum( counts[[1]] * log(updownbranch[[1]]) ) + sum( counts[[2]] * log(updownbranch[[2]]) ) )
}

# using only nonoverlapping counts, plus priors -- indep't poisson.
mmeans <- c( rep(mmean,nmuts-1), cpgmean )
pmeans <- rep( pmean, nfreqs )
lud <- function (params) {
    # params are: tlen[2]/tlen[1], tlen[1]*mutrates, basefreqs[-1]
    branchlens <- c(1,params[1])
    mutrates <- params[1+(1:nmuts)]
    initfreqs <- params[1+nmuts+(1:(nfreqs-1))]
    initfreqs <- c(1-sum(initfreqs),initfreqs)
    patfreqs <- initfreqs[patcomp]
    dim(patfreqs) <- dim(patcomp)
    patfreqs <- rowSums( patfreqs )
    # params are: mutrates*tlen
    if (any(mutrates<0) | any(initfreqs<0)) {
        return( -Inf )
    } else {
        # only do in one direction... ?
        updownbranch <- getupdowntrans( genmatrix, projmatrix, mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=branchlens )
        # return (positive) log-posterior
        return( 
                (-1)*sum(updownbranch[nonoverlapping[[1]]]) 
                + sum( nov.counts[[1]] * log(updownbranch[nonoverlapping[[1]]]) ) 
                - (branchlens[2]-1)/tvar
                - sum(mmeans*mutrates) 
                + sum( (pmeans-1)*log(initfreqs) )
            )
    }
}

# point estimates
initpar <- c( 2 * runif(1), 2 * runif( nmuts ) * mean(mutrates) * tlen[1], runif(length(initfreqs)-1)/length(initfreqs) ) # random init
truth <- c( tlen[2]/tlen[1], mutrates * tlen[1], initfreqs[-1] )  # truth
lbs <- c( 1e-6, rep(1e-6,nmuts), rep(1e-6,length(initfreqs)-1) )
ubs <- c( 20, rep(20,nmuts), rep(1,length(initfreqs)) )
cheating.ans <- optim( par=truth, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(trace=3,fnscale=likfun(truth)) )
random.ans <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(trace=3,fnscale=likfun(truth)) )

estimates <- data.frame( rbind(init=initpar, ans=random.ans$par, cheating=cheating.ans$par, truth=truth ) )
colnames(estimates) <- paste("muttime",seq_along(mutrates),sep='')
estimates$likfun <- apply( estimates, 1, likfun )
write.table( estimates, file=resultsfile, quote=FALSE, sep="\t" )

# bayesian
mrun <- metrop( lud, initial=random.ans$par, nbatch=nbatches, blen=blen, scale=stepscale )

# look at observed/expected counts
all.expected <- lapply( 1:nrow(estimates), function (k) {
            x <- unlist(estimates[k,])
            list( predictcounts( win, lwin, rwin, initcounts=rowSums(counts[[1]]), mutrates=x[1:nmuts], selcoef=numeric(0), genmatrix=genmatrix, projmatrix=projmatrix ) )
    } )
names(all.expected) <- rownames(estimates)

# look at observed/expected counts in smaller windows
cwin <- min(2,win); lrcwin <- min(1,lwin,rwin)
subcounts <- projectcounts( lwin=lwin, countwin=cwin, lcountwin=lrcwin, rcountwin=lrcwin, counts=counts[[1]] )
all.subexpected <- lapply( all.expected, function (x) { list( projectcounts( lwin=lwin, countwin=cwin, lcountwin=lrcwin, rcountwin=lrcwin, counts=x[[1]] ) ) } )

save( opt, counts, genmatrix, projmatrix, subtransmatrix, lud, likfun, truth, cheating.ans, random.ans, estimates, initpar, nonoverlapping, nov.counts, mmeans, all.expected, cwin, subcounts, all.subexpected, mrun, win, lwin, rwin, nmuts, file=datafile )

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
