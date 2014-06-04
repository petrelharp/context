#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Infer mutation rates.
"

option_list <- list(
        make_option( c("-u","--indir"), type="character", default=NULL, help="Directory to look for input in, and write output files to." ),
        make_option( c("-i","--infile"), type="character", default=NULL, help="Table of count data."),
        make_option( c("-v","--revfile"), type="character", default=NULL, help="Table of count data in the reverse orientation"),
        make_option( c("-j","--jobid"), type="character", default=formatC(1e6*runif(1),width=6,format="d",flag="0"), help="Unique job id. [default random]"),
        make_option( c("-w","--win"), type="integer", default=1, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--lwin"), type="integer", default=1, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rwin"), type="integer", default=1, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-n","--nbatches"), type="integer", default=20, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-s","--stepscale"), type="numeric", default=1e-4, help="Scale of proposal steps for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-a","--optimscale"), type="numeric", default=1e-2, help="Scale of proposal steps for optimization algorithm. [default \"%default\"]" ),
        make_option( c("-m","--mmean"), type="double", default=1, help="Prior mean on single base mutation rates. [default \"%default\"]" ),
        make_option( c("-c","--cpgmean"), type="double", default=1, help="Prior variance on CpG rate. [default \"%default\"]" ),
        make_option( c("-p","--pprior"), type="double", default=1, help="Parameter for Dirichlet prior on base frequencies. [default \"%default\"]" ),
        make_option( c("-t","--tprior"), type="double", default=.5, help="Parameter for Beta prior on branch length. [default \"%default\"]" ),
        make_option( c("-d","--boundary"), type="character", default="none", help="Boundary conditions for generator matrix. [default \"%default\"]"),
        make_option( c("-y","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases in computing generator matrix. [default \"%default\"]" ),
        make_option( c("-g","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix, or TRUE [default] to look for one. (otherwise, will compute)"),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$infile) & is.null(opt$indir)) { stop("No input file.  Run\n  bcells-inference.R -h\n for help.\n") }
attach(opt)
options(error=traceback)

winlen <- lwin+win+rwin

if (gmfile=="TRUE") { 
    gmfile <- paste(paste("genmatrices/genmatrix",winlen,boundary,meanboundary,sep="-"),".RData",sep='') 
}

if (substr(indir,nchar(indir),nchar(indir)) %in% c("/","\\")) { indir <- substr(indir,1,nchar(indir)-1) }
if (is.null(opt$infile)) { infile <- paste(indir,"/", winlen,".",win,".counts",sep='') }
if (is.null(opt$revfile)) { revfile <- paste(dirname(infile),"/rev.",basename(infile),sep='') }
if (!file.exists(infile) | !file.exists(revfile)) { stop("Cannot read file ", infile) }

basedir <- paste(infile,"-results",sep='')
if (!file.exists(basedir)) { dir.create(basedir) }
basename <- paste(basedir,"/win-",win,"-",lwin,"-",rwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
resultsfile <- paste( basename ,"-results.tsv",sep='')
plotfile <- paste( basename ,"-plot",sep='')

if (logfile=="") {
    logfile <- paste(basename,".Rout",sep="")
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    if (!interactive()) { sink(file=logcon, type="message") }
    sink(file=logcon, type="output", split=interactive())   # send both to log file
}

scriptdir <- "../"
source(paste(scriptdir,"context-inference-fns.R",sep=''))
# source(paste(scriptdir,"sim-context-fns.R",sep=''))

require(mcmc)

if (file.exists(gmfile)) {
    load(gmfile)
} else {
    stop(paste("Cannot read",gmfile,"."))
}
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin )

# read in counts (produced with count-paired-tuples.py)
counts <- lapply( list(infile,revfile), function (ifile) {
        count.table <- read.table(ifile,header=TRUE,stringsAsFactors=FALSE)
        counts <- Matrix(0,nrow=nrow(genmatrix),ncol=ncol(projmatrix))
        rownames(counts) <- rownames(genmatrix)
        colnames(counts) <- colnames(projmatrix)
        stopifnot( all( count.table$reference %in% rownames(genmatrix) ) & all(count.table$derived %in% colnames(projmatrix)) ) 
        counts[cbind( match(count.table$reference,rownames(genmatrix)), match(count.table$derived,colnames(projmatrix)) )] <- count.table$count
        return(counts)
    } )
initcounts <- lapply( counts, rowSums )

# simple point estimates for starting positions
simple.counts <- lapply(counts, countmuts,mutpats=mutpats,lwin=lwin)
adhoc <- lapply( simple.counts, function (x) x[1,]/x[2,] )

# move from base frequencies (what we estimate) to pattern frequencies
nmuts <- length(mutpats)
nfreqs <- length(bases)
npats <- nrow(genmatrix)
patcomp <- apply( do.call(rbind, strsplit(rownames(genmatrix),'') ), 2, match, bases )  # which base is at each position in each pattern
likfun <- function (params) {
    # params are: tlen[1]/sum(tlen), sum(tlen)*mutrates, initfreqs
    branchlens <- c(params[1],1-params[1])
    mutrates <- params[1+(1:nmuts)]
    initfreqs <- params[1+nmuts+(1:nfreqs)]
    initfreqs <- initfreqs/sum(initfreqs)
    patfreqs <- initfreqs[patcomp]
    dim(patfreqs) <- dim(patcomp)
    patfreqs <- apply( patfreqs, 1, prod )
    # these are collapsed transition matrix
    updownbranch <- list(  # note "up" branch is from simpler summaries
            getupdowntrans( genmatrix, projmatrix, mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=rev(branchlens) ),
            getupdowntrans( genmatrix, projmatrix, mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=branchlens )
        )
    # if (any(sapply(updownbranch,function(x) any(!is.numeric(x))|any(x<0)))) { browser() }
    # return negative log-likelihood plus a penalty to keep initfreqs summing to (almost) 1
    return( 
                (-1) * ( sum( counts[[1]] * log(updownbranch[[1]]) ) + sum( counts[[2]] * log(updownbranch[[2]]) ) ) 
                + 100*(sum(initfreqs)-1)^2
            )
}

# only nonoverlapping counts, plus priors -- indep't trials. (multinomial)
mmeans <- c( rep(mmean,nmuts-1), cpgmean )
ppriors <- rep( pprior, nfreqs )
tpriors <- rep( tprior, 2 )
lud <- function (params) {
    # params are: tlen[1]/sum(tlen), sum(tlen)*mutrates, initfreqs[-length(initfreqs)]
    branchlens <- c(params[1],1-params[1])
    mutrates <- params[1+(1:nmuts)]
    initfreqs <- params[1+nmuts+(1:(nfreqs-1))]
    initfreqs <- c(initfreqs,1-sum(initfreqs))
    patfreqs <- initfreqs[patcomp]
    dim(patfreqs) <- dim(patcomp)
    patfreqs <- apply( patfreqs, 1, prod )
    if (any(mutrates<0) | any(initfreqs<0)) {
        return( -Inf )
    } else {
        # only do in one direction... ?
        updownbranch <- getupdowntrans( genmatrix, projmatrix, mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=rev(branchlens) )
        # return (positive) log-posterior
        return( 
                # (-1)*sum(updownbranch[nonoverlapping[[1]]]) +
                sum( counts[[1]] * log(updownbranch) ) 
                + sum( (tpriors-1)*log(branchlens) )
                - sum(mmeans*mutrates) 
                + sum( (ppriors-1)*log(initfreqs) )
            )
    }
}

# MLE point estimates
#   params are: tlen[1]/sum(tlen), sum(tlen)*mutrates, initfreqs
initpar <- c( runif(1), adhoc[[1]], rep(.25,4) ) # random init
names(initpar) <- c( "rel.branchlen", paste("tmut:", unlist( sapply( sapply( mutpats, lapply, paste, collapse="->" ), paste, collapse=" | " ) ) ), paste("init",bases,sep='.') )
lbs <- c( 1e-6, rep(0,nmuts), rep(1e-6,length(bases)) )
ubs <- c( 1, rep(20,nmuts), rep(1,length(bases)) )
# # do in parallel:
# cheating.parjob <- mcparallel( optim( par=truth, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(trace=3,fnscale=likfun(truth)) ) ),
# random.parjob <- mcparallel( optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(trace=3,fnscale=likfun(truth)) ) )
# cheating.ans <- mccollect( cheating.parjob, wait=TRUE )
# mle <- mccollect( random.parjob, wait=TRUE )
mle <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(trace=3,fnscale=abs(likfun(initpar)),parscale=rep(optimscale,length(initpar))) )
# renormalize the initial frequencies
mle.par <- mle$par; mle.par[1+nmuts+(1:nfreqs)] <- mle.par[1+nmuts+(1:nfreqs)] / sum( mle.par[1+nmuts+(1:nfreqs)] )

estimates <- data.frame( rbind(init=initpar, mle=mle.par) )
colnames(estimates) <- c( "branchlen", paste("tmut:", unlist( sapply( sapply( mutpats, lapply, paste, collapse="->" ), paste, collapse=" | " ) ) ), paste("init",bases,sep='.') )
estimates$likfun <- apply( estimates, 1, likfun )
write.table( estimates, file=resultsfile, quote=FALSE, sep="\t" )

# bayesian
#  note we deal with initial freqs not summing to 1 differently from in optim( ) -- need to remove last entry.
#  mrun.parjob <- mcparallel( metrop( lud, initial=mle.par[-length(mle.par)], nbatch=nbatches, blen=blen, scale=stepscale ) )
#  mrun <- mccollect(mcrun.parjob)
mrun <- metrop( lud, initial=mle.par[-length(mle.par)], nbatch=nbatches, blen=blen, scale=stepscale )

# look at observed/expected counts
all.expected <- lapply( 1:nrow(estimates), function (k) {
            x <- unlist(estimates[k,])
            branchlens <- c(x[1],1-x[1])
            mutrates <- x[1+(1:nmuts)]
            initfreqs <- x[1+nmuts+(1:nfreqs)]
            initfreqs <- initfreqs/sum(initfreqs)
            names(initfreqs) <- bases
            list( 
                    predicttreecounts( win, lwin, rwin, initcounts=rowSums(counts[[1]]), mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), genmatrix=genmatrix, projmatrix=projmatrix, initfreqs=initfreqs, tlens=rev(branchlens) ),
                    predicttreecounts( win, lwin, rwin, initcounts=rowSums(counts[[2]]), mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), genmatrix=genmatrix, projmatrix=projmatrix, initfreqs=initfreqs, tlens=branchlens )
                )
    } )
names(all.expected) <- rownames(estimates)

# look at observed/expected counts in smaller windows
cwin <- min(2,win); lrcwin <- min(1,lwin,rwin)
subcounts <- lapply( counts, function (x) 
        projectcounts( lwin=lwin, countwin=cwin, lcountwin=lrcwin, rcountwin=lrcwin, counts=x ) )
all.subexpected <- lapply( all.expected, lapply, function (x)
        projectcounts( lwin=lwin, countwin=cwin, lcountwin=lrcwin, rcountwin=lrcwin, counts=x ) )

save( opt, counts, genmatrix, projmatrix, lud, likfun, mle, estimates, initpar, mmeans, ppriors, tpriors, all.expected, cwin, subcounts, all.subexpected, mrun, win, lwin, rwin, nmuts, nfreqs, npats, patcomp, file=datafile )

# plot (long) counts
pdf(file=paste(plotfile,"-longcounts.pdf",sep=''),width=10, height=8, pointsize=10)
layout(matrix(1:sum(sapply(counts,ncol)),nrow=2))
for (j in seq_along(counts)) {
    for (k in 1:ncol(counts[[j]])) {
        lord <- order( all.expected[["mle"]][[j]][,k] )
        plot( counts[[j]][lord,k], type='n', xaxt='n', xlab='', ylim=range(c(unlist(all.expected[["mle"]]),unlist(as.matrix(counts[[j]])))), 
            ylab='counts', main=colnames(counts[[j]])[k] )
        axis(1,at=1:nrow(counts[[j]]),labels=rownames(counts[[j]])[lord],las=3)
        points( counts[[j]][lord,k], pch=j )
        lines(all.expected[["mle"]][[j]][lord,k],col='green', lty=j, lwd=2)
        lines(all.expected[["init"]][[j]][lord,k],col='red', lty=j, lwd=2)
        if (k==1) legend("topleft",legend=c("estimated","initial"),lty=1,col=c("green","red"))
    }
}
dev.off()

# plot (shorter) counts 
pdf(file=paste(plotfile,"-shortcounts.pdf",sep=''),width=10, height=8, pointsize=10)
layout(matrix(1:sum(sapply(subcounts,ncol)),nrow=2))
cols <- rainbow(2+length(all.expected))[1:length(all.expected)]
for (j in seq_along(subcounts)) {
    for (k in 1:ncol(subcounts[[j]])) {
        lord <- order( all.subexpected[["mle"]][[j]][,k] )
        plot( subcounts[[j]][lord,k], xaxt='n', xlab='', main=colnames(subcounts[[j]])[k], log='y' )
        axis(1,at=1:nrow(subcounts[[j]]),labels=rownames(subcounts[[j]])[lord],las=3)
        invisible( lapply(seq_along(all.subexpected),function(y) { lines(all.subexpected[[y]][[j]][lord,k],col=cols[y]) } ) )
        legend("topleft",legend=names(all.subexpected),lty=1,col=cols)
    }
}
dev.off()

# residuals of (shorter) counts 
pdf(file=paste(plotfile,"-shortresids.pdf",sep=''),width=10, height=8, pointsize=10)
layout(matrix(seq_along(subcounts)))
cols <- rainbow(2+length(all.expected))[1:length(all.expected)]
all.subresids <- lapply( all.subexpected, function (x) mapply(function(u,v) (u-v)/sqrt(v),x,subcounts) )
for (j in seq_along(counts)) {
    z <- sapply( lapply( all.subresids[c("mle","init")], "[[", j ), as.vector )
    rownames(z) <- paste( rownames(subcounts[[j]])[row(subcounts[[j]])], colnames(subcounts[[j]])[col(subcounts[[j]])], sep="->" )
    plot( 0, type='n', xlab="", ylab="normalized residuals", xlim=c(0,ncol(z)+1), ylim=range(z), xaxt='n' )
    text( jitter(col(z),fac=2), z, labels=rownames(z), col=col(z) )
    axis(1, at=1:ncol(z), labels=colnames(z) )
}
dev.off()


print(format(Sys.time(),"%Y-%m-%d-%H-%M"))
