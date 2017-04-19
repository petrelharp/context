#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Run mcmc longer.\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file containing simulation." ),
        make_option( c("-w","--shortwin"), type="integer", default=1, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--leftwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rightwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-n","--nbatches"), type="integer", default=1000, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-c","--stepscale"), type="numeric", default=1e-4, help="Scale of proposal steps for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-s","--restart"), action="store_true", default=FALSE, help="Start a whole new MCMC run?" ),
        make_option( c("-d","--boundary"), type="character", default="none", help="Boundary conditions for generator matrix. [default \"%default\"]"),
        make_option( c("-y","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases in computing generator matrix. [default \"%default\"]" ),
        make_option( c("-g","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix, or TRUE [default] to look for one. (otherwise, will compute)"),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)

if (interactive()) { nbatches <- 100; blen <- 10; restart <- FALSE; gmfile <- TRUE }

if (!'infile'%in%names(opt)) { stop("Run\n  cpg-mcmc.R -h\n for help.") }

library(contextual)
library(contextutils)
library(simcontext)

library(mcmc)

basedir <- gsub(".RData","",infile,fixed=TRUE)
load(infile)
if (!file.exists(basedir)) {
    dir.create(basedir)
}



basename <- paste(basedir,"/win-",leftwin,"-",shortwin,"-",rightwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
plotfile <- paste( basename ,"-plot",sep='')
mcmcdatafiles <- list.files(path=basedir,pattern="-mcmc.*RData",full.names=TRUE)
mcmcnum <- 1+max(c(0,as.numeric(gsub(".*-mcmc-","",gsub(".RData","",mcmcdatafiles)))),na.rm=TRUE)

if (logfile=="") { logfile <- paste(basename,"-mcmc-run-",mcmcnum,".Rout",sep='') }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    if (!interactive()) { sink(file=logcon, type="message") }
    sink(file=logcon, type="output", split=interactive()) 
}

load(datafile)  # has mrun and previous things
if (length(mcmcdatafiles)>0) { load(grep(paste("-mcmc-",mcmcnum-1,".RData",sep=''),mcmcdatafiles,fixed=TRUE,value=TRUE)) }

########

# Inference.
longwin <- leftwin+shortwin+rightwin
if (gmfile=="TRUE") { gmfile <- paste(paste("genmatrices/genmatrix",longwin,boundary,meanboundary,sep="-"),".RData",sep='') }
if (file.exists(gmfile)) {
    load(gmfile) 
} else {
    if (meanboundary>0) {
        genmatrix <- meangenmatrix( leftwin=1, rightwin=1, patlen=longwin, mutpats=mutpats, selpats=list(), mutrates=mutrates*tlen, selcoef=numeric(0), boundary=boundary )
    } else {
        genmatrix <- makegenmatrix( patlen=longwin, mutpats=mutpats, selpats=list(), mutrates=mutrates*tlen, selcoef=numeric(0), boundary=boundary )
    }
}
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin )
counts <- list( counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[1]]$initseq, simseqs[[1]]$finalseq, leftwin=leftwin ) )
# want only patterns with leftmost possible position changed
nonoverlapping <- leftchanged(rownames(counts[[1]]),colnames(counts[[1]]),leftwin=leftwin,shortwin=shortwin)
nov.counts <- counts[[1]][nonoverlapping]
initcounts <- rowSums(counts[[1]])
nmuts <- length(mutpats)
# using only nonoverlapping counts, plus priors -- indep't poisson.
lud <- function (mutrates) {
    # params are: mutrates*tlen, selcoef
    if (any(mutrates<0)) {
        return( -Inf )
    } else {
        genmatrix@x <- update(genmatrix,mutrates=mutrates,selcoef=numeric(0))
        # this is collapsed transition matrix
        meancounts <- initcounts * computetransmatrix( genmatrix, projmatrix )
        # return (positive) log-posterior
        return( 
            (-1)*sum(meancounts[nonoverlapping]) 
            + sum( nov.counts * log(meancounts[nonoverlapping]) ) 
            - sum(mmeans*mutrates) )
    }
}

if (restart) {
    mrun <- metrop( mrun, initial=random.ans$par, nbatch=nbatches, blen=blen, scale=stepscale )
} else {
    mrun <- metrop( mrun, nbatch=nbatches, blen=blen, scale=stepscale )
}


save( leftwin, shortwin, rightwin, lud, mrun, initcounts, nov.counts, nonoverlapping, file=paste(basename,"-mcmc-",mcmcnum,".RData",sep='') )

pdf(file=paste(plotfile,"-mcmc-",mcmcnum,".pdf",sep=''),width=6, height=4, pointsize=10)
matplot( mrun$batch, type='l', col=1:length(truth) )
abline(h=truth, col=adjustcolor(1:length(truth),.5), lwd=2)
abline(h=estimates["ans",], col=1:length(truth) )
legend("topright", col=c(1:length(truth),1,adjustcolor(1,.5)), lwd=c(rep(1,length(truth)),1,2),legend=c(colnames(estimates)[1:length(truth)],"point estimate","truth"))
dev.off()

mutlabels <- paste("mut:", unlist( sapply( sapply( mutpats, lapply, paste, collapse="->" ), paste, collapse=" | " ) ) )
pdf(file=paste(plotfile,"-mcmc-",mcmcnum,"-pairwise.pdf",sep=''),width=6,height=6,pointsize=10)
pairs( rbind( mrun$batch, truth ), col=c( rep(adjustcolor("black",.1),nrow(mrun$batch)), adjustcolor("red",1)), pch=20, cex=c( rep(.25,nrow(mrun$batch)), 2), labels=mutlabels, gap=.1 )
dev.off()
