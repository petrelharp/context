#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Run mcmc longer.\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file containing simulation." ),
        make_option( c("-w","--win"), type="integer", default=2, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--lwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-n","--nbatches"), type="integer", default=1000, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-s","--restart"), action="store_true", default=FALSE, help="Start a whole new MCMC run?" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)

options(error=traceback)

if (is.null(infile) | is.null(nbatches)) { cat("Run\n  ising-inference.R -h\n for help.") }

if (logfile=="" & !interactive()) { logfile <- gsub(".RData",".Rout",infile,fixed=TRUE) }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message") 
}

scriptdir <- "../"
source(paste(scriptdir,"codon-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))

require(mcmc)

load(infile)
basedir <- gsub(".RData","",infile,fixed=TRUE)
if (!file.exists(basedir)) {
    dir.create(basedir)
}

# set-up
bases <- c("X","O")

basename <- paste(basedir,"/win-",lwin,"-",win,"-",rwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
plotfile <- paste( basename ,"-plot",sep='')
mcmcplotfiles <- list.files(path=basedir,pattern="-mcmc*.pdf")
mcmcdatafiles <- list.files(path=basedir,pattern="-mcmc*.RData")
mcmcnum <- 1+max(c(0,(-1)*as.numeric(gsub(".*-mcmc","",gsub(".RData","",mcmcdatafiles)))),na.rm=TRUE)

load(datafile)  # has mrun and previous things

########

# Inference.
winlen <- lwin+win+rwin
genmatrix <- meangenmatrix( lwin=1, rwin=1, patlen=winlen, mutpats=mutpats, selpats=selpats, mutrates=mutrates*tlen, selcoef=selcoef, boundary="none" )
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin )
counts <- list( counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[1]]$initseq, simseqs[[1]]$finalseq, lwin=lwin ) )
# want only patterns with leftmost possible position changed
nonoverlapping <- leftchanged(rownames(counts[[1]]),colnames(counts[[1]]),lwin=lwin,win=win)
nov.counts <- counts[[1]][nonoverlapping]
initcounts <- rowSums(counts[[1]])
nmuts <- length(mutpats)
nsel <- length(selpats)
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

if (restart) {
    next.mrun <- metrop( mrun, initial=random.ans$par, nbatch=nbatches, blen=blen, scale=1e-2 )
} else {
    next.mrun <- metrop( mrun, nbatch=nbatches, blen=blen, scale=1e-2 )
}


pdf(file=paste(plotfile,"-mcmc-",mcmcnum,".pdf",sep=''),width=6, height=4, pointsize=10)
matplot( mrun$batch, type='l', lty=c(rep(1,nmuts),rep(2,nsel)), col=1:length(truth) )
abline(h=truth, col=adjustcolor(1:length(truth),.5), lty=c(rep(1,nmuts),rep(2,nsel)), lwd=2)
abline(h=estimates["ans",], col=1:length(truth), lty=c(rep(1,nmuts),rep(2,nsel)) )
legend("topright",lty=c(rep(1,nmuts),rep(2,nsel),1,1), col=c(1:length(truth),1,adjustcolor(1,.5)), lwd=c(rep(1,length(truth)),1,2),legend=c(colnames(estimates)[1:length(truth)],"point estimate","truth"))
dev.off()

save( next.mrun, file=paste(basename,"-mcmc-",mcmcnum,".RData",sep='') )
