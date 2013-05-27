#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Run mcmc longer.\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file containing simulation." ),
        make_option( c("-w","--win"), type="integer", default=3, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--lwin"), type="integer", default=3, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rwin"), type="integer", default=3, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-n","--nbatches"), type="integer", default=1000, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-c","--stepscale"), type="numeric", default=3e-3, help="Scale of proposal steps for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-s","--restart"), action="store_true", default=FALSE, help="Start a whole new MCMC run?" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)

if (interactive()) { win <- lwin <- rwin <- 3; nbatches <- 100; blen <- 10; restart <- FALSE }

options(error=traceback)

if (is.null(infile) | is.null(nbatches)) { cat("Run\n  tasep-inference.R -h\n for help.") }

basedir <- gsub(".RData","",infile,fixed=TRUE)
load(infile)
if (!file.exists(basedir)) {
    dir.create(basedir)
}

scriptdir <- "../"
source(paste(scriptdir,"codon-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))
require(mcmc)

# set-up
bases <- c("X","O")

basename <- paste(basedir,"/win-",lwin,"-",win,"-",rwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
plotfile <- paste( basename ,"-plot",sep='')
mcmcdatafiles <- list.files(path=basedir,pattern="-mcmc.*RData",full.names=TRUE)
mcmcnum <- 1+max(c(0,as.numeric(gsub(".*-mcmc-","",gsub(".RData","",mcmcdatafiles)))),na.rm=TRUE)

if (logfile=="" & !interactive()) { logfile <- paste(basedir,"/mcmc-run-",mcmcnum,".Rout",sep='') }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message") 
    sink(file=logcon, type="output") 
}

load(datafile)  # has mrun and previous things
if (length(mcmcdatafiles)) { load(grep(paste("-mcmc-",mcmcnum-1,".RData",sep=''),mcmcdatafiles,fixed=TRUE,value=TRUE)) }

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
lud <- function (mutrates) {
    # params are: mutrates*tlen, selcoef
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

if (restart) {
    mrun <- metrop( mrun, initial=random.ans$par, nbatch=nbatches, blen=blen, scale=stepscale )
} else {
    mrun <- metrop( mrun, nbatch=nbatches, blen=blen, scale=stepscale )
}


save( lwin, win, rwin, lud, mrun, file=paste(basename,"-mcmc-",mcmcnum,".RData",sep='') )

pdf(file=paste(plotfile,"-mcmc-",mcmcnum,".pdf",sep=''),width=6, height=4, pointsize=10)
matplot( mrun$batch, type='l', lty=c(rep(1,nmuts),rep(2,nsel)), col=1:length(truth) )
abline(h=truth, col=adjustcolor(1:length(truth),.5), lty=c(rep(1,nmuts),rep(2,nsel)), lwd=2)
abline(h=estimates["ans",], col=1:length(truth), lty=c(rep(1,nmuts),rep(2,nsel)) )
legend("topright",lty=c(rep(1,nmuts),rep(2,nsel),1,1), col=c(1:length(truth),1,adjustcolor(1,.5)), lwd=c(rep(1,length(truth)),1,2),legend=c(colnames(estimates)[1:length(truth)],"point estimate","truth"))
dev.off()

