#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Run mcmc longer.\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file containing simulation." ),
        make_option( c("-w","--shortwin"), type="integer", default=3, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--leftwin"), type="integer", default=3, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rightwin"), type="integer", default=3, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-n","--nbatches"), type="integer", default=1000, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-c","--stepscale"), type="numeric", default=3e-3, help="Scale of proposal steps for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-f","--shift"), type="integer", default=NULL, help="Number of bases to skip between window starting positions. [default equal to shortwin]" ),
        make_option( c("-s","--restart"), action="store_true", default=FALSE, help="Start a whole new MCMC run?" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$shift)) { opt$shift <- opt$shortwin }
attach(opt)

if (interactive()) { shortwin <- leftwin <- rightwin <- 3; nbatches <- 100; blen <- 10; restart <- FALSE }

if (!'infile'%in%names(opt)) { stop("Run\n  tasep-inference.R -h\n for help.") }

options(error=traceback)

basedir <- gsub(".RData","",infile,fixed=TRUE)
load(infile)
if (!file.exists(basedir)) {
    dir.create(basedir)
}

scriptdir <- "../"
source(paste(scriptdir,"context-inference-fns.R",sep=''),chdir=TRUE)
source(paste(scriptdir,"sim-context-fns.R",sep=''),chdir=TRUE)
require(mcmc)

# set-up
bases <- c("X","O")

basename <- paste(basedir,"/win-",leftwin,"-",shortwin,"-",rightwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
plotfile <- paste( basename ,"-plot",sep='')
mcmcdatafiles <- list.files(path=basedir,pattern="-mcmc.*RData",full.names=TRUE)
mcmcnum <- 1+max(c(0,as.numeric(gsub(".*-mcmc-","",gsub(".RData","",mcmcdatafiles)))),na.rm=TRUE)

if (logfile=="" & !interactive()) { logfile <- paste(basename,"-mcmc-run-",mcmcnum,".Rout",sep='') }
if (logfile != '' & !is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message") 
    sink(file=logcon, type="output") 
}

load(datafile)  # has mrun and previous things
if (length(mcmcdatafiles)) { load(grep(paste("-mcmc-",mcmcnum-1,".RData",sep=''),mcmcdatafiles,fixed=TRUE,value=TRUE)) }

########

# Inference.
longwin <- leftwin+shortwin+rightwin
genmatrix <- meangenmatrix( leftwin=1, rightwin=1, patlen=longwin, mutpats=mutpats, selpats=selpats, mutrates=mutrates*tlen, selcoef=selcoef, boundary="none" )
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin )
counts <- list( counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[1]]$initseq, simseqs[[1]]$finalseq, leftwin=leftwin, shift=shortwin ) )
# look only every (shortwin)-th window
initcounts <- lapply( counts[[1]], rowSums )
nmuts <- length(mutpats)
nsel <- length(selpats)
# using only nonoverlapping counts, plus priors -- indep't poisson.
lud <- function (mutrates,k=1,...) {
    # params are: mutrates*tlen, selcoef
    if (any(mutrates<0)) {
        return( -Inf )
    } else {
        genmatrix@x <- update(genmatrix,mutrates=mutrates,selcoef=numeric(0))
        # this is collapsed transition matrix
        meancounts <- initcounts[[k]] * computetransmatrix( genmatrix, projmatrix )
        # return (positive) log-posterior
        return( (-1)*sum(meancounts) + sum( counts[[1]][[k]] * log(meancounts) ) - sum(mmean*mutrates) )
    }
}

# choose randomly which window position to use?
shiftpos <- sample( 1:length(counts[[1]]), 1 )
if (restart) {
    initval <- runif(nmuts) * 2 
    mrun <- metrop( lud, initial=initval, nbatch=nbatches, blen=blen, scale=stepscale, k=shiftpos )
} else {
    mrun <- metrop( mrun, nbatch=nbatches, blen=blen, scale=stepscale, k=shiftpos )
}


save( leftwin, shortwin, rightwin, lud, mrun, file=paste(basename,"-mcmc-",mcmcnum,".RData",sep='') )

pdf(file=paste(plotfile,"-mcmc-",mcmcnum,".pdf",sep=''),width=6, height=4, pointsize=10)
matplot( mrun$batch, type='l', lty=c(rep(1,nmuts),rep(2,nsel)), col=1:length(truth) )
abline(h=truth, col=adjustcolor(1:length(truth),.5), lty=c(rep(1,nmuts),rep(2,nsel)), lwd=2)
abline(h=estimates["ans",], col=1:length(truth), lty=c(rep(1,nmuts),rep(2,nsel)) )
legend("topright",lty=c(rep(1,nmuts),rep(2,nsel),1,1), col=c(1:length(truth),1,adjustcolor(1,.5)), lwd=c(rep(1,length(truth)),1,2),legend=c(colnames(estimates)[1:length(truth)],"point estimate","truth"))
dev.off()

