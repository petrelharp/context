#!/usr/bin/Rscript --vanilla
require(optparse)
options(error=traceback)

usage <- "\
Run mcmc longer.\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file containing simulation." ),
        make_option( c("-n","--nbatches"), type="integer", default=1000, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-r","--restart"), action="store_true", default=FALSE, help="Start a whole new MCMC run?" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)

if (is.null(infile) | is.null(nbatches)) { cat("Run\n  ising-inference.R -h\n for help.") }

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

load(datafile)  # has mrun and everything else we might need in it

next.mrun <- metrop( mrun, initial=random.ans$par, nbatch=nbatches, blen=blen, scale=1e-2 )

pdf(file=paste(plotfile,"-mcmc-",mcmcnum,".pdf",sep=''),width=6, height=4, pointsize=10)
matplot( mrun$batch, type='l', lty=c(rep(1,nmuts),rep(2,nsel)), col=1:length(truth) )
abline(h=truth, col=adjustcolor(1:length(truth),.5), lty=c(rep(1,nmuts),rep(2,nsel)), lwd=2)
abline(h=estimates["ans",], col=1:length(truth), lty=c(rep(1,nmuts),rep(2,nsel)) )
legend("topright",lty=c(rep(1,nmuts),rep(2,nsel),1,1), col=c(1:length(truth),1,adjustcolor(1,.5)), lwd=c(rep(1,length(truth)),1,2),legend=c(colnames(estimates)[1:length(truth)],"point estimate","truth"))
dev.off()

save( next.mrun, file=paste(basename,"-mcmc-",mcmcnum,".RData",sep='') )
