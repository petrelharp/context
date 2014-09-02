#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Run mcmc longer.\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file containing prior MCMC run." ),
        make_option( c("-n","--nbatches"), type="integer", default=1000, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-c","--stepscale"), type="character", default="1e-4", help="Scale of proposal steps for mutation*time parameters for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-e","--initscale"), type="character", default="3e-3", help="Scale of proposal steps for initial frequencies for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-f","--treescale"), type="character", default="1e-4", help="Scale of proposal steps for relative tree branch lengths for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-z","--id"), type="integer", default=sprintf("%04.0f",10000*runif(1)), help="Identifier for this MCMC run."), 
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.character(opt$stepscale)) { opt$stepscale <- eval(parse(text=opt$stepscale)) }
if (is.character(opt$initscale)) { opt$initscale <- eval(parse(text=opt$initscale)) }
if (is.character(opt$treescale)) { opt$treescale <- eval(parse(text=opt$treescale)) }
attach(opt)

if (interactive()) { nbatches <- 10; blen <- 1 }

if (!'infile'%in%names(opt)) { stop("Run\n  cpg-inference.R -h\n for help.") }

# options(error=traceback)

scriptdir <- "../"
source(paste(scriptdir,"context-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))

basename <- paste( gsub("[0-9]*.RData","",infile,fixed=FALSE), "-mcmc-run-", id, sep='' )
outfile <- paste( basename, ".RData", sep='' )
plotfile <- paste( basename, "-plot",sep='')

if (logfile=="") { logfile <- paste(basename,".Rout",sep='') }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="output", split=interactive()) 
    if (!interactive()) { sink(file=logcon, type="message",append=TRUE) }
}

load(infile)
# if( !exists("which.taxa") ) { which.taxa <- c("human","chimpanzee") }

########

# Inference.
longwin <- leftwin+shortwin+rightwin
oldwin <- c(leftwin,shortwin,rightwin)
load(gmfile)   # contains leftwin, shortwin, rightwin: check compatability
stopifnot( all( oldwin == c(leftwin,shortwin,rightwin) ) )
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

countfile <- paste(paste("countdata/counts",leftwin,shortwin,rightwin,sep="-"),".RData",sep="")
load(countfile)

# initial parameters
initparams <- rowMeans( sapply( frame.mrun, "[[", "final" ) )

# set up root distribution
nmuts <- length(mutpats)
nfreqs <- length(bases)
npats <- nrow(genmatrix)

#####
# set up mcmc
require(mcmc)
patfreqmat <- log(initparams[1+nmuts+(1:(nfreqs-1))])[patcomp]
dim(patfreqmat) <- dim(patcomp)
lud <- function (params,which.frame) {
    # params are: tlen[1]/sum(tlen), sum(tlen)*mutrates, initfreqs[-length(initfreqs)]
    branchlens <- c(params[1],1-params[1])
    mutrates <- params[1+(1:nmuts)]
    initfreqs <- params[1+nmuts+(1:(nfreqs-1))]
    initfreqs <- c(initfreqs,1-sum(initfreqs))
    patfreqmat[] <- log(initfreqs)[patcomp]
    patfreqs <- exp( colSums( patfreqmat ) )
    if (any(mutrates<0) | any(initfreqs<0) | any(branchlens<0) ) {
        return( -Inf )
    } else {
        # only do in one direction... ?
        updownbranch <- getupdowntrans( genmatrix, projmatrix, mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=rev(branchlens) )
        # return (positive) log-posterior
        return( 
                sum( counts[[ which.taxa[[1]] ]][[ which.taxa[[2]] ]][[which.frame]] * log(updownbranch) )
                + sum( (tpriors-1)*log(branchlens) )
                - sum(mmeans*mutrates) 
                + sum( (ppriors-1)*log(initfreqs) )
            )
    }
}

## SEEMS TO WORK WELL WITHOUT THIS
# if (length(stepscale) == 1) {
#     # steps for tlen and freqs should be larger than for mutations
#     stepscale <- stepscale * c( rel.tlen=1, 6e6*rep(1e-8,nmuts)/30, rep(1,nfreqs-1) )  # reasonable for hu-ch?
# }

frame.mrun <- mclapply( seq_along( counts[[1]][[1]] ), function (which.frame) {
            metrop( lud, initial=initparams, nbatch=nbatches, blen=blen, scale=stepscale, which.frame=which.frame )
        }, mc.cores=3 )

save( opt, which.taxa, leftwin, shortwin, rightwin, boundary, meanboundary, mmeans, ppriors, tpriors, patcomp, gmfile, frame.mrun, file=outfile )

if (FALSE) {
    # rejigger this
pdf(file=paste(plotfile,".pdf",sep=''),width=6, height=4, pointsize=10)
subs <- seq.int(1, nrow(mrun$batch), by=max(1,floor(nrow(mrun$batch)/1000)) )
matplot( mrun$batch[subs,], type='l', col=1:length(truth) )
abline(h=truth[-length(truth)], col=adjustcolor(1:length(truth),.5), lwd=2)
abline(h=estimates["ans",], col=1:(length(truth)-1) )
legend("topright", col=c(1:length(truth),1,adjustcolor(1,.5)), lwd=c(rep(1,length(truth)-1),1,2),legend=c(colnames(estimates)[1:(length(truth)-1)],"point estimate","truth"))
dev.off()

pdf(file=paste(plotfile,"-pairwise.pdf",sep=''),width=6,height=6,pointsize=10)
mutlabels <- c("branchlen", paste("mut:", unlist( sapply( sapply( mutpats, lapply, paste, collapse="->" ), paste, collapse=" | " ) ) ), names(initfreqs) )
subs <- seq.int(1, nrow(mrun$batch), by=max(1,floor(nrow(mrun$batch)/1000)) )
pairs( rbind( mrun$batch[subs,], truth[-length(truth)] ), col=c( rep(adjustcolor("black",.1),length(subs)), adjustcolor("red",1)), pch=20, cex=c( rep(.25,length(subs)), 2), labels=mutlabels, gap=.1 )
dev.off()
}
