#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Run mcmc longer.\
"

option_list <- list(
        make_option( c("-u","--indir"), type="character", default=NULL, help="Directory to look for input in, and write output files to." ),
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file containing simulation." ),
        make_option( c("-v","--revfile"), type="character", default=NULL, help="Table of count data in the reverse orientation"),
        make_option( c("-j","--jobid"), type="character", default=formatC(1e6*runif(1),width=6,format="d",flag="0"), help="Unique job id. [default random]"),
        make_option( c("-w","--shortwin"), type="integer", default=1, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--leftwin"), type="integer", default=1, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rightwin"), type="integer", default=1, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-m","--mmean"), type="double", default=1, help="Prior mean on single base mutation rates. [default \"%default\"]" ),
        make_option( c("-c","--cpgmean"), type="double", default=1, help="Prior variance on CpG rate. [default \"%default\"]" ),
        make_option( c("-p","--pprior"), type="double", default=1, help="Parameter for Dirichlet prior on base frequencies. [default \"%default\"]" ),
        make_option( c("-t","--tprior"), type="double", default=.5, help="Parameter for Beta prior on branch length. [default \"%default\"]" ),
        make_option( c("-n","--nbatches"), type="integer", default=1000, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-z","--stepscale"), type="character", default="1e-4", help="Scale of proposal steps for mutation*time parameters for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-e","--initscale"), type="character", default="3e-3", help="Scale of proposal steps for initial frequencies for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-s","--restart"), action="store_true", default=FALSE, help="Start a whole new MCMC run?" ),
        make_option( c("-d","--boundary"), type="character", default="none", help="Boundary conditions for generator matrix. [default \"%default\"]"),
        make_option( c("-y","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases in computing generator matrix. [default \"%default\"]" ),
        make_option( c("-g","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix, or TRUE [default] to look for one. (otherwise, will compute)"),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.character(opt$stepscale)) { opt$stepscale <- eval(parse(text=opt$stepscale)) }
if (is.character(opt$initscale)) { opt$initscale <- eval(parse(text=opt$initscale)) }
if (is.null(opt$infile) & is.null(opt$indir)) { stop("No input file.  Run\n  bcells-inference.R -h\n for help.\n") }
attach(opt)

if (interactive()) { nbatches <- 100; blen <- 10; restart <- FALSE }

scriptdir <- "../"
source(paste(scriptdir,"context-inference-fns.R",sep=''))
# source(paste(scriptdir,"sim-context-fns.R",sep=''))
require(mcmc)

longwin <- leftwin+shortwin+rightwin

if (substr(indir,nchar(indir),nchar(indir)) %in% c("/","\\")) { indir <- substr(indir,1,nchar(indir)-1) }
if (is.null(opt$infile)) { infile <- paste(indir,"/", longwin,".",shortwin,".counts",sep='') }
if (is.null(opt$revfile)) { revfile <- paste(dirname(infile),"/rev.",basename(infile),sep='') }
if (!file.exists(infile) | !file.exists(revfile)) { stop("Cannot read file ", infile) }

basedir <- paste(infile,"-dual-results",sep='')
if (!file.exists(basedir)) { dir.create(basedir) }
basename <- paste(basedir,"/win-",shortwin,"-",leftwin,"-",rightwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
plotfile <- paste( basename ,"-plot",sep='')
mcmcdatafiles <- list.files(path=basedir,pattern="-mcmc.*RData",full.names=TRUE)
mcmcnum <- 1+max(c(0,as.numeric(gsub(".*-mcmc-","",gsub(".RData","",mcmcdatafiles)))),na.rm=TRUE)

if (logfile=="") { logfile <- paste(basename,"-mcmc-run-",mcmcnum,".Rout",sep='') }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="output", split=interactive()) 
    if (!interactive()) { sink(file=logcon, type="message",append=TRUE) }
}

load(datafile)  # has mrun and previous things
if (length(mcmcdatafiles)>0) { load(grep(paste("-mcmc-",mcmcnum-1,".RData",sep=''),mcmcdatafiles,fixed=TRUE,value=TRUE)) }

########

# Inference.
if (gmfile=="TRUE") { gmfile <- paste(paste("genmatrices/genmatrix",longwin,boundary,meanboundary,sep="-"),".RData",sep='') }
if (file.exists(gmfile)) {
    load(gmfile) 
} else {
    stop("Cannot find gmfile.")
}
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin )

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

# move from base frequencies (what we estimate) to pattern frequencies
nmuts <- length(mutpats)
nfreqs <- length(bases)
npats <- nrow(genmatrix)
patcomp <- apply( do.call(rbind, strsplit(rownames(genmatrix),'') ), 2, match, bases )  # which base is at each position in each pattern
mmeans <- c( rep(mmean,nmuts-1), cpgmean )
ppriors <- rep( pprior, nfreqs )
lud <- function (params) {
    # params are: tlen[1]/sum(tlen), sum(tlen)*mutrates1, sum(tlen)*mutrates2, initfreqs[-length(initfreqs)]
    branchlens <- c(.5,.5)
    mutrates <- list( params[(1:nmuts)],  params[nmuts+(1:nmuts)] )
    initfreqs <- params[2*nmuts+(1:(nfreqs-1))]
    initfreqs <- c(initfreqs,1-sum(initfreqs))
    if (any(unlist(mutrates)<0) | any(initfreqs<0)) { return( -Inf ) } 
    patfreqs <- initfreqs[patcomp]
    dim(patfreqs) <- dim(patcomp)
    patfreqs <- apply( patfreqs, 1, prod )
    updownbranch <- list(  # note "up" branch is from simpler summaries
            getupdowntrans( genmatrix, projmatrix, mutrates=mutrates, selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=rev(branchlens) ),
            getupdowntrans( genmatrix, projmatrix, mutrates=mutrates, selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=branchlens )
        )
    # return (positive) log-posterior
    return( 
            # (-1)*sum(updownbranch[nonoverlapping[[1]]]) +
            sum( counts[[1]] * log(updownbranch[[1]]) ) 
            + sum( counts[[2]] * log(updownbranch[[2]]) ) 
            - sum(mmeans*unlist(mutrates)) 
            + sum( (ppriors-1)*log(initfreqs) )
        )
}

# simple point estimates for starting positions
adhoc <- lapply(counts, countmuts,mutpats=mutpats,leftwin=leftwin)
adhoc <- unlist( sapply( adhoc, function (x) x[1,]/x[2,] ) )

# construct simplified proposal distribution for mcmc steps:
scalemat <- c( rep_len( stepscale, nmuts ), rep_len( stepscale, nmuts ), rep_len( initscale, nfreqs-1 ) )

if (restart) {
    # mrun <- metrop( lud, initial=random.ans$par[-length(random.ans$par)], nbatch=nbatches, blen=blen, scale=scalemat )
    # if (is.finite(lud(mle$par[-length(mle$par)]))) {
        mrun <- metrop( lud, initial=mle$par[-length(mle$par)], nbatch=nbatches, blen=blen, scale=scalemat )
    # } else {  # the point estimate broke; cheat to get a good starting point
    #     mrun <- metrop( lud, initial=truth, nbatch=nbatches, blen=blen, scale=scalemat )
    # }
} else {
    mrun <- metrop( mrun, nbatch=nbatches, blen=blen, scale=scalemat )
}
colnames(mrun$batch) <- colnames(estimates[1:ncol(mrun$batch)])

save( leftwin, shortwin, rightwin, lud, mrun, file=paste(basename,"-mcmc-",mcmcnum,".RData",sep='') )

npar <- ncol(estimates)-1
mlepar <- estimates["mle",][1:npar]

pdf(file=paste(plotfile,"-mcmc-",mcmcnum,".pdf",sep=''),width=10, height=8, pointsize=10)
subs <- seq.int(1, nrow(mrun$batch), by=max(1,floor(nrow(mrun$batch)/1000)) )
matplot( mrun$batch[subs,], type='l', col=1:npar )
abline(h=mlepar, col=adjustcolor(1:npar,.5), lwd=2)
dev.off()

pdf(file=paste(plotfile,"-mcmc-",mcmcnum,"-pairwise.pdf",sep=''),width=10,height=8,pointsize=10)
mutlabels <- colnames(estimates)[-npar]
subs <- seq.int(1, nrow(mrun$batch), by=max(1,floor(nrow(mrun$batch)/1000)) )
pairs( rbind( mrun$batch[subs,], mlepar[-npar] ), col=c( rep(adjustcolor("black",.1),length(subs)), adjustcolor("red",1)), pch=20, cex=c( rep(.25,length(subs)), 2), labels=mutlabels, gap=.1 )
dev.off()
