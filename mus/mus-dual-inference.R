#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Infer mutation rates.  Allow different parameters on each branch.
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

basedir <- paste(infile,"-dual-results",sep='')
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
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

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
adhoc <- lapply(counts, countmuts,mutpats=mutpats,lwin=lwin)
adhoc <- lapply( adhoc, function (x) x[1,]/x[2,] )

# move from base frequencies (what we estimate) to pattern frequencies
nmuts <- length(mutpats)
nfreqs <- length(bases)
npats <- nrow(genmatrix)
patcomp <- apply( do.call(rbind, strsplit(rownames(genmatrix),'') ), 2, match, bases )  # which base is at each position in each pattern
likfun <- function (params) {
    # params are: sum(tlen)*mutrates2, sum(tlen)*mutrates2, , initfreqs
    branchlens <- c(.5,.5)
    mutrates <- list( params[(1:nmuts)],  params[nmuts+(1:nmuts)] )
    initfreqs <- params[2*nmuts+(1:nfreqs)]
    initfreqs <- initfreqs/sum(initfreqs)
    patfreqs <- initfreqs[patcomp]
    dim(patfreqs) <- dim(patcomp)
    patfreqs <- apply( patfreqs, 1, prod )
    # these are collapsed transition matrix
    updownbranch <- list(  # note "up" branch is from simpler summaries
            getupdowntrans( genmatrix, projmatrix, mutrates=mutrates, selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=rev(branchlens) ),
            getupdowntrans( genmatrix, projmatrix, mutrates=mutrates, selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=branchlens )
        )
    # if (any(sapply(updownbranch,function(x) any(!is.numeric(x))|any(x<0)))) { browser() }
    # return negative log-likelihood plus a penalty to keep initfreqs summing to (almost) 1
    return( 
                (-1) * ( sum( counts[[1]] * log(updownbranch[[1]]) ) 
                        + sum( counts[[2]] * log(updownbranch[[2]]) ) ) 
                + 100*(sum(initfreqs)-1)^2
            )
}

# only nonoverlapping counts, plus priors -- indep't trials. (multinomial)
mmeans <- c( rep(mmean,nmuts-1), cpgmean )
ppriors <- rep( pprior, nfreqs )
lud <- function (params) {
    # params are: tlen[1]/sum(tlen), sum(tlen)*mutrates1, sum(tlen)*mutrates2, initfreqs[-length(initfreqs)]
    branchlens <- c(.5,.5)
    mutrates <- list( params[(1:nmuts)],  params[nmuts+(1:nmuts)] )
    initfreqs <- params[nmuts+(1:(nfreqs-1))]
    initfreqs <- c(initfreqs,1-sum(initfreqs))
    patfreqs <- initfreqs[patcomp]
    dim(patfreqs) <- dim(patcomp)
    patfreqs <- apply( patfreqs, 1, prod )
    if (any(unlist(mutrates)<0) | any(initfreqs<0)) {
        return( -Inf )
    } else {
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
}

# MLE point estimates
#   params are: tlen[1]/sum(tlen), sum(tlen)*mutrates, initfreqs
initpar <- c( adhoc[[1]], adhoc[[2]], rep(.25,4) ) # random init
names(initpar) <- c( 
        paste( unlist( sapply( sapply( mutpats, lapply, paste, collapse="->", sep='' ), paste, collapse="|", sep='' ) ), ".ab", sep='' ),
        paste( unlist( sapply( sapply( mutpats, lapply, paste, collapse="->", sep='' ), paste, collapse="|", sep='' ) ), ".ab", sep='' ),
        paste("init",bases,sep='.') )
lbs <- c( rep(0,2*nmuts), rep(1e-6,length(bases)) )
ubs <- c( rep(20,2*nmuts), rep(1,length(bases)) )
# # do in parallel:
mle <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(trace=3,fnscale=abs(likfun(initpar)),parscale=rep(optimscale,length(initpar)),maxit=500) )
# renormalize the initial frequencies
mle.par <- mle$par; mle.par[2*nmuts+(1:nfreqs)] <- mle.par[2*nmuts+(1:nfreqs)] / sum( mle.par[2*nmuts+(1:nfreqs)] )

estimates <- data.frame( rbind(init=initpar, mle=mle.par) )
colnames(estimates) <- c( 
        paste( unlist( sapply( sapply( mutpats, lapply, paste, collapse="->", sep='' ), paste, collapse="|", sep='' ) ), ".ab", sep='' ),
        paste( unlist( sapply( sapply( mutpats, lapply, paste, collapse="->", sep='' ), paste, collapse="|", sep='' ) ), ".ab", sep='' ),
        paste("init",bases,sep='.') )
estimates$likfun <- apply( estimates, 1, likfun )
write.table( estimates, file=resultsfile, quote=FALSE, sep="\t" )

# bayesian
#  note we deal with initial freqs not summing to 1 differently from in optim( ) -- need to remove last entry.
#  mrun.parjob <- mcparallel( metrop( lud, initial=mle.par[-length(mle.par)], nbatch=nbatches, blen=blen, scale=stepscale ) )
#  mrun <- mccollect(mcrun.parjob)
mrun <- metrop( lud, initial=mle.par[-length(mle.par)], nbatch=nbatches, blen=blen, scale=stepscale )

save( opt, counts, genmatrix, projmatrix, subtransmatrix, lud, likfun, mle, estimates, initpar, mmeans, ppriors, mrun, win, lwin, rwin, nmuts, nfreqs, npats, patcomp, file=datafile )

print(format(Sys.time(),"%Y-%m-%d-%H-%M"))
