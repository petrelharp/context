#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Infer parameters from output of sim-tree-cpg.R .\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default="mtCDNApri-sub.nuc", help=".nuc file with data. [default \"%default\"]"),
        make_option( c("-f","--basedir"), type="character", default="", help="Directory for output. [default infile with suffix removed]"),
        make_option( c("-w","--win"), type="integer", default=2, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--lwin"), type="integer", default=1, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rwin"), type="integer", default=1, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-x","--maxit"), type="integer", default=0, help="Maximum number of iterates in optim for point estimate. [default \"%default\"]" ),
        make_option( c("-n","--nbatches"), type="integer", default=200, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=1, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-s","--stepscale"), type="numeric", default=1e-4, help="Scale of proposal steps for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-m","--mmean"), type="double", default=1, help="Prior mean on triplet time-scaled mutation rates. [default \"%default\"]" ),
        make_option( c("-c","--cpgmean"), type="double", default=1, help="Prior variance on CpG time-scaled rate. [default \"%default\"]" ),
        make_option( c("-p","--pprior"), type="double", default=1, help="Parameter for Dirichlet prior on base frequencies. [default \"%default\"]" ),
        make_option( c("-v","--tprior"), type="double", default=.5, help="Parameter for Beta prior on branch length. [default \"%default\"]" ),
        make_option( c("-d","--boundary"), type="character", default="none", help="Boundary conditions for generator matrix. [default \"%default\"]"),
        make_option( c("-y","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases in computing generator matrix. [default \"%default\"]" ),
        make_option( c("-z","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix, or TRUE [default] to look for one. (otherwise, will compute)"),
        make_option( c("-e","--countfile"), type="character", default="", help="Record tuple counts in this file. [default: countdata/counts-(lwin)-(win)-(rwin).RData]" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)

winlen <- lwin+win+rwin

if (gmfile=="TRUE") { gmfile <- paste(paste("genmatrices/genmatrix",winlen,boundary,meanboundary,sep="-"),".RData",sep='') }

if (countfile=="") { countfile <- paste(paste("countdata/counts",lwin,win,rwin,sep="-"),".RData",sep="") }

if (basedir=="") { basedir <- gsub(".nuc","",infile,fixed=TRUE) }


# took mtCDNApri.nuc from paml/examples/mtCDNA
#  and: removed extra spaces on first line (only one space between numbers)
#       made each sequence be on a single line
#       removed trailing text

require(Biostrings)

scriptdir <- "../"
source(paste(scriptdir,"codon-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))

if (!file.exists(basedir)) { dir.create(basedir) }
basename <- paste(basedir,"/win-",lwin,"-",win,"-",rwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')

if (logfile=="") {
    logfile <- gsub(".RData",".Rout",datafile,fixed=TRUE)
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message", split=interactive()) 
    sink(file=logcon, type="output", split=interactive())   # send both to log file
}


if (file.exists(gmfile)) {
    load(gmfile)
} else {
    stop("Precompute generator matrix file.")
}
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

# all pairwise counts
if (!file.exists(countfile)) {
    mtCDNA <- as( readDNAMultipleAlignment(infile, format="phylip" ), "DNAStringSet" )
    mtCDNA <- mtCDNA[ c("human", "chimpanzee", "bonobo", "gorilla") ]

    require(parallel)
    counts <- mclapply( names( mtCDNA ), function (x) {
            x <- lapply( names( mtCDNA ), function (y) {
                    counttrans( rownames(projmatrix), colnames(projmatrix), mtCDNA[[x]],  mtCDNA[[y]],  lwin=lwin, shift=3 ) 
                } ) 
            names(x) <- names( mtCDNA )
            return(x)
        }, mc.cores=6 )
    names(counts) <- names( mtCDNA )
    if (!file.exists(dirname(countfile))) { dir.create(dirname(countfile)) }
    save( counts, file=countfile )
} else {
    load( countfile )
}

# not used in likelihood:
initcounts <- lapply( lapply(counts,"[[",1), lapply, rowSums )  # initial counts for (x,y) only depends on x

# set up root distribution
nmuts <- length(mutpats)
nfreqs <- length(bases)
npats <- nrow(genmatrix)
patcomp <- apply( do.call(rbind, strsplit(rownames(genmatrix),'') ), 2, match, bases )  # which base is at each position in each pattern
patcomp <- t(patcomp) # column-wise more efficient

which.taxa <- c("human","chimpanzee")
# which.frame will denote which element of counts[[x]][[y]] to use

# composite likelihood with taxa x->y + y->x
likfun <- function (params,which.frame) {
    # params are: tlen[1]/sum(tlen), sum(tlen)*mutrates, initfreqs
    branchlens <- c(params[1],1-params[1])
    mutrates <- params[1+(1:nmuts)]
    initfreqs <- params[1+nmuts+(1:nfreqs)]
    initfreqs <- initfreqs/sum(initfreqs)
    patfreqs <- log(initfreqs)[patcomp]
    dim(patfreqs) <- dim(patcomp)
    patfreqs <- exp( colSums( patfreqs ) )
    # these are collapsed transition matrix
    updownbranch <- list(  # note "up" branch is from simpler summaries
            getupdowntrans( genmatrix, projmatrix, mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=rev(branchlens) ),
            getupdowntrans( genmatrix, projmatrix, mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=branchlens )
        )
    # if (any(sapply(updownbranch,function(x) any(!is.numeric(x))|any(x<0)))) { browser() }
    # return negative log-likelihood plus a penalty to keep initfreqs summing to (almost) 1
    return( 
                (-1) * ( 
                    sum( counts[[which.taxa[[1]]]][[which.taxa[[2]]]][[which.frame]] * log(updownbranch[[1]]) ) + 
                    sum( counts[[which.taxa[[2]]]][[which.taxa[[1]]]][[which.frame]] * log(updownbranch[[2]]) ) ) +
                100*(sum(initfreqs)-1)^2
            )
}

initparams <- c( rel.tlen=0.5, 6e6*rep(1e-8,nmuts)/30, rep(1/nfreqs,nfreqs) )  # reasonable for hu-ch?
stopifnot( all( sapply( seq_along( counts[[1]][[1]] ), function (which.frame) { is.finite( likfun(initparams,which.frame=which.frame) ) } ) ) )

lbs <- c( 1e-6, rep(0,nmuts), rep(1e-6,nfreqs) )
ubs <- c( 1, rep(20,nmuts), rep(1,nfreqs) )

if (maxit>0) {
    frame.ans <- mclapply( seq_along( counts[[1]][[1]] ), function (which.frame) {
                ans <- optim( par=initparams, fn=likfun, which.frame=which.frame, method="L-BFGS-B", lower=lbs, upper=ubs, 
                        control=list(fnscale=abs(likfun(initparams,which.frame=1)), maxit=maxit, parscale=initparams ) )
             }, mc.cores=3 )
} else  {
    frame.ans <- NULL
}


#####
# set up mcmc
require(mcmc)
mmeans <- c( rep(mmean,nmuts-1), cpgmean )
ppriors <- rep( pprior, nfreqs )
tpriors <- rep( tprior, 2 )
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

# bayesian
#  note we deal with initial freqs not summing to 1 differently from in optim( ) -- need to remove last entry.
frame.mrun <- mclapply( seq_along( counts[[1]][[1]] ), function (which.frame) {
            metrop( lud, initial=initparams[-length(initparams)], nbatch=nbatches, blen=blen, scale=stepscale, which.frame=which.frame )
        }, mc.cores=3 )

save( opt, which.taxa, lwin, rwin, win, winlen, boundary, meanboundary, mmeans, ppriors, tpriors, nmuts, nfreqs, npats, patcomp, gmfile, projmatrix, subtransmatrix, counts, initcounts, frame.ans, frame.mrun, file=datafile )
