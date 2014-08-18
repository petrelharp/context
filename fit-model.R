#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Infer parameters from paired counts file, which records instances of Tmer transitions.
"

option_list <- list(
    # input/output
        make_option( c("-i","--infile"), type="character", help="Input file with tuple counts, tab-separated, with header 'reference', 'derived', 'count'. [default, looks in basedir]" ),
        make_option( c("-o","--outfile"), type="character", help="File to save results to.  [default: base of infile + base of genmatrix + jobid + .RData]"),
        make_option( c("-u","--basedir"), type="character", default=NULL, help="Directory to put output in. [default: same as infile]"),
        make_option( c("-l","--leftwin"), type="integer", help="Size of left-hand context." ),
        make_option( c("-m","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix."),
        make_option( c("-x","--maxit"), type="integer", default=100, help="Number of iterations of optimization to run for. [default=%default]"),
        make_option( c("-j","--jobid"), type="character", default=formatC(1e6*runif(1),width=6,format="d",flag="0"), help="Unique job id. [default random]")
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$infile)) { stop("No input file.  Run\n  bcells-inference.R -h\n for help.\n") }
if (is.null(opt$basedir)) { opt$basedir <- dirname(opt$infile) }
if (is.null(opt$outfile)) { opt$outfile <- paste( opt$basedir, "/", gsub("\\.[^.]*","",basename(opt$infile) ), "-", gsub("\\.[^.]*","",basename(opt$gmfile) ), "-", opt$jobid, ".RData", sep='' ) }
print(opt) # this will go in the pbs log
options(error = quote({dump.frames(to.file = TRUE); q()}))

source("../context-inference-fns.R")

attach(opt)

# load generator matrix
stopifnot(file.exists(opt$gmfile))
load(opt$gmfile)  # provides 'genmatrix'

# read in counts
counts <- read.counts(opt$infile, opt$leftwin, bases=genmatrix@bases, longpats=rownames(genmatrix) )
stopifnot( all( rownames(counts) == rownames(genmatrix) ) )

projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin(counts), fpatterns=colnames(counts) )

# get ad hoc initial guesses at parameters
adhoc <- countmuts(counts=counts,mutpats=mutpats,leftwin=leftwin(counts))
adhoc.mutrates <- adhoc[1,]/adhoc[2,]
adhoc.mutrates <- ifelse( is.finite(adhoc.mutrates) & adhoc.mutrates > 0, adhoc.mutrates, 1e-4 )


# Compute (quasi)-likelihood function using all counts -- multinomial as described in eqn:comp_like.
likfun <- function (params){
    # params are: mutrates*tlen, shape
    genmatrix@x <- update(genmatrix,mutrates=params[1:nmuts(genmatrix)],selcoef=params[seq(1+nmuts(genmatrix),length.out=nsel(genmatrix))])
    # this is collapsed transition matrix
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, tlen=1, time="fixed") # shape=params[length(params)], time="gamma" )
    # return POSITIVE log-likelihood
    ans <- sum( counts@counts * log(subtransmatrix) )
    if (!is.finite(ans)) print(paste("Warning: non-finite likelihood with params:",params))
    return(ans)
}

initpar <- adhoc.mutrates
stopifnot( length(initpar) == nmuts(genmatrix)+nsel(genmatrix) )
lbs <- rep(1e-6,nmuts(genmatrix))
ubs <- rep(2,nmuts(genmatrix))
parscale <- 1e-3 * rep(mean(adhoc.mutrates),length(adhoc.mutrates))

baseval <- likfun(initpar)
stopifnot( is.finite(baseval) )
optim.results <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(fnscale=(-1)*abs(baseval), parscale=parscale, maxit=maxit) )

model <- new( "context",
             counts=counts,
             genmatrix=genmatrix,
             projmatrix=projmatrix,
             mutrates=adhoc.mutrates,
             selcoef=numeric(nsel(genmatrix)),
             params=numeric(0),
             results=optim.results,
             likfun=likfun
         )

save(model,file=outfile)
