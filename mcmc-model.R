#!/usr/bin/Rscript --vanilla
require(optparse)
require(jsonlite)
require(mcmc)

usage <- "\
Sample from the posterior on the parameters given data,\
beginning from a model fit already via likelihood or previous MCMC run. \
\
Config file gives prior means on model parameters, for instance:\
{ mutpriors: [ .01 ] } \
"

option_list <- list(
    # input/output
        make_option( c("-i","--infile"), type="character", help="Input file with tuple counts, tab-separated, with header 'reference', 'derived', 'count'. [default, looks in basedir]" ),
        make_option( c("-o","--outfile"), type="character", help="File to save results to.  [default: base of infile + 'mcmc' + jobid + .RData]"),
        make_option( c("-u","--basedir"), type="character", default=NULL, help="Directory to put output in. [default: same as infile]"),
        make_option( c("-c","--priorfile"), type="character", help="JSON config file giving prior parameters."),
        make_option( c("-b","--nbatches"), type="integer", default=100, help="Number of MCMC batches to run. [default=%default]"),
        make_option( c("-l","--blen"), type="integer", default=100, help="Length of each MCMC batch. [default=%default]"),
        make_option( c("-s","--stepscale"), type="numeric", help="Size of MCMC steps."),
        make_option( c("-j","--jobid"), type="character", default=formatC(1e6*runif(1),width=6,format="d",flag="0"), help="Unique job id. [default random]")
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$infile)) { stop("No input file.  Run\n  bcells-inference.R -h\n for help.\n") }
if (is.null(opt$basedir)) { opt$basedir <- dirname(opt$infile) }
if (is.null(opt$outfile)) { opt$outfile <- paste( opt$basedir, "/", gsub("\\.[^.]*","",basename(opt$infile) ), "-mcmc-", opt$jobid, ".RData", sep='' ) }
print(opt) # this will go in the pbs log
options(error = quote({dump.frames(to.file = TRUE); q()}))

source("../context-inference-fns.R")

# load generator matrix
stopifnot(file.exists(opt$infile))
load(opt$infile)  # provides 'model'

genmatrix <- model@genmatrix
projmatrix <- model@projmatrix
counts <- model@data

# read in config file
if (!is.null(opt$priorfile)) {
    prior <- fromJSON(opt$priorfile,simplifyMatrix=FALSE)
    mutprior <- prior$mutprior
    if (nsel(model)>0) { selprior <- prior$selprior } else { selprior <- numeric(0) }
} else { 
    mutprior <- model@mutprior
    selprior <- model@selprior
}

# (quasi)-log-posterior 
likfun <- function (params){
    # params are: mutrates*tlen
    if (any(params<0)) { return( -Inf ) }
    mutrates <- params[1:nmuts(genmatrix)]
    selcoef <- params[seq( nmuts(genmatrix), length.out=nsel(genmatrix) )]
    genmatrix@x <- update(genmatrix, mutrates=mutrates, selcoef=selcoef )
    # this is collapsed transition matrix
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, tlen=1, time="fixed") # shape=params[length(params)], time="gamma" )
    # return POSITIVE log-likelihood 
    ans <- sum( counts@counts * log(subtransmatrix) ) 
    if (!is.finite(ans)) print(params)
    return( ans - sum(mutrates/mutprior) - sum(selcoef/selprior) )
}

initpar <- coef(model)
baseval <- likfun(initpar)
stopifnot( is.finite(baseval) )
if (is.null(opt$stepscale)) { opt$stepscale <- mean(initpar)/100 }

mrun <- metrop( likfun, initial=initpar, nbatch=opt$nbatches, blen=opt$blen, scale=opt$stepscale )

model <- new( "contextMCMC",
             data=model@data,
             genmatrix=model@genmatrix,
             projmatrix=model@projmatrix,
             mutrates=mrun$final[1:nmuts(genmatrix)],
             selcoef=mrun$final[seq(nmuts(genmatrix),length.out=nsel(genmatrix))],
             params=numeric(0),
             results=as.list(mrun),
             likfun=likfun,
             mutprior=mutprior,
             selprior=selprior
         )


save(model,file=outfile)
