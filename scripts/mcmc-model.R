#!/usr/bin/env Rscript
library(optparse)

invocation <- commandArgs()

usage <- "\
Sample from the posterior on the parameters given data,\
beginning from a model fit already via likelihood or previous MCMC run. \
\
Config file gives prior means on model parameters, for instance:\
{ mutrates.prior: [ .01 ] } \
\
Also, the scale on which the MCMC tries to move around.  These are by default divided by 20, so the same scaling parameters can be used for the MCMC as for the optimization.
"

option_list <- list(
    # input/output
        make_option( c("-i","--infile"), type="character", help="Input file with previously fit 'context' model object."),
        make_option( c("-o","--outfile"), type="character", help="File to save results to.  [default: base of infile + 'mcmc' + jobid + .RData]"),
        make_option( c("-u","--basedir"), type="character", default=NULL, help="Directory to put output in. [default: same as infile]"),
        make_option( c("-c","--configfile"), type="character", help="JSON config file giving prior parameters."),
        make_option( c("-t","--tlen"), type="numeric", default=1, help="Guess at time quantity to scale initial values of mutation parameters by. [default=%default]"),
        make_option( c("-s","--scalefac"), type="numeric", default=.05, help="Multiply the scale factors in the config file by this much for the MCMC steps. [default=%default]"),
        make_option( c("-b","--nbatches"), type="integer", default=100, help="Number of MCMC batches to run for (results will be means of each batch). [default=%default]"),
        make_option( c("-l","--blen"), type="integer", default=1, help="Length of each MCMC batch. [default=%default]"),
        make_option( c("-j","--jobid"), type="character", default=formatC(1e6*runif(1),width=6,format="d",flag="0"), help="Unique job id. [default random]")
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$infile)) { stop("No input file.  Run\n  Rscript mcmc-model.R -h\n for help.\n") }
if (is.null(opt$basedir)) { opt$basedir <- dirname(opt$infile) }
if (is.null(opt$outfile)) { opt$outfile <- paste( opt$basedir, "/", gsub("(-mcmc-[0-9]*)*\\.[^.]*","",basename(opt$infile) ), "-mcmc-", opt$jobid, ".RData", sep='' ) }
print(opt) # this will go in the pbs log

library(jsonlite)
library(mcmc)

library(contextual)
library(contextutils)

options(error = print.and.dump)

# load previously fit model
stopifnot(file.exists(opt$infile))
load(opt$infile)  # provides 'model'

genmatrix <- model@genmatrix
projmatrix <- model@projmatrix
counts <- model@counts

# shorter counts for full likelihood
stopifnot(shortwin(counts) > 1)
counts.0 <- projectcounts(counts, new.shortwin=shortwin(counts) - 1)
projmatrix.0 <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin(counts), fpatterns=colnames(counts.0) )

# read in config file
prior.config <- read.config(opt$configfile)  # returns NULL if not present
if (is.null(prior.config$mutrates.prior)) { prior.config$mutrates.prior <- model@mutrates.prior }
if (is.null(prior.config$selcoef.prior)) { 
    if (nsel(model)==0) { 
        prior.config$selcoef.prior <- numeric(0) 
    } else  { 
        prior.config$selcoef.prior <- model@selcoef.prior 
    }
}
if (is.null(prior.config$fixfn.params.prior)) { 
    if (length(fixparams(model))==0) { 
        prior.config$fixfn.params.prior <- numeric(0) 
    } else  { 
        prior.config$fixfn.params.prior <- model@fixfn.params.prior 
    }
}
prior.config <- fill.default.config( prior.config )

if (length(prior.config$mutrates) != length(model@mutrates) 
        || length(prior.config$selcoef) != length(model@selcoef) 
        || length(prior.config$fixfn.params) != length(model@params) ) {
    stop(sprintf("Configuration in %s does not match that in already-fit model of %s .", opt$configfile, opt$infile))
}

# scale tuning parameters
if (FALSE) {  # this leads to trouble
    if (is.null(prior.config$mutrates.scale)) {
        prior.config$mutrates.scale <- 1e-3 * rep( mean(prior.config$mutrates),nmuts(genmatrix) ) 
    }
    if (is.null(prior.config$selcoef.scale)) {
        prior.config$selcoef.scale <- rep(.001,nsel(genmatrix))
    }
    if (is.null(prior.config$fixfn.params.scale)) {
        prior.config$fixfn.params.scale <- 0.1*prior.config$fixfn.params
    }
}
parscale <- with(prior.config, unlist( c(mutrates.scale*opt$tlen, selcoef.scale, fixfn.params.scale) ) )
names(parscale) <- c( mutnames(genmatrix@mutpats), selnames(genmatrix@selpats), fixparams(genmatrix) )
parscale <- parscale * opt$scalefac

# skip these parameters
use.par <- ( parscale!=0 )
params <- coef(model)
initpar <- params[use.par]

# The log-posterior function.
likfun <- function (sub.params){
    # params are: mutrates*tlen
    params[use.par] <- sub.params
    mutrates <- params[1:nmuts(genmatrix)]
    selcoef <- params[seq( 1+nmuts(genmatrix), length.out=nsel(genmatrix) )]
    fparams <- params[seq( 1+nmuts(genmatrix)+nsel(genmatrix), length.out=length(fixparams(genmatrix)) )]
    names(fparams) <- fixparams(model)
    if (any(mutrates<0)) { return( -Inf ) }
    genmatrix@x <- do.call( update_x, 
                           c( list(
                                   G=genmatrix, 
                                   mutrates=mutrates, 
                                   selcoef=selcoef ), 
                             as.list(fparams) ) )
    # this is collapsed transition matrix
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, tlen=1, time="fixed") # shape=params[length(params)], time="gamma" )
    subtransmatrix.0 <- computetransmatrix( genmatrix, projmatrix.0, tlen=1, time="fixed") # shape=params[length(params)], time="gamma" )
    # return POSITIVE log-likelihood
    num <- sum( counts@counts * log(subtransmatrix) )
    den <- sum( counts.0@counts * log(subtransmatrix.0) )
    ans <- num - den
    if (!is.finite(ans)) { return( -Inf ) }
    else { return( ans - sum(mutrates/prior.config$mutrates.prior) - sum((selcoef/prior.config$selcoef.prior)^2) - sum((fparams/prior.config$fixfn.params.prior)^2) ) }
}

baseval <- likfun(initpar)
stopifnot( is.finite(baseval) )


mrun <- metrop( likfun, initial=initpar, nbatch=opt$nbatches, blen=opt$blen, scale=parscale[use.par] )

mrun$use.par <- use.par
mrun$parscale <- parscale
mrun$initpar <- initpar
mrun.final.par <- params
mrun.final.par[use.par] <- mrun$final


model <- new( "contextMCMC",
             counts=model@counts,
             genmatrix=model@genmatrix,
             projmatrix=model@projmatrix,
             mutrates=mrun.final.par[1:nmuts(genmatrix)],
             selcoef=mrun.final.par[seq(1+nmuts(genmatrix),length.out=nsel(genmatrix))],
             params=mrun.final.par[seq(1+nmuts(genmatrix)+nsel(genmatrix),length.out=length(fixparams(model)))],
             results=unclass(mrun),
             likfun=likfun,
             mutrates.prior=prior.config$mutrates.prior,
             selcoef.prior=prior.config$selcoef.prior,
             fixfn.params.prior=prior.config$fixfn.params.prior,
             invocation=invocation
         )

# set this up so that we can call likfun again in the future, directly
likfun.env <- new.env()
assign("use.par",use.par,envir=likfun.env)
assign("params",params,envir=likfun.env)
assign("prior.config",prior.config,envir=likfun.env)
assign("genmatrix",model@genmatrix,envir=likfun.env)
assign("projmatrix",model@projmatrix,envir=likfun.env)
assign("counts",model@counts,envir=likfun.env)
assign("projmatrix.0",projmatrix.0,envir=likfun.env)
assign("counts.0",counts.0,envir=likfun.env)
environment(model@likfun) <- likfun.env



save(model,file=opt$outfile)

