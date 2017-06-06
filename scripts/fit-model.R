#!/usr/bin/env Rscript
library(optparse)


invocation <- commandArgs()

usage <- "\
Infer parameters from paired counts file, which records instances of Tmer transitions.
\
The config file is a JSON file defining initial values for: \
    mutrates : numeric vector of same length as mutpats giving mutation rates per generation \
    selcoef : numeric vector of same length as selpats giving selection coefficients \
    fixfn.params: named list of additional parameters to pass to fixfn \
as well as the appropriate scale to search for parameter values over, by \
    mutrate.scale : numeric vector for mutrates \
    selcoef.scale : numeric vector for selcoef \
    fixfn.params.scale : numeric vector for fixfn.params \
\
Parameters whose scale is set to zero *will be regarded as fixed.* \
\
If any of these are missing, possibly inappropriate values will be guessed. \
"

option_list <- list(
    # input/output
        make_option( c("-i","--infile"), type="character", help="Input file with tuple counts, tab-separated, with header 'reference', 'derived', 'count'. [default, looks in basedir]" ),
        make_option( c("-c","--configfile"), type="character", help="Config file with initial guesses at parameter values and relevant scales for each."),
        make_option( c("-t","--tlen"), type="numeric", default=1, help="Guess at time quantity to scale initial values of mutation parameters by. [default=%default]"),
        make_option( c("-o","--outfile"), type="character", help="File to save results to.  [default: base of infile + base of genmatrix + jobid + .RData]"),
        make_option( c("-u","--basedir"), type="character", default=NULL, help="Directory to put output in. [default: same as infile]"),
        make_option( c("-m","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix."),
        make_option( c("-x","--maxit"), type="integer", default=100, help="Number of iterations of optimization to run for. [default=%default]"),
        make_option( c("-n","--nonnegative"), action="store_true", default=FALSE, help="Constrain mutation rates to be nonnegative? [default=%default]"),
        make_option( c("-j","--jobid"), type="character", default=formatC(1e6*runif(1),width=6,format="d",flag="0"), help="Unique job id. [default random]"),
        make_option( c("-z","--seed"), type="integer", help="Seed for pseudorandom number generator; an integer. [default: does not meddle]")
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$infile) || is.null(opt$configfile)) { stop("No input file.  Run\n  fit-model.R -h\n for help.\n") }
if (!file.exists(opt$infile)) { stop("Cannot read input file.") }
if ((!is.null(opt$configfile)) && (!file.exists(opt$configfile)) ) { stop("Could not find config file `", opt$configfile, "`.") }
if (is.null(opt$basedir)) { opt$basedir <- dirname(opt$infile) }
if (is.null(opt$outfile)) { opt$outfile <- paste( opt$basedir, "/", gsub("\\.[^.]*","",basename(opt$infile) ), "-", gsub("\\.[^.]*","",basename(opt$gmfile) ), "-", if (opt$nonnegative) { "nonnegative" } else {"unconstrained"}, "-", opt$jobid, ".RData", sep='' ) }
if ( !is.null(opt$seed) ) { set.seed(opt$seed) }
print(opt) # this will go in the pbs log

library(contextual)
library(contextutils)

options(error = print.and.dump)

# read configuration
init.config <- fill.default.config( read.config(opt$configfile) )  # returns NULL if it is NULL
if (!is.null(init.config$tree)) { warning("Appears to be a tree-config file; should be using fit-tree-model.R?") }
if ( is.null(init.config$mutrates) || 
        is.null(init.config$mutrates.scale) || 
        is.null(init.config$selcoef) || 
        is.null(init.config$selcoef.scale) || 
        is.null(init.config$fixfn.params) || 
        is.null(init.config$fixfn.params.scale) ) {
    stop("Need to specify starting parameter values and scalings in config file.")
}

# load generator matrix
stopifnot(file.exists(opt$gmfile))
load(opt$gmfile)  # provides 'genmatrix'
check.genmatrix(init.config,genmatrix)

# read in counts
counts <- read.counts(opt$infile, bases=genmatrix@bases, longpats=rownames(genmatrix) )
stopifnot( all( rownames(counts) == rownames(genmatrix) ) )

# shorter counts for full likelihood
stopifnot(shortwin(counts) > 1)
counts.0 <- projectcounts(counts, new.shortwin=shortwin(counts) - 1)

projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin(counts), fpatterns=colnames(counts) )
projmatrix.0 <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin(counts), fpatterns=colnames(counts.0) )

if (FALSE) {  ## DO NOT DO THIS; move elsewhere
    if (is.null(init.config$mutrates)) {
        # get ad hoc initial guesses at parameters
        adhoc <- countmuts(counts=counts,mutpats=genmatrix@mutpats,leftwin=leftwin(counts))
        init.config$mutrates <- adhoc[1,]/adhoc[2,]
        init.config$mutrates <- ifelse( is.finite(init.config$mutrates) & init.config$mutrates > 0, init.config$mutrates, 1e-4 )
    }
    if (is.null(init.config$selcoef)) {
        # start selection parameters at zero, why not?
        init.config$selcoef <- rep(1e-6,nsel(genmatrix))
        names(init.config$selcoef) <- selnames(genmatrix@selpats)
    }
    if (is.null(init.config$fixfn.params)) {
        # how to choose additional parameters?? Need guidance here.
        init.config$fixfn.params <- rep(1,length(fixparams(genmatrix)))
        names(init.config$fixfn.params) <- fixparams(genmatrix)
        init.config$fixfn.params[names(init.config$fixfn.params)=="Ne"] <- 200 + 2000*runif(1)
    }
    # scale tuning parameters
    if (is.null(init.config$mutrates.scale)) {
        init.config$mutrates.scale <- 1e-3 * rep( mean(init.config$mutrates),nmuts(genmatrix) )
    }
    if (is.null(init.config$selcoef.scale)) {
        init.config$selcoef.scale <- rep(.001,nsel(genmatrix))
    }
    if (is.null(init.config$fixfn.params.scale)) {
        init.config$fixfn.params.scale <- 0.1*init.config$fixfn.params
    }
}

initpar <- with(init.config, unlist(c(mutrates*opt$tlen,selcoef,fixfn.params)) )
parscale <- with(init.config, unlist( c(mutrates.scale*opt$tlen, selcoef.scale, fixfn.params.scale) ) )
names(initpar) <- names(parscale) <- c( mutnames(genmatrix@mutpats), selnames(genmatrix@selpats), fixparams(genmatrix) )

# skip these parameters
use.par <- ( parscale!=0 )
params <- initpar

# The log-likelihood function.
likfun <- function (sub.params){
    # params are: mutrates, selcoef, fixparams
    params[use.par] <- sub.params
    fparams <- params[seq( 1+nmuts(genmatrix)+nsel(genmatrix), 
                           length.out=length(fixparams(genmatrix)) )]
    names(fparams) <- fixparams(genmatrix)
    genmatrix@x <- do.call( update_x, 
                           c( list( 
                                   G=genmatrix,mutrates=params[1:nmuts(genmatrix)],
                                   selcoef=params[seq(1+nmuts(genmatrix),length.out=nsel(genmatrix))]), 
                              as.list(fparams) ) )
    # this is collapsed transition matrix
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, tlen=1, time="fixed") # shape=params[length(params)], time="gamma" )
    subtransmatrix.0 <- computetransmatrix( genmatrix, projmatrix.0, tlen=1, time="fixed") # shape=params[length(params)], time="gamma" )
    # return POSITIVE log-likelihood
    num <- sum( counts@counts * log(subtransmatrix) )
    den <- sum( counts.0@counts * log(subtransmatrix.0) )
    ans <- num - den
    if (!is.finite(ans)) { print(paste("Warning: non-finite likelihood with params:",paste(params,collapse=", "))) }
    return(ans)
}

stopifnot( length(initpar) == nmuts(genmatrix)+nsel(genmatrix)+length(fixparams(genmatrix)) )
mut.lb <- if (opt$nonnegative) { 1e-6 } else { -Inf }
lbs <- c( rep(mut.lb,nmuts(genmatrix)), rep(-5,nsel(genmatrix)), rep(-Inf,length(fixparams(genmatrix))) )
ubs <- c( rep(Inf,nmuts(genmatrix)), rep(5,nsel(genmatrix)), rep(Inf,length(fixparams(genmatrix))) )

baseval <- likfun(initpar[use.par])
if ( ! is.finite(baseval) ) { stop("Likelihood is not finite at initial values.") }
# # adjust parscale if necessary
# for (k in seq_along(initpar)) {
#     while ( parscale[k] > 1e-8 && (!is.finite(likfun(initpar+ifelse(seq_along(parscale)==k,parscale,0))) ) || (!is.finite(likfun(initpar+ifelse(seq_along(parscale)==k,parscale,0))) ) ) {
#         cat("Reducing parscale for ", names(initpar)[k], ".\n")
#         parscale[k] <- parscale[k]/4
#     }
# }

optim.results <- optim( par=initpar[use.par], fn=likfun, method="L-BFGS-B", lower=lbs[use.par], upper=ubs[use.par], control=list(fnscale=(-1)*abs(baseval), parscale=parscale[use.par], maxit=opt$maxit) )

# save some more things in optim.results
optim.results$use.par <- use.par
optim.results$parscale <- parscale
optim.results$initpar <- initpar

optim.par <- initpar
optim.par[use.par] <- optim.results$par
optim.results$par <- optim.par


model <- new( "context",
             counts=counts,
             genmatrix=genmatrix,
             projmatrix=projmatrix,
             mutrates=optim.results$par[1:nmuts(genmatrix)],
             selcoef=optim.results$par[seq(nmuts(genmatrix)+1,length.out=nsel(genmatrix))],
             params=optim.results$par[seq(1+nmuts(genmatrix)+nsel(genmatrix),length.out=length(fixparams(genmatrix)))],
             results=optim.results,
             likfun=likfun,
             invocation=invocation
         )

# set this up so that we can call likfun again in the future, directly
likfun.env <- new.env()
assign("use.par",use.par,envir=likfun.env)
assign("params",optim.results$par,envir=likfun.env)
assign("genmatrix",model@genmatrix,envir=likfun.env)
assign("projmatrix",model@projmatrix,envir=likfun.env)
assign("counts",model@counts,envir=likfun.env)
environment(model@likfun) <- likfun.env

save(model,file=opt$outfile)
