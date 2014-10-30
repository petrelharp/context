#!/usr/bin/Rscript --vanilla
require(optparse)


invocation <- commandArgs()

usage <- "\
Infer parameters from paired counts file, which records instances of Tmer transitions.
"

option_list <- list(
    # input/output
        make_option( c("-i","--infile"), type="character", help="Input file with tuple counts, tab-separated, with header 'reference', 'derived', 'count'. [default, looks in basedir]" ),
        make_option( c("-c","--configfile"), type="character", help="Config file with initial guesses at parameter values. [default: makes cheap guesses]"), 
        make_option( c("-o","--outfile"), type="character", help="File to save results to.  [default: base of infile + base of genmatrix + jobid + .RData]"),
        make_option( c("-u","--basedir"), type="character", default=NULL, help="Directory to put output in. [default: same as infile]"),
        make_option( c("-l","--leftwin"), type="integer", help="Size of left-hand context." ),
        make_option( c("-m","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix."),
        make_option( c("-x","--maxit"), type="integer", default=100, help="Number of iterations of optimization to run for. [default=%default]"),
        make_option( c("-j","--jobid"), type="character", default=formatC(1e6*runif(1),width=6,format="d",flag="0"), help="Unique job id. [default random]")
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$infile)) { stop("No input file.  Run\n  bcells-inference.R -h\n for help.\n") }
if (!file.exists(opt$infile)) { stop("Cannot read input file.") }
if ((!is.null(opt$configfile)) && (!file.exists(opt$configfile)) ) { stop("Could not find config file `", opt$configfile, "`.") }
if (is.null(opt$basedir)) { opt$basedir <- dirname(opt$infile) }
if (is.null(opt$outfile)) { opt$outfile <- paste( opt$basedir, "/", gsub("\\.[^.]*","",basename(opt$infile) ), "-", gsub("\\.[^.]*","",basename(opt$gmfile) ), "-", opt$jobid, ".RData", sep='' ) }
print(opt) # this will go in the pbs log

source("../context-inference-fns.R")
options(error = print.and.dump)

# load generator matrix
stopifnot(file.exists(opt$gmfile))
load(opt$gmfile)  # provides 'genmatrix'

# read in counts
counts <- read.counts(opt$infile, leftwin=opt$leftwin, bases=genmatrix@bases, longpats=rownames(genmatrix) )
stopifnot( all( rownames(counts) == rownames(genmatrix) ) )

projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin(counts), fpatterns=colnames(counts) )

if (!is.null(opt$configfile)) {
    init.config <- read.config(opt$configfile)
    initpar <- with( init.config, c( mutrates, selcoef, unlist(fixfn.params) ) )
    parscale <- c( 1e-3 * rep( mean(init.config$mutrates),nmuts(genmatrix)), rep(max(1e-4,mean(abs(init.config$selcoef))),nsel(genmatrix)), 0.1*unlist(init.config$fixfn.params) )
    names(initpar) <- c( mutnames(genmatrix@mutpats), selnames(genmatrix@selpats), fixparams(genmatrix) )
} else {
    # get ad hoc initial guesses at parameters
    adhoc <- countmuts(counts=counts,mutpats=genmatrix@mutpats,leftwin=leftwin(counts))
    adhoc.mutrates <- adhoc[1,]/adhoc[2,]
    adhoc.mutrates <- ifelse( is.finite(adhoc.mutrates) & adhoc.mutrates > 0, adhoc.mutrates, 1e-4 )
    # start selection parameters at zero, why not?
    adhoc.selcoef <- rep(1e-6,nsel(genmatrix))
    names(adhoc.selcoef) <- selnames(genmatrix@selpats)
    # how to choose additional parameters?? Need guidance here.
    adhoc.fixparams <- rep(1,length(fixparams(genmatrix)))
    names(adhoc.fixparams) <- fixparams(genmatrix)
    adhoc.fixparams[names(adhoc.fixparams)=="Ne"] <- 200 + 2000*runif(1)
    initpar <- c( adhoc.mutrates, adhoc.selcoef, adhoc.fixparams )
    parscale <- c( 1e-3 * rep( mean(adhoc.mutrates),nmuts(genmatrix)), rep(.05,nsel(genmatrix)), 0.1*adhoc.fixparams )
    names(initpar) <- c( mutnames(genmatrix@mutpats), selnames(genmatrix@selpats), fixparams(genmatrix) )
}

# Compute (quasi)-likelihood function using all counts -- multinomial as described in eqn:comp_like.
likfun <- function (params){
    # params are: mutrates, selcoef, fixparams
    fparams <- params[seq( 1+nmuts(genmatrix)+nsel(genmatrix), length.out=length(fixparams(genmatrix)) )]
    names(fparams) <- fixparams(genmatrix)
    genmatrix@x <- do.call( update, c( list( G=genmatrix,mutrates=params[1:nmuts(genmatrix)],selcoef=params[seq(1+nmuts(genmatrix),length.out=nsel(genmatrix))]), as.list(fparams) ) )
    # this is collapsed transition matrix
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, tlen=1, time="fixed") # shape=params[length(params)], time="gamma" )
    # return POSITIVE log-likelihood
    ans <- sum( counts@counts * log(subtransmatrix) )
    if (!is.finite(ans)) { print(paste("Warning: non-finite likelihood with params:",paste(params,collapse=", "))) }
    return(ans)
}

stopifnot( length(initpar) == nmuts(genmatrix)+nsel(genmatrix)+length(fixparams(genmatrix)) )
lbs <- c( rep(1e-6,nmuts(genmatrix)), rep(-5,nsel(genmatrix)), rep(-Inf,length(fixparams(genmatrix))) )
ubs <- c( rep(2,nmuts(genmatrix)), rep(5,nsel(genmatrix)), rep(Inf,length(fixparams(genmatrix))) )

baseval <- likfun(initpar)
stopifnot( is.finite(baseval) )
# adjust parscale if necessary
for (k in seq_along(initpar)) {
    while ( parscale[k] > 1e-8 && (!is.finite(likfun(initpar+ifelse(seq_along(parscale)==k,parscale,0))) ) || (!is.finite(likfun(initpar+ifelse(seq_along(parscale)==k,parscale,0))) ) ) {
        cat("Reducing parscale for ", names(initpar)[k], ".\n")
        parscale[k] <- parscale[k]/4
    }
}

optim.results <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(fnscale=(-1)*abs(baseval), parscale=parscale, maxit=opt$maxit) )


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

save(model,file=opt$outfile)
