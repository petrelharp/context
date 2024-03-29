#!/usr/bin/env Rscript
library(optparse)

invocation <- commandArgs()

usage <- "\
Infer parameters from paired counts file, which records instances of Tmer transitions.
\
See model-desc.md for a description of the config file.
\
Parameters whose scale is set to zero *will be regarded as fixed.* \
"

option_list <- list(
    # input/output
        make_option( c("-i","--infile"), type="character", help="Input file with tuple counts, tab-separated, with header 'reference', 'derived', 'count'. [default, looks in basedir]" ),
        make_option( c("-c","--configfile"), type="character", help="File with model configuration (optional if genmatrix provided, not on a tree)."),
        make_option( c("-o","--outfile"), type="character", help="File to save results to.  [default: base of infile + base of genmatrix + jobid + .RData]"),
        make_option( c("-u","--basedir"), type="character", default=NULL, help="Directory to put output in. [default: same as infile]"),
        make_option( c("-x","--maxit"), type="integer", default=100, help="Number of iterations of optimization to run for. [default=%default]"),
        make_option( c("-n","--unconstrained"), action="store_true", default=FALSE, help="Do not constrain mutation rates to be nonnegative? [default=%default]"),
        make_option( c("-j","--jobid"), type="character", default=formatC(1e6*runif(1),width=6,format="d",flag="0"), help="Unique job id. [default random]"),
        make_option( c("-z","--seed"), type="integer", help="Seed for pseudorandom number generator; an integer. [default: does not meddle]")
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
opt$nonnegative <- !opt$unconstrained
if (is.null(opt$infile)) { stop("No input file.  Run\n  fit-tree-model.R -h\n for help.\n") }
if (is.null(opt$configfile)) { stop("No config file.  Run\n  fit-tree-model.R -h\n for help.\n") }
if ((!is.null(opt$configfile)) && (!file.exists(opt$configfile)) ) { stop("Could not find config file `", opt$configfile, "`.") }
if (!file.exists(opt$infile)) { stop("Cannot read input file.") }
if (is.null(opt$basedir)) { opt$basedir <- dirname(opt$infile) }
if (is.null(opt$outfile)) { 
    model.id <- if (is.null(opt$configfile)) { gsub("\\.[^.]*","",basename(opt$gmfile) ) } else { gsub("\\.[^.]*","",basename(opt$configfile) ) }
    opt$outfile <- paste( opt$basedir, "/", gsub("\\.[^.]*","",basename(opt$infile) ), "-", model.id, "-", if (opt$nonnegative) { "nonnegative" } else {"unconstrained"}, "-", opt$jobid, ".RData", sep='' ) 
}
print(opt) # this will go in the pbs log

library(contextual)
library(contextutils)

options(error = print.and.dump)

# read in config
config <- parse.models( treeify.config( read.config(opt$configfile) ) )
counts.config <- read.config.counts(opt$infile)

# find the right generator matrix files
for (mm in config$.models) { 
    config[[mm]]$genmatrix <- file.path(dirname(opt$configfile), gsub("%",counts.config$longwin,config[[mm]]$genmatrix,fixed=TRUE)) 
}


# which models go with which edges
modelnames <- config.dereference( config, nodenames(config$tree) )
# load generator matrices
genmatrices <- lapply( selfname(config$.models), function (mm) {
        if (!file.exists(config[[mm]]$genmatrix)) { stop(paste("Can't find file", config[[mm]]$genmatrix), ".") }
        load( config[[mm]]$genmatrix )
        check.genmatrix( config[[mm]], genmatrix )
        return(genmatrix)
    } )
longpats <- rownames(genmatrices[[1]])
stopifnot( all( sapply( genmatrices, function(gm) { all(rownames(gm)==longpats) } ) ) )

# read in counts
counts <- read.counts(opt$infile, bases=config$bases, longpats=rownames(genmatrices[[1]]) )
stopifnot( all( rownames(counts) == longpats ) )
counts.0 <- projectcounts(counts, new.shortwin=shortwin(counts) - 1)

projmatrix <- collapsepatmatrix( ipatterns=longpats, leftwin=leftwin(counts), 
                                shortwin=shortwin(counts), bases=counts@bases )
projmatrix.0 <- collapsepatmatrix( ipatterns=longpats, leftwin=leftwin(counts), 
                                  shortwin=shortwin(counts.0), bases=counts@bases )

# parameters are ordered as: initfreqs, tlens, ( mutrates, selcoef, fixparams ) x genmatrices
# this will parse a vector of params and give the ones requested
.nparams <- sapply( genmatrices, function (x) { ( nmuts(x) + nsel(x) + length(fixparams(x)) ) } )
.param.info <- lapply( seq_along(genmatrices), function (k) {
        nbases <- length(genmatrices[[k]]@bases)
        ntimes <- nrow(config$tree$edge)
        nothers <- if (k>1) { cumsum(.nparams)[k-1] } else { 0 }
        nmuts <- nmuts(genmatrices[[k]])
        nsel <- nsel(genmatrices[[k]])
        nfixparams <- length(fixparams(genmatrices[[k]]))
        list( 
            "initfreqs" = 1:nbases,
            "tlens" = nbases + (1:ntimes),
            "mutrates" = nbases + ntimes + nothers + seq(1,length.out=nmuts),
            "selcoef" = nbases + ntimes + nothers + nmuts + seq(1,length.out=nsel),
            "fixparams" = nbases + ntimes + nothers + nmuts + nsel + seq(1,length.out=nfixparams)
        ) } )
names(.param.info) <- names(genmatrices)
.param.map <- function (gm=1,type,params) { params[.param.info[[gm]][[type]]] }

initparam.list <- c( list( config$initfreqs, config$tree$edge.length ), lapply( config[config$.models], function (x) {
                    y <- unlist( x[ c("mutrates","selcoef","fixfn.params") ] )
                    names(y) <- c( mutnames(x$mutpats), selnames(x$selpats), names(x$fixfn.params) )
                    return(y)
    } ) )
names(initparam.list[[1]]) <- paste("init",config$bases,sep='.')
names(initparam.list[[2]]) <- paste("tlen",nodenames(config$tree)[config$tree$edge[,2]],sep='.')
initparams <- unlist( initparam.list )

parscale.list <- c( list( config$initfreqs.scale, config$tlen.scale ), lapply( config[config$.models], function (x) {
                    y <- unlist( x[ c("mutrates.scale","selcoef.scale","fixfn.params.scale") ] )
                    names(y) <- c( mutnames(x$mutpats), selnames(x$selpats), names(x$fixfn.params) )
                    return(y)
    } ) )
names(parscale.list[[1]]) <- paste("init",config$bases,sep='.')
names(parscale.list[[2]]) <- paste("tlen",nodenames(config$tree)[config$tree$edge[,2]],sep='.')
parscale <- unlist( parscale.list )

# skip these parameters
use.par <- ( parscale!=0 )
params <- initparams

# compute root distribution
initfreq.index  <- product.index( longpats=longpats, bases=config$bases ) # which base is at each position in each pattern
initfreqs <- .param.map( type="initfreqs", params=initparams )
initfreqs <- initfreqs/sum(initfreqs)
root.distrn <- get.root.distrn( initfreqs, initfreq.index )

# create setup to efficiently re-compute below
# note that normalize=FALSE can be used -- the answer is the same whether normalize=TRUE or FALSE
#  and this skips the row normalization step.
peel.setup <- peel.transmat( tree=config$tree, rowtaxon=rowtaxon(counts), coltaxa=coltaxa(counts), modelnames=modelnames, genmatrices=genmatrices, 
                            projmatrix=projmatrix, root.distrn=root.distrn, tlens=config$tree$edge.length, return.list=TRUE, normalize=FALSE )
peel.setup.0 <- peel.transmat( tree=config$tree, rowtaxon=rowtaxon(counts), coltaxa=coltaxa(counts), modelnames=modelnames, genmatrices=genmatrices, 
                            projmatrix=projmatrix.0, root.distrn=root.distrn, tlens=config$tree$edge.length, return.list=TRUE, normalize=FALSE )

# reorder columns of counts to match the order we compute transmat in
counts <- reorder.counts( counts, nodenames(config$tree)[peel.setup$col.order[[peel.setup$row.node]]] )
counts.0 <- reorder.counts( counts.0, nodenames(config$tree)[peel.setup.0$col.order[[peel.setup.0$row.node]]] )

# Compute likelihood function using all counts
likfun <- function (sub.params){
    params[use.par] <- sub.params
    # params are: initfreqs, tlens, ( mutrates, selcoef, fixparams ) x genmatrices
    # First, update genmatrices:
    mutrates.list <- lapply( seq_along(genmatrices), .param.map, type="mutrates", params=params )
    selcoef.list <- lapply( seq_along(genmatrices), .param.map, type="selcoef", params=params )
    fixparam.list <- lapply( seq_along(genmatrices), function (k) { x <- .param.map(k, type="fixparams", params=params); names(x) <- fixparams(genmatrices[[k]]); x } )
    for (k in seq_along(genmatrices)) {
        genmatrices[[k]]@x <- do.call( update_x, c( list( G=genmatrices[[k]],mutrates=mutrates.list[[k]],selcoef=selcoef.list[[k]] ), as.list(fixparam.list[[k]]) ) )
    }
    # get branch lengths
    tlens <- .param.map( type="tlens", params=params )
    # compute root distribution
    initfreqs <- .param.map( type="initfreqs", params=params )
    initfreqs <- initfreqs/sum(initfreqs)
    root.distrn <- get.root.distrn( initfreqs, initfreq.index )
    # Now, peel 
    transmat <- peel.transmat.compute( setup=peel.setup, genmatrices=genmatrices, root.distrn=root.distrn, tlens=tlens, return.list=FALSE, normalize=FALSE )
    transmat.0 <- peel.transmat.compute( setup=peel.setup.0, genmatrices=genmatrices, root.distrn=root.distrn, tlens=tlens, return.list=FALSE, normalize=FALSE )
    # return POSITIVE log-likelihood
    num <- sum( counts@counts * log(transmat) )
    den <- sum( counts.0@counts * log(transmat.0) )
    ans <- num-den
    if (!is.finite(ans)) { print(paste("Warning: non-finite likelihood with params:",params)) }
    return(ans)
}

mut.lb <- if (opt$nonnegative) { 1e-6 } else { -Inf }
lbs <- unlist( c( 
                rep(1e-6, length(config$bases)), 
                rep(1e-6, length(config$tree$edge.length)), 
                lapply( genmatrices, function (gm) { c( 
                                rep(mut.lb, nmuts(gm)), 
                                rep(-5, nsel(gm)), 
                                rep(-Inf, length(fixparams(gm))) ) } ) ) )
ubs <- unlist( c( 
                rep(1, length(config$bases)), 
                rep(Inf, length(config$tree$edge.length)), 
                lapply( genmatrices, function (gm) { c( 
                                rep(Inf, nmuts(gm)), 
                                rep(5, nsel(gm)), 
                                rep(Inf, length(fixparams(gm))) ) } ) ) )

stopifnot( all( initparams >= lbs ) && all( initparams <= ubs ) )

likfun.time <- system.time( { baseval <- likfun(initparams[use.par]) } )
cat("Time to evaluate likelihood:\n")
print(likfun.time)
stopifnot( is.finite(baseval) )

optim.results <- optim( par=initparams[use.par], fn=likfun, method="L-BFGS-B", lower=lbs[use.par], upper=ubs[use.par], control=list(fnscale=(-1)*abs(baseval), parscale=parscale[use.par], maxit=opt$maxit) )

# put back in parameters not optimized over into results
optim.par <- initparams
optim.par[use.par] <- optim.results$par
optim.results$par <- optim.par

# and into the tree
config$tree$edge.length <- .param.map(type="tlens",params=optim.par)

model <- new( "contextTree",
             counts=counts,
             tree=config$tree,
             initfreqs=.param.map(type="initfreqs",params=optim.results$par),
             models=lapply( config$.models, function (mname) {
                         # remove prepended model names from names
                         nn <- function (x) { names(x) <- gsub(paste0("^",mname,"."),"",names(x)); x }
                         new( "contextModel",
                             genmatrix=genmatrices[[mname]],
                             projmatrix=projmatrix,
                             mutrates=nn(.param.map(mname,"mutrates",optim.results$par)),
                             selcoef=nn(.param.map(mname,"selcoef",optim.results$par)),
                             params=nn(.param.map(mname,"fixparams",optim.results$par))
                         )
                     } ),
             modelnames=modelnames,
             results=optim.results,
             likfun=likfun,
             invocation=invocation
         )

# set this up so that we can call likfun again in the future, directly
# ... except it DOESN'T WORK; do it better
likfun.env <- new.env()
assign("use.par",use.par,envir=likfun.env)
assign("params",params,envir=likfun.env)
assign(".param.map",.param.map,envir=likfun.env)
assign(".param.info",.param.info,envir=likfun.env)
assign("peel.setup",peel.setup,envir=likfun.env)
assign("peel.setup.0",peel.setup.0,envir=likfun.env)
assign("modelnames",modelnames,envir=likfun.env)
assign("genmatrices",lapply(model@models,slot,"genmatrix"),envir=likfun.env)
assign("initfreq.index",initfreq.index,envir=likfun.env)
assign("counts",model@counts,envir=likfun.env)
assign("counts.0",counts.0,envir=likfun.env)
environment(model@likfun) <- likfun.env


cat("Saving output to ", opt$outfile, " .\n")
save(model,file=opt$outfile)
