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

genmatrices <- lapply(model@models, slot, "genmatrix")
counts <- model@counts
longpats <- rownames(counts)

# shorter counts for full likelihood
stopifnot(shortwin(counts) > 1)
counts.0 <- projectcounts(counts, new.shortwin = shortwin(counts) - 1)
projmatrix <- collapsepatmatrix( ipatterns=longpats, leftwin=leftwin(counts), 
                                shortwin=shortwin(counts), bases=counts@bases )
projmatrix.0 <- collapsepatmatrix( ipatterns=longpats, leftwin=leftwin(counts), 
                                  shortwin=shortwin(counts.0), bases=counts@bases )


# read in config file
prior.config <- parse.models( treeify.config( read.config(opt$configfile) ) )

for (model_name in names(model@models)) {
    if (! (model_name %in% names(prior.config))) {
        stop(model_name, "not in config.")
    }
    if ( length(prior.config[[model_name]]$mutrates) != length(model@models[[model_name]]@mutrates) 
         || length(prior.config[[model_name]]$selcoef) != length(model@models[[model_name]]@selcoef) 
         || length(prior.config[[model_name]]$fixfn.params) != length(model@models[[model_name]]@params) ) {
        stop(sprintf("Configuration in %s does not match that in already-fit model of %s .", opt$configfile, opt$infile))
    }
    if ( length(prior.config[[model_name]]$mutrates.prior.mean) != length(model@models[[model_name]]@mutrates) 
         || length(prior.config[[model_name]]$mutrates.prior.sd) != length(model@models[[model_name]]@mutrates) 
         || length(prior.config[[model_name]]$selcoef.prior.mean) != length(model@models[[model_name]]@selcoef) 
         || length(prior.config[[model_name]]$selcoef.prior.sd) != length(model@models[[model_name]]@selcoef) 
         || length(prior.config[[model_name]]$fixfn.params.prior.mean) != length(model@models[[model_name]]@params)
         || length(prior.config[[model_name]]$fixfn.params.prior.sd) != length(model@models[[model_name]]@params) ) {
        stop(sprintf("Configuration in %s does not have prior mean and SDs specified.", opt$configfile))
    }
}


# hackily,
# parameters are ordered as: initfreqs, tlens, ( mutrates, selcoef, fixparams ) x genmatrices
# this will parse a vector of params and give the ones requested
.nparams <- sapply( genmatrices, function (x) { ( nmuts(x) + nsel(x) + length(fixparams(x)) ) } )
.param.info <- lapply( seq_along(genmatrices), function (k) {
        nbases <- length(genmatrices[[k]]@bases)
        ntimes <- nrow(model@tree$edge)
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

initparam.list <- c(list( model@initfreqs, model@tree$edge.length ), 
                     lapply( model@models, function (x) {
                            c(x@mutrates, x@selcoef, x@params)
                    } ) )
initparams <- unlist( initparam.list )

parscale.list <- c(list( prior.config$initfreqs.scale, prior.config$tlen.scale ), 
                   lapply( names(model@models), function (model_name) {
                         x <- prior.config[[model_name]]
                         y <- unlist( x[ c("mutrates.scale","selcoef.scale","fixfn.params.scale") ] )
                         names(y) <- names(initparam.list[[model_name]])
                         return(y)
                } ) )
names(parscale.list[[1]]) <- paste("init",prior.config$bases,sep='.')
names(parscale.list[[2]]) <- paste("tlen",nodenames(prior.config$tree)[prior.config$tree$edge[,2]],sep='.')
parscale <- unlist( parscale.list )
parscale <- parscale * opt$scalefac


prior.mean.list <- c(list( prior.config$initfreqs.prior.mean, prior.config$tlen.prior.mean ), 
                     lapply( names(model@models), function (model_name) {
                         x <- prior.config[[model_name]]
                         y <- unlist( x[ c("mutrates.prior.mean","selcoef.prior.mean","fixfn.params.prior.mean") ] )
                         names(y) <- names(initparam.list[[model_name]])
                         return(y)
                    } ) )
names(prior.mean.list[[1]]) <- paste("init",prior.config$bases,sep='.')
names(prior.mean.list[[2]]) <- paste("tlen",nodenames(prior.config$tree)[prior.config$tree$edge[,2]],sep='.')
prior.mean <- unlist( prior.mean.list )

prior.sd.list <- c(list( prior.config$initfreqs.prior.sd, prior.config$tlen.prior.sd ), 
                     lapply( names(model@models), function (model_name) {
                         x <- prior.config[[model_name]]
                         y <- unlist( x[ c("mutrates.prior.sd","selcoef.prior.sd","fixfn.params.prior.sd") ] )
                         names(y) <- names(initparam.list[[model_name]])
                         return(y)
                    } ) )
names(prior.sd.list[[1]]) <- paste("init",prior.config$bases,sep='.')
names(prior.sd.list[[2]]) <- paste("tlen",nodenames(prior.config$tree)[prior.config$tree$edge[,2]],sep='.')
prior.sd <- unlist( prior.sd.list )


# skip these parameters
use.par <- ( parscale!=0 )
params <- initparams

# get prior parameters (all half-Gaussian)

# compute root distribution
initfreq.index  <- product.index( longpats=longpats, bases=prior.config$bases ) # which base is at each position in each pattern
initfreqs <- .param.map( type="initfreqs", params=initparams )
initfreqs <- initfreqs/sum(initfreqs)
root.distrn <- get.root.distrn( initfreqs, initfreq.index )

# create setup to efficiently re-compute below
# note that normalize=FALSE can be used -- the answer is the same whether normalize=TRUE or FALSE
#  and this skips the row normalization step.
peel.setup <- peel.transmat( tree=model@tree, rowtaxon=rowtaxon(counts), coltaxa=coltaxa(counts), modelnames=model@modelnames, genmatrices=genmatrices, 
                            projmatrix=projmatrix, root.distrn=root.distrn, tlens=model@tree$edge.length, return.list=TRUE, normalize=FALSE )
peel.setup.0 <- peel.transmat( tree=model@tree, rowtaxon=rowtaxon(counts), coltaxa=coltaxa(counts), modelnames=model@modelnames, genmatrices=genmatrices, 
                            projmatrix=projmatrix.0, root.distrn=root.distrn, tlens=model@tree$edge.length, return.list=TRUE, normalize=FALSE )

# reorder columns of counts to match the order we compute transmat in
counts <- reorder.counts( counts, nodenames(model@tree)[peel.setup$col.order[[peel.setup$row.node]]] )
counts.0 <- reorder.counts( counts.0, nodenames(model@tree)[peel.setup.0$col.order[[peel.setup.0$row.node]]] )


# The log-posterior function.
likfun <- function (sub.params){
    # params are: mutrates*tlen
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
    if (any(transmat<0) || any(transmat.0<0)) { return( -Inf ) }
    # return POSITIVE log-likelihood
    num <- sum( counts@counts * log(transmat) )
    den <- sum( counts.0@counts * log(transmat.0) )
    ans <- num-den
    # priors: all (half-)Gaussian
    if (!is.finite(ans)) { return( -Inf ) }
    else { return( ans - sum(((params[use.par]-prior.mean[use.par])/prior.sd[use.par])^2) ) }
}

likfun.time <- system.time( { baseval <- likfun(initparams[use.par]) } )
cat("Time to evaluate likelihood:\n")
print(likfun.time)
stopifnot( is.finite(baseval) )

mrun <- metrop( likfun, initial=initparams[use.par], nbatch=opt$nbatches, 
                blen=opt$blen, scale=parscale[use.par] )

mrun$use.par <- use.par
mrun$parscale <- parscale
mrun$initpar <- initparams
mrun.final.par <- params
mrun.final.par[use.par] <- mrun$final

model@tree$edge.length <- .param.map(type="tlens",params=mrun.final.par)
model@initfreqs <- .param.map(type="initfreqs",params=mrun.final.par)
final.mutrates.list <- lapply( seq_along(model@models), .param.map, type="mutrates", params=mrun.final.par )
final.selcoef.list <- lapply( seq_along(model@models), .param.map, type="selcoef", params=mrun.final.par )
final.fixparam.list <- lapply( seq_along(model@models), function (k) { x <- .param.map(k, type="fixparams", params=mrun.final.par); names(x) <- fixparams(model@models[[k]]@genmatrix); x } )
# remove prepended model names from names
nn <- function (x) { names(x) <- gsub(paste0("^",mname,"."),"",names(x)); x }
for (k in seq_along(model@models)) {
    mname <- names(model@models)[k]
    model@models[[k]]@mutrates=nn(.param.map(mname,"mutrates",mrun.final.par))
    model@models[[k]]@selcoef=nn(.param.map(mname,"selcoef",mrun.final.par))
    model@models[[k]]@params=nn(.param.map(mname,"fixparams",mrun.final.par))
    model@models[[k]]@genmatrix@x <- do.call( update_x, 
                                             c( list(G=model@models[[k]]@genmatrix, 
                                                     mutrates=final.mutrates.list[[k]],
                                                     selcoef=final.selcoef.list[[k]] ), 
                                               as.list(final.fixparam.list[[k]]) ) )
}
model@results <- unclass(mrun)
model@likfun <- likfun
model@invocation <- invocation

# Note that we are NOT saving a number of useful things, 
# for instance, the priors.


save(model, file=opt$outfile)

cat(sprintf("Done! Output to %s", opt$outfile))

