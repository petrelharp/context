###
# EVERYTHING
#
# i.e. works through the important bits, by hand, in R
# of the various scripts, namely
#   read in config
#   simulate sequence
#   count tuples
#   fit model
#   compute residuals

library(contextual)
library(simcontext)

library(MASS)
library(numDeriv)

modelfile <- "less-simple-shape-model.json"

model.config <- read.config(modelfile)

config <- treeify.config(model.config,tlen=10)
config <- parse.models(config)
simseq.tree(200,config,count.trans=TRUE)

## sim-seq
config <- treeify.config(model.config,tlen=1)
config <- parse.models(config)
simseqs <- simseq.tree(1e3,config,count.trans=TRUE)


## options
opt <- list( longwin=8, shortwin=3, leftwin=2)

## count-seq
longpats <- getpatterns(opt$longwin,config$bases)
shortpats <- getpatterns(opt$shortwin,config$bases)
counts <- counttrans.list( list(longpats,shortpats), 
    simseqs=simseqs, 
    leftwin=opt$leftwin, bases=config$bases,
   shift=0 )

## make-genmatrix
mm <- "tip"
opt$boundary <- "none"
genmatrix <- makegenmatrix(
                patlen=opt$longwin, 
                mutpats=config[[mm]]$mutpats, 
                selpats=config[[mm]]$selpats, 
                selfactors=config[[mm]]$selfactors, 
                boundary=opt$boundary, 
                bases=config[[mm]]$bases, 
                fixfn=config[[mm]]$fixfn,
                Ne=config[[mm]]$fixfn.params$Ne
            )
genmatrix@x <- do.call( update_x, c( list(G=genmatrix, mutrates=config[[mm]]$mutrates, selcoef=config[[mm]]$selcoef), config[[mm]]$fixfn.params ) )


## fit-model
opt$tlen <- 1

projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin(counts), fpatterns=colnames(counts) )

initpar <- with(model.config, unlist(c(mutrates*opt$tlen,selcoef,fixfn.params)) )
parscale <- with(model.config, unlist( c(mutrates.scale*opt$tlen, selcoef.scale, fixfn.params.scale) ) )
names(initpar) <- names(parscale) <- c( mutnames(genmatrix@mutpats), selnames(genmatrix@selpats), fixparams(genmatrix) )

# skip these parameters
use.par <- ( parscale!=0 )
params <- initpar

likfun <- function (sub.params){
    # params are: mutrates, selcoef, fixparams
    params[use.par] <- sub.params
    fparams <- params[seq( 1+nmuts(genmatrix)+nsel(genmatrix), length.out=length(fixparams(genmatrix)) )]
    names(fparams) <- fixparams(genmatrix)
    genmatrix@x <- do.call( update_x, c( list( G=genmatrix,mutrates=params[1:nmuts(genmatrix)],selcoef=params[seq(1+nmuts(genmatrix),length.out=nsel(genmatrix))]), as.list(fparams) ) )
    # this is collapsed transition matrix
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, tlen=1, time="fixed") # shape=params[length(params)], time="gamma" )
    # return POSITIVE log-likelihood
    ans <- sum( counts@counts * log(subtransmatrix) )
    if (!is.finite(ans)) { print(paste("Warning: non-finite likelihood with params:",paste(params,collapse=", "))) }
    return(ans)
}

stopifnot( length(initpar) == nmuts(genmatrix)+nsel(genmatrix)+length(fixparams(genmatrix)) )
lbs <- c( rep(1e-8,nmuts(genmatrix)), rep(-5,nsel(genmatrix)), rep(-Inf,length(fixparams(genmatrix))) )
ubs <- c( rep(2,nmuts(genmatrix)), rep(5,nsel(genmatrix)), rep(Inf,length(fixparams(genmatrix))) )

baseval <- likfun(initpar[use.par])
stopifnot( is.finite(baseval) )

optim.results <- optim( par=initpar[use.par], fn=likfun, method="L-BFGS-B", lower=lbs[use.par], upper=ubs[use.par], control=list(fnscale=(-1)*abs(baseval), parscale=parscale[use.par] ), hessian=TRUE )

# save some more things in optim.results
optim.results$use.par <- use.par
optim.results$parscale <- parscale
optim.results$initpar <- initpar

optim.par <- initpar
optim.par[use.par] <- optim.results$par
optim.results$par <- optim.par

fit.model <- new( "context",
             counts=counts,
             genmatrix=genmatrix,
             projmatrix=projmatrix,
             mutrates=optim.results$par[1:nmuts(genmatrix)],
             selcoef=optim.results$par[seq(nmuts(genmatrix)+1,length.out=nsel(genmatrix))],
             params=optim.results$par[seq(1+nmuts(genmatrix)+nsel(genmatrix),length.out=length(fixparams(genmatrix)))],
             results=optim.results,
             likfun=likfun,
             invocation="by hand"
         )



## compute-resids


residframe <- computeresids (fit.model,
        pretty            = TRUE,
        in_longwin        = opt$longwin,
        in_shortwin       = opt$shortwin,
        in_leftwin        = opt$leftwin,
        counts            = counts,
        genmatrixfile     = opt$genmatrixfile
    )

## look at likelihood surface

# compare likelihood profiles
time <- as.numeric(opt$tlen)
timevec <- c( rep(as.numeric(opt$tlen),nmuts(fit.model)), rep(1,length(coef(fit.model))-nmuts(fit.model)) )
true.model <- fit.model
true.model@mutrates <- config$tip$mutrates*time
true.model@selcoef <- config$tip$selcoef

compare.params <- data.frame( 
        fit=coef(fit.model)/timevec,
        simulated=coef(true.model)/timevec,
        fixed=!fit.model@results$use.par
    )

# Hessian of the likelihood surface
likfun.hess <- hessian(function(x)fit.model@likfun(x*timevec[fit.model@results$use.par]),method.args=list(eps=min(coef(model)/10)),x=coef(model)[fit.model@results$use.par])
colnames( likfun.hess ) <- rownames( likfun.hess ) <- names(coef(fit.model))[fit.model@results$use.par]
# asymptotic covariance matrix
est.covmat <- ginv(likfun.hess)

f <- function (x) { fit.model@likfun((x*timevec)[fit.model@results$use.par]) }

# profiles
for (k in (1:nrow(compare.params))[fit.model@results$use.par]) {
    pvals <- seq( .9*min(compare.params[k,c("fit","simulated")]), 1.1*max(compare.params[k,c("fit","simulated")]), length.out=10 )
    inds <- 1:nrow(compare.params)
    inds[k] <- 1+nrow(compare.params)
    lvals <- sapply( pvals, function (x) { f( c(compare.params$fit,x)[inds] ) } )
    l2vals <- sapply( pvals, function (x) { f( c(compare.params$simulated,x)[inds] ) } )
    plot( pvals, lvals, type='b', main=paste("likelihood surface for", compare.params$name[k]), ylab="log-likelihood", xlab='value', ylim=range(lvals,l2vals,finite=TRUE) )
    lines( pvals, l2vals, col='red', type='b' )
    abline(v=compare.params[k,c("fit","simulated")], col=c("black","red"))
}

# contours
for (j in (1:nrow(compare.params))[fit.model@results$use.par]) {
    for (k in (1:nrow(compare.params))[fit.model@results$use.par]) {
        inds <- 1:nrow(compare.params)
        inds[j] <- 1+nrow(compare.params)
        inds[k] <- 2+nrow(compare.params)
        pvals <- seq( .9*min(compare.params[j,c("fit","simulated")]), 1.1*max(compare.params[j,c("fit","simulated")]), length.out=8 )
        qvals <- seq( .9*min(compare.params[k,c("fit","simulated")]), 1.1*max(compare.params[k,c("fit","simulated")]), length.out=8 )
        pqvals <- expand.grid( pvals, qvals )
        lvals <- apply( pqvals, 1, function (x) { f( c(compare.params$fit,x)[inds] ) } )
        cols <- colorspace::diverge_hcl(64)
        plot( pqvals, cex=3, pch=20, col=cols[cut(lvals,breaks=length(cols))] )
    }
}

