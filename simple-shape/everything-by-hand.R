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

source("../context-inference-fns.R",chdir=TRUE)
source("../sim-context-fns.R",chdir=TRUE)


model.config <- read.config("simple-shape-model.json")

config <- treeify.config(model.config,tlen=10)
config <- parse.models(config)
simseq.tree(200,config,count.trans=TRUE)

## sim-seq
config <- treeify.config(model.config,tlen=1)
config <- parse.models(config)
simseqs <- simseq.tree(1e3,config,count.trans=TRUE)


## count-seq
opt <- list( longwin=8, shortwin=3, leftwin=2)
longpats <- getpatterns(opt$longwin,config$bases)
shortpats <- getpatterns(opt$shortwin,config$bases)
counts <- counttrans.list( list(longpats,shortpats), 
    simseqs=simseqs, 
    leftwin=opt$leftwin, bases=config$bases,
   shift=0 )

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
genmatrix@x <- do.call( update, c( list(G=genmatrix, mutrates=config[[mm]]$mutrates, selcoef=config[[mm]]$selcoef), config[[mm]]$fixfn.params ) )


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
    genmatrix@x <- do.call( update, c( list( G=genmatrix,mutrates=params[1:nmuts(genmatrix)],selcoef=params[seq(1+nmuts(genmatrix),length.out=nsel(genmatrix))]), as.list(fparams) ) )
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

optim.results <- optim( par=initpar[use.par], fn=likfun, method="L-BFGS-B", lower=lbs[use.par], upper=ubs[use.par], control=list(fnscale=(-1)*abs(baseval), parscale=parscale[use.par] ) )

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
    genmatrixfile     = opt$genmatrixfile)

