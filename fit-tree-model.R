#!/usr/bin/Rscript --vanilla
require(optparse)

invocation <- commandArgs()

usage <- "\
Infer parameters from paired counts file, which records instances of Tmer transitions.
"

option_list <- list(
    # input/output
        make_option( c("-i","--infile"), type="character", help="Input file with tuple counts, tab-separated, with header 'reference', 'derived', 'count'. [default, looks in basedir]" ),
        make_option( c("-l","--leftwin"), type="integer", help="Size of left-hand context for count file." ),
        make_option( c("-o","--outfile"), type="character", help="File to save results to.  [default: base of infile + base of genmatrix + jobid + .RData]"),
        make_option( c("-u","--basedir"), type="character", default=NULL, help="Directory to put output in. [default: same as infile]"),
        make_option( c("-c","--configfile"), type="character", help="File with model configuration (optional if genmatrix provided, not on a tree)."),
        make_option( c("-w","--longwin"), type="numeric", help="Size of long window. (optional if genmatrix specified, not on a tree)"),
        make_option( c("-s","--shortwin"), type="integer", help="Size of short window. (optional if genmatrix specified, not on a tree)" ),
        make_option( c("-x","--maxit"), type="integer", default=100, help="Number of iterations of optimization to run for. [default=%default]"),
        make_option( c("-j","--jobid"), type="character", default=formatC(1e6*runif(1),width=6,format="d",flag="0"), help="Unique job id. [default random]")
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$infile)) { stop("No input file.  Run\n  bcells-inference.R -h\n for help.\n") }
if (is.null(opt$config)) { stop("No config file.  Run\n  bcells-inference.R -h\n for help.\n") }
if (!file.exists(opt$infile)) { stop("Cannot read input file.") }
if (is.null(opt$basedir)) { opt$basedir <- dirname(opt$infile) }
if (is.null(opt$outfile)) { opt$outfile <- paste( opt$basedir, "/", gsub("\\.[^.]*","",basename(opt$infile) ), "-", gsub("\\.[^.]*","",basename(opt$gmfile) ), "-", opt$jobid, ".RData", sep='' ) }
print(opt) # this will go in the pbs log

source("../context-inference-fns.R")

options(error = print.and.dump)

config <- parse.models( treeify.config( read.config(opt$configfile) ) )
if (is.null(opt$longwin)) { stop("Must specify longwin and shortwin with tree-based models.") }
for (mm in config$.models) { config[[mm]]$genmatrix <- gsub("%",opt$longwin,config[[mm]]$genmatrix,fixed=TRUE) }

# load generator matrices
genmatrices <- lapply( selfname(config$.models), function (mm) {
        if (!file.exists(config[[mm]]$genmatrix)) { stop(paste("Can't find file", config[[mm]]$genmatrix), ".") }
        env <- new.env()
        load( config[[mm]]$genmatrix, envir=env )
        with(env,genmatrix)
    } )
longpats <- rownames(genmatrices[[1]])
stopifnot( all( sapply( genmatrices, function(gm) { all(rownames(gm)==longpats) } ) ) )

# read in counts
counts <- read.counts(opt$infile, leftwin=opt$leftwin, bases=genmatrices[[1]]@bases, longpats=rownames(genmatrices[[1]]) )
stopifnot( all( rownames(counts) == longpats ) )

projmatrix <- collapsepatmatrix( ipatterns=longpats, leftwin=leftwin(counts), fpatterns=colnames(counts) )

# order parameters (mutrates, selcoef, fixfn.params), then for each generator matrix
.param.map <- function (gm,type,params) {
    nparams <- sapply( genmatrices, function (x) { c( nmuts(x), nsel(x), length(fixparams(x)) ) } )
    p <- params[ do.call( seq, as.list( c(0,cumsum(nparams))[ (gm-1)*3 + pmatch(type,c("mutrates","selcoef","fixparams")) + c(0,1) ] + c(1,0) ) ) ]
    if ( pmatch(type,c("mutrates","selcoef","fixparams"))==3 ) { names(p) <- fixparams(genmatrices[[gm]]) }
    return(p)
}
initparam.list <- lapply( config[config$.models], "[", c("mutrates","selcoef","fixfn.params") )

cherry.transmats <- function (m1,m2,x) {
    mm <- m1[,rep(1:ncol(m1),ncol(m2))] * m2[,rep(1:ncol(m2),each=ncol(m1))] 
    if (!missing(x)) { mm <- sweep( mm, 1, x, "*" ) }
    return(mm)
}
twig.transmat <- function (m,x) { sweep(m,1,x,"*") }

# indices of nodes:
root.node <- get.root(config$tree)
row.node <- match( rowtaxa(counts), nodenames(config$tree) )
col.nodes <- match( coltaxa(counts), nodenames(config$tree) )

# find path from root to rowtaxa:
downpath <- c(row.node)
while( ! root.node %in% downpath ) { downpath <- c( get.parent(downpath[1],config$tree), downpath ) }
# list of active nodes that come off of the path from root to rowtaxa:
up.twigs <- lapply( selfname( col.nodes[ get.parent(col.nodes,config$tree) %in% downpath ]), function (x) projmatrix )
# list of other active nodes
up.active <- lapply( selfname( setdiff( col.nodes, up.twigs )), function (x) projmatrix )
while (length(up.active)>0) {
    up.cherries <- get.cherries(names(up.active),config$tree)
    uppair <- up.cherries[1,]
    # do something
    cat( "up: ", uppair, "\n" )
    transmat1 <- computetransmatrix( genmatrices[[uppair[1]]], up.active[[uppair[1]]], tlen=tlens[uppair[1]], time="fixed")
    transmat2 <- computetransmatrix( genmatrices[[uppair[2]]], up.active[[uppair[2]]], tlen=tlens[uppair[2]], time="fixed")
    newmat <- cherry.transmats( transmat1, transmat2 )
    # remove cherry
    up.active <- up.active[-uppair]
    newly.active <- get.parent( uppair[1], config$tree )
    if (get.parent(newly.active,config$tree) %in% downpath) {
        up.twigs <- c( up.twigs, list(newmat) )
        names(up.twigs)[length(up.twigs)] <- newly.active
    } else {
        up.active <- c( up.active, list(newmat) )
        names(up.active)[length(up.active)] <- newly.active
    }
}
# reorder up.twigs to match downpath
up.twigs <- up.twigs[ match( get.parent(up.twigs,config$tree), downpath ) ]
# start things off at the root (which is downpath[1])
cat("root: ", downpath[1], up.twigs[1], "\n" )
full.transmat <- twig.transmat( up.twigs[1], root.distrn )
# now move back down
for (k in seq_along(up.twigs[-1])) {
    cat("down: ", downpath[k+1], up.twigs[k+1], "\n" )
    # XXX
}
# and finish off
cat("final: ", row.node, "\n")


# Compute (quasi)-likelihood function using all counts -- multinomial as described in eqn:comp_like.
likfun <- function (params){
    # params are: mutrates, selcoef, fixparams
    # First, update genmatrices:
    mutrates.list <- lapply( seq_along(genmatrices), .param.map, type="mut", params=params )
    selcoef.list <- lapply( seq_along(genmatrices), .param.map, type="sel", params=params )
    fixparam.list <- lapply( seq_along(genmatrices), .param.map, type="fix", params=params )
    for (k in seq_along(genmatrices)) {
        genmatrices[[k]]@x <- do.call( update, c( list( G=genmatrices[[k]],mutrates=mutrates.list[[k]],selcoef=selcoef.list[[k]] ), as.list(fixparam.list[[k]]) ) )
    }
    # Now, peel XXX



    # this is collapsed transition matrix
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, tlen=1, time="fixed") # shape=params[length(params)], time="gamma" )
    # return POSITIVE log-likelihood
    ans <- sum( counts@counts * log(subtransmatrix) )
    if (!is.finite(ans)) { print(paste("Warning: non-finite likelihood with params:",params)) }
    return(ans)
}

initparams <- unlist( initparam.list )
stopifnot( length(initparams) == nmuts(genmatrix)+nsel(genmatrix)+length(fixparams(genmatrix)) )
lbs <- c( rep(1e-6,nmuts(genmatrix)), rep(-5,nsel(genmatrix)), rep(-Inf,length(fixparams(genmatrix))) )
ubs <- c( rep(2,nmuts(genmatrix)), rep(5,nsel(genmatrix)), rep(Inf,length(fixparams(genmatrix))) )
parscale <- c( 1e-3 * rep( mean(adhoc.mutrates),nmuts(genmatrix)), rep(.05,nsel(genmatrix)), rep(1,length(fixparams(genmatrix))) )

baseval <- likfun(initparams)
stopifnot( is.finite(baseval) )
optim.results <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(fnscale=(-1)*abs(baseval), parscale=parscale, maxit=opt$maxit) )


model <- new( "context",
             counts=counts,
             genmatrix=genmatrix,
             projmatrix=projmatrix,
             mutrates=optim.results$par[1:nmuts(genmatrix)],
             selcoef=optim.results$par[seq(nmuts(genmatrix)+1,length.out=nsel(genmatrix))],
             params=optim.results$par[seq(1+nmuts(genmatrix)+nsel(genmatrix),length.out=length(adhoc.fixparams))],
             results=optim.results,
             likfun=likfun,
             invocation=invocation
         )

save(model,file=opt$outfile)
