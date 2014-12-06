#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Infer parameters from paired counts file.
"

option_list <- list(
    # input/output
        make_option( c("-u","--basedir"), type="character", default=NULL, help="Directory to look for input in, and write output files to." ),
        make_option( c("-i","--infile"), type="character", default=NULL, help="Input file with tuple counts, tab-separated, with header 'reference', 'derived', 'count'. [default, looks in basedir]" ),
        make_option( c("-g","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix, or TRUE [default] to look for one. (otherwise, will compute)"),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" ),
        make_option( c("-j","--jobid"), type="character", default=formatC(1e6*runif(1),width=6,format="d",flag="0"), help="Unique job id. [default random]"),
    # context and pattern size
        make_option( c("-w","--shortwin"), type="integer", default=1, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--leftwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rightwin"), type="integer", default=2, help="Size of right-hand context. [default \"%default\"]" ),
        make_option( c("-k","--patlen"), type="integer", default=1, help="Include mutation rates for all tuples of this length. [default \"%default\"]" ),
        make_option( c("-d","--boundary"), type="character", default="none", help="Boundary conditions for generator matrix. [default \"%default\"]"),
        make_option( c("-y","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases in computing generator matrix. [default \"%default\"]" ),
    # mcmc run lengths
        make_option( c("-n","--nbatches"), type="integer", default=20, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-s","--stepscale"), type="numeric", default=1e-3, help="Scale of proposal steps for Metropolis algorithm relative to point estimates. [default \"%default\"]" ),
    # prior parameters
        make_option( c("-t","--shapemean"), type="double", default=1, help="Prior mean on shape parameter for gamma time. [default \"%default\"]" ),
        make_option( c("-m","--singlemean"), type="double", default=1, help="Prior mean on single base mutation rates. [default \"%default\"]" ),
        make_option( c("-z","--doublemean"), type="double", default=1, help="Prior mean on double base mutation rates. [default \"%default\"]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$infile) & is.null(opt$basedir)) { stop("No input file.  Run\n  bcells-inference.R -h\n for help.\n") }
print(opt) # this will go in the pbs log
attach(opt)
options(error = quote({dump.frames(to.file = TRUE); q()}))

longwin <- leftwin+shortwin+rightwin

if (substr(basedir,nchar(basedir),nchar(basedir)) %in% c("/","\\")) { basedir <- substr(basedir,1,nchar(basedir)-1) }
if (is.null(opt$infile)) { infile <- paste(basedir,"/", basedir, ".tuples.",longwin,".",leftwin,".counts",sep='') }
if (!file.exists(infile)) { stop("Cannot read file ", infile) }


if (gmfile=="TRUE") { gmfile <- paste(paste("genmatrices/genmatrix",longwin,boundary,meanboundary,patlen,sep="-"),".RData",sep='') }

if (logfile=="") {
    logfile <- paste(infile,".Rout",sep='')
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message")
    sink(file=logcon, type="output", split=interactive())   # send both to log file
}

scriptdir <- "../"
source(paste(scriptdir,"context-inference-fns.R",sep=''))

require(mcmc)

basename <- paste(basedir,"/win-",leftwin,"-",shortwin,"-",rightwin,"-",patlen,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
resultsfile <- paste( basename ,"-results.tsv",sep='')
plotfile <- paste( basename ,"-plot",sep='')

date()
cat("basename: ", basename, "\n")
print(opt)


if (file.exists(gmfile)) {
    load(gmfile)
} else {
    stop("Can't find generator matrix in ", gmfile, " -- provide file name exactly?")
}
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, shortwin=shortwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE, time="gamma" )

# read in counts (produced with count-paired-tuples.py)
count.table <- read.table(infile,header=TRUE,stringsAsFactors=FALSE)
counts <- Matrix(0,nrow=nrow(genmatrix),ncol=ncol(projmatrix))
rownames(counts) <- rownames(genmatrix)
colnames(counts) <- colnames(projmatrix)
stopifnot( all( count.table$reference %in% rownames(genmatrix) ) & all(count.table$derived %in% colnames(projmatrix)) ) 
counts[cbind( match(count.table$reference,rownames(genmatrix)), match(count.table$derived,colnames(projmatrix)) )] <- count.table$count
initcounts <- rowSums(counts)

# ad-hoc estimate
adhoc <- countmuts(counts=counts,mutpats=mutpats,leftwin=leftwin)
adhoc <- adhoc[1,]/adhoc[2,]

nmuts <- length(mutpats)
nsel <- length(selpats)
stopifnot(nsel==0)
# (quasi)-likelihood function using all counts -- multinomial
likfun <- function (params) {
    # params are: mutrates*tlen, shape
    genmatrix@x <- update(genmatrix,mutrates=params[1:nmuts],selcoef=numeric(0))
    # this is collapsed transition matrix
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, tlen=1, time="fixed") # shape=params[length(params)], time="gamma" )
    # return NEGATIVE log-likelihood 
    ans <- (-1) * sum( counts * log(subtransmatrix) ) + (params[length(params)]-1)^2
    if (!is.finite(ans)) { print(params) }
    return(ans)
}

mutlens <- sapply( mutpats, function (x) { nchar(x[[1]][1]) } )
# prior parameters
priormeans <- rep(NA,length(mutpats))
priormeans[mutlens==1] <- singlemean
priormeans[mutlens==2] <- doublemean
# posterior: multinomial (but, overlap!)
lud <- function (params) {
    # params are: mutrates*tlen, shape
    shape <- params[length(params)]
    mutrates <- params[1:nmuts]
    if (any(mutrates<0)) {
        return( -Inf )
    } else {
        genmatrix@x <- update(genmatrix,mutrates=mutrates,selcoef=numeric(0))
        subtransmatrix <- computetransmatrix( genmatrix, projmatrix, shape=shape, time="gamma" )
        # return POSITIVE log-posterior
        return( sum( counts * log(subtransmatrix) ) - sum(mutrates/priormeans) - shape/shapemean )
    }
}


# point estimates
initpar <- c(adhoc,1)
lbs <- c( rep(1e-6,nmuts), .01 )
ubs <- c( rep(2,nmuts), 5 )
parscale <- c( rep(mean(adhoc),length(adhoc)), 1 )*1e-3

baseval <- likfun(initpar)
stopifnot( is.finite(baseval) )
adhoc.ans <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(fnscale=abs(baseval), parscale=parscale, maxit=100) )
if (adhoc.ans$convergence!=0) { warning("optim failed to converge.") }
point.estimate <- adhoc.ans$par
names(point.estimate) <- c( sapply(mutpats,function(x){paste(sapply(x,paste,collapse='->'),collapse='/')}), "shape" )

# ###
# # check out likelihood function
# likfun(initpar)
# likfun(initpar*1.2)
# gradest(likfun,initpar*1.05)
# pcounts <- function (params) { predictcounts(shortwin, leftwin, rightwin, initcounts=rowSums(counts), mutrates=params[1:nmuts], selcoef=numeric(0), scale=params[length(params)], genmatrix=genmatrix, projmatrix=projmatrix, time="gamma" ) }
# ip <- pcounts(initpar)
# ip2 <- pcounts(10*initpar)
# layout(t(1:2))
# plot(as.vector(counts),as.vector(ip)); abline(0,1)
# plot(as.vector(counts),as.vector(ip2)); abline(0,1)
# likfun.parts <- function (params) {
#     genmatrix@x <- update(genmatrix,mutrates=params[1:nmuts],selcoef=numeric(0))
#     subtransmatrix <- computetransmatrix( genmatrix, projmatrix, shape=params[length(params)], time="gamma" )
#     return( (-1) * ( counts * log(subtransmatrix) ) )
# }
# plot( likfun.parts(initpar), likfun.parts(initpar*10) )
# plot( rowSums(likfun.parts(initpar)), rowSums(likfun.parts(initpar*10)) )

# bayesian
#  mrun.parjob <- mcparallel( metrop( lud, initial=random.ans.par[-length(random.ans.par)], nbatch=nbatches, blen=blen, scale=stepscale ) )
#  mrun <- mccollect(mcrun.parjob)
stepscale <- c( rep(mean(point.estimate[1:nmuts]),nmuts), point.estimate[nmuts+1] ) * stepscale
mrun <- metrop( lud, initial=point.estimate, nbatch=nbatches, blen=blen, scale=stepscale )

# look at observed/expected counts
pcounts <- function (params) { predictcounts(shortwin, leftwin, rightwin, initcounts=rowSums(counts), mutrates=params[1:nmuts], selcoef=numeric(0), scale=params[length(params)], genmatrix=genmatrix, projmatrix=projmatrix, time="gamma" ) }
expected <- pcounts(point.estimate)

# look at observed/expected counts in smaller windows
cwin <- min(2,shortwin); lrcwin <- min(1,leftwin,rightwin)
subcounts <- projectcounts( counts, new.leftwin=lrcwin, new.shortwin=cwin, new.longwin=lrcwin+cwin+lrcwin )
subexpected <- projectcounts( expected, new.leftwin=lrcwin, new.shortwin=cwin, new.longwin=lrcwin+cwin+lrcwin )

date()
cat("done with computation.\n")
cat("saving to: ", datafile, "\n")

save( opt, counts, genmatrix, projmatrix, subtransmatrix, lud, likfun, adhoc.ans, point.estimate, initpar, singlemean, doublemean, shapemean, expected, cwin, subcounts, subexpected, mrun, shortwin, leftwin, rightwin, longwin, patlen, nmuts, file=datafile )

# plot (long) counts
pdf(file=paste(plotfile,"-longcounts.pdf",sep=''),width=11, height=8, pointsize=10)
layout(matrix(1:min(ncol(counts),6),nrow=2))
for (k in 1:ncol(counts)) {
    lord <- order( expected[,k] )
    plot( counts[lord,k], type='n', xaxt='n', xlab='', ylim=range(c(unlist(as.matrix(counts)),unlist(as.matrix(expected)))), ylab='counts', main=colnames(counts)[k] )
    axis(1,at=1:nrow(counts),labels=rownames(counts)[lord],las=3)
    points( counts[lord,k] )
    lines(expected[lord,k],col='red' )
    if (k==1) legend("topleft",legend=c("fitted"),lty=1,col=c("red"))
}
dev.off()

# plot (shorter) counts 
pdf(file=paste(plotfile,"-shortcounts.pdf",sep=''),width=6, height=4, pointsize=10)
layout(matrix(1:min(ncol(subcounts),6),nrow=2))
par(mar=c(4,2,1,1)+.1)
for (k in 1:ncol(subcounts)) {
    lord <- order( subexpected[,k] )
    plot( subcounts[lord,k], xaxt='n', xlab='', main=colnames(subcounts)[k], log='y', ylim=1+range(c(as.matrix(subcounts),as.matrix(subexpected))), ylab='' )
    axis(1,at=1:nrow(subcounts),labels=rownames(subcounts)[lord],las=3)
    lines(subexpected[lord,k],col='red')
    legend("topleft",legend='fitted',lty=1,col='red')
}
dev.off()

# residuals of (shorter) counts 
pdf(file=paste(plotfile,"-shortresids.pdf",sep=''),width=6, height=4, pointsize=10)
par(mar=c(4,4,2,1)+.1)
cols <- "red"
subresids <- (subcounts-subexpected)/sqrt(subexpected)
z <- t( t( as.vector( subresids ) ) )
rownames(z) <- paste( rownames(subcounts)[row(subcounts)], colnames(subcounts)[col(subcounts)], sep="->" )
plot( seq_along(z), z, type='n', xlim=c(1,length(z))+c(-1,+1)*length(z)/20 )
text( rank(z), z, labels=rownames(z) )
abline(h=c(-2,0,2), lty=c(3,1,3) )
dev.off()

# observed vs expected
pdf(file=paste(plotfile,"-obs-exp.pdf",sep=''),width=6, height=4, pointsize=10)
layout(t(1:2))
par(mar=c(4,4,2,1)+.1)
plot( as.vector(expected), as.vector(counts), log='xy', xlab="true expected counts", ylab="counts" )
abline(0,1)
plot( as.vector(expected), as.vector(counts), xlab="true expected counts", ylab="counts" )
abline(0,1)
dev.off()


print(format(Sys.time(),"%Y-%m-%d-%H-%M"))
sink(NULL); close(logcon)
