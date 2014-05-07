#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Infer parameters from paired counts file.
"

option_list <- list(
        make_option( c("-u","--basedir"), type="character", default=NULL, help="Directory to look for input in, and write output files to." ),
        make_option( c("-i","--infile"), type="character", default=NULL, help="Input file with tuple counts, tab-separated, with header 'reference', 'derived', 'count'. [default, looks in basedir]" ),
        make_option( c("-w","--win"), type="integer", default=1, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--lwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-k","--patlen"), type="integer", default=1, help="Include mutation rates for all tuples of this length. [default \"%default\"]" ),
        make_option( c("-n","--nbatches"), type="integer", default=20, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-s","--stepscale"), type="numeric", default=1e-4, help="Scale of proposal steps for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-m","--mmean"), type="double", default=1, help="Prior mean on single base mutation rates. [default \"%default\"]" ),
        make_option( c("-c","--cpgmean"), type="double", default=1, help="Prior variance on CpG rate. [default \"%default\"]" ),
        make_option( c("-z","--gcbias"), type="double", default=1, help="Prior mean for strength of GC bias. [default \"%default\"]" ),
        make_option( c("-p","--pprior"), type="double", default=1, help="Parameter for Dirichlet prior on base frequencies. [default \"%default\"]" ),
        make_option( c("-v","--tprior"), type="character", default="c(1,1,1,1)", help="Parameters for Dirichlet priors on relative branch lengths. [default \"%default\"]" ),
        make_option( c("-d","--boundary"), type="character", default="none", help="Boundary conditions for generator matrix. [default \"%default\"]"),
        make_option( c("-y","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases in computing generator matrix. [default \"%default\"]" ),
        make_option( c("-g","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix, or TRUE [default] to look for one. (otherwise, will compute)"),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.character(opt$tprior)) { opt$tprior <- eval(parse(text=opt$tprior)) }
if (length(opt$tprior)==1) { opt$tprior <- rep(opt$tprior,4) }
if (is.null(opt$infile) & is.null(opt$basedir)) { cat("Run\n  bcells-inference.R -h\n for help.\n") }
attach(opt)
options(error=traceback)

winlen <- lwin+win+rwin

if (is.null(opt$infile)) { infile <- paste(basedir,"/tuples.",winlen,".",lwin,".counts",sep='') }

if (gmfile=="TRUE") { gmfile <- paste(paste("genmatrices/genmatrix",winlen,boundary,meanboundary,patlen,sep="-"),".RData",sep='') }

if (logfile!="") {
    logfile <- paste(infile,".Rout",sep='')
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message", split=interactive()) 
    sink(file=logcon, type="output", split=interactive())   # send both to log file
}

scriptdir <- "../"
source(paste(scriptdir,"codon-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))

require(mcmc)

basename <- paste(basedir,"/win-",lwin,"-",win,"-",rwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
resultsfile <- paste( basename ,"-results.tsv",sep='')
plotfile <- paste( basename ,"-plot",sep='')

if (file.exists(gmfile)) {
    load(gmfile)
} else {
    if (meanboundary>0) {
        genmatrix <- meangenmatrix( lwin=1, rwin=1, patlen=winlen, mutpats=mutpats, selpats=selpats, mutrates=rep(1,length(mutpats)), selcoef=rep(1,length(selpats)), boundary=boundary )
    } else {
        genmatrix <- makegenmatrix( patlen=winlen, mutpats=mutpats, selpats=selpats, mutrates=rep(1,length(mutpats)), selcoef=rep(1,length(selpats)), boundary=boundary )
    }
}
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

# read in counts (produced with count-paired-tuples.py)
count.table <- read.table(infile,header=TRUE,stringsAsFactors=FALSE)
counts <- Matrix(0,nrow=nrow(genmatrix),ncol=ncol(projmatrix))
rownames(counts) <- rownames(genmatrix)
colnames(counts) <- colnames(projmatrix)
stopifnot( all( count.table$reference %in% rownames(genmatrix) ) & all(count.table$derived %in% colnames(projmatrix)) ) 
counts[cbind( match(count.table$reference,rownames(genmatrix)), match(count.table$derived,colnames(projmatrix)) )] <- count.table$count

# ad-hoc estimate
adhoc <- countmuts(counts=counts,mutpats=mutpats,lwin=lwin)
adhoc <- adhoc[1,]/adhoc[2,]

nmuts <- length(mutpats)
nsel <- length(selpats)
stopifnot(nsel==0)
# (quasi)-likelihood function using all counts -- binomial
likfun <- function (params) {
    cat(".")
    # params are: mutrates*tlen
    genmatrix@x <- update(genmatrix,mutrates=params,selcoef=numeric(0))
    # this is collapsed transition matrix
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix )
    # return negative log-likelihood 
    ans <- (-1) * sum( counts * log(subtransmatrix) )
    if (!is.finite(ans)) { print(params) }
    return(ans)
}


# point estimates
initpar <- .1 + runif( nmuts )/4
lbs <- rep(1e-8,nmuts)
ubs <- rep(20,nmuts)

stopifnot( is.finite(likfun(initpar)) )

require(parallel)
pjob <- mcparallel( { optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs ) } )
random.ans <- mccollect(pjob)[[1]]

random.ans <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(trace=3) )
stopifnot(random.ans$convergence==0)

point.estimate <- random.ans$par
names(point.estimate) <- sapply(mutpats,function(x){paste(sapply(x,paste,collapse='->'),collapse='/')})

# is it symmetric under reverse-complement?
revmutrates <- point.estimate[reverse.complement(mutpats)]
plot(point.estimate,revmutrates)

predicted.counts <- predictcounts(win=win,lwin=lwin,rwin=rwin,initcounts=rowSums(counts),mutrates=point.estimate,selcoef=numeric(0),genmatrix=genmatrix,projmatrix=projmatrix)
p.vals <- ppois( as.matrix(counts), lambda=as.matrix(predicted.counts) )

# look at residuals
layout(t(1:2))
plot( (predicted.counts - counts)/sqrt(predicted.counts), p.vals, cex=pmin(3,predicted.counts/mean(predicted.counts))  )
plot(as.vector(predicted.counts),as.vector(counts),xlab='predicted',ylab='observed',cex=pmin(3,abs(log(p.vals/(1-p.vals))/10)),col=1+(p.vals>0),log='xy')
abline(0,1)


# bayesian
#  note we deal with initial freqs not summing to 1 differently from in optim( ) -- need to remove last entry.
#  mrun.parjob <- mcparallel( metrop( lud, initial=random.ans.par[-length(random.ans.par)], nbatch=nbatches, blen=blen, scale=stepscale ) )
#  mrun <- mccollect(mcrun.parjob)
mrun <- metrop( lud, initial=random.ans.par[-length(random.ans.par)], nbatch=nbatches, blen=blen, scale=stepscale )

# look at observed/expected counts
if (is.null(names(initfreqs))) { names(initfreqs) <- bases }
all.expected <- lapply( 1:nrow(estimates), function (k) {
            x <- unlist(estimates[k,])
            branchlens <- c(x[1],1-x[1])
            mutrates <- x[1+(1:nmuts)]
            initfreqs <- x[1+nmuts+(1:nfreqs)]
            initfreqs <- initfreqs/sum(initfreqs)
            list( 
                    predicttreecounts( win, lwin, rwin, initcounts=rowSums(counts[[1]]), mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), genmatrix=genmatrix, projmatrix=projmatrix, initfreqs=initfreqs, tlens=rev(branchlens) ),
                    predicttreecounts( win, lwin, rwin, initcounts=rowSums(counts[[2]]), mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), genmatrix=genmatrix, projmatrix=projmatrix, initfreqs=initfreqs, tlens=branchlens )
                )
    } )
names(all.expected) <- rownames(estimates)

# look at observed/expected counts in smaller windows
cwin <- min(2,win); lrcwin <- min(1,lwin,rwin)
subcounts <- lapply( counts, function (x) 
        projectcounts( lwin=lwin, countwin=cwin, lcountwin=lrcwin, rcountwin=lrcwin, counts=x ) )
all.subexpected <- lapply( all.expected, lapply, function (x)
        projectcounts( lwin=lwin, countwin=cwin, lcountwin=lrcwin, rcountwin=lrcwin, counts=x ) )

save( opt, counts, genmatrix, projmatrix, subtransmatrix, lud, likfun, truth, cheating.ans, random.ans, estimates, initpar, nonoverlapping, nov.counts, mmeans, ppriors, tpriors, all.expected, cwin, subcounts, all.subexpected, mrun, win, lwin, rwin, nmuts, nfreqs, npats, patcomp, file=datafile )

# plot (long) counts
pdf(file=paste(plotfile,"-longcounts.pdf",sep=''),width=6, height=4, pointsize=10)
layout(matrix(1:sum(sapply(counts,ncol)),nrow=2))
for (j in seq_along(counts)) {
    for (k in 1:ncol(counts[[j]])) {
        lord <- order( all.expected[["truth"]][[j]][,k] )
        plot( counts[[j]][lord,k], type='n', xaxt='n', xlab='', ylim=range(c(unlist(all.expected[["truth"]]),unlist(as.matrix(counts[[j]])),unlist(all.expected[["ans"]]))), ylab='counts', main=colnames(counts[[j]])[k] )
        axis(1,at=1:nrow(counts[[j]]),labels=rownames(counts[[j]])[lord],las=3)
        points( counts[[j]][lord,k], pch=j )
        lines(all.expected[["truth"]][[j]][lord,k],col='red', lty=j)
        lines(all.expected[["ans"]][[j]][lord,k],col='green', lty=j, lwd=2)
        lines(all.expected[["cheating"]][[j]][lord,k],col='blue',lty=j)
        if (k==1) legend("topleft",legend=c("expected","estimated","cheating"),lty=1,col=c("red","green","blue"))
    }
}
dev.off()

# plot (shorter) counts 
pdf(file=paste(plotfile,"-shortcounts.pdf",sep=''),width=6, height=4, pointsize=10)
layout(matrix(1:sum(sapply(subcounts,ncol)),nrow=2))
cols <- rainbow(2+length(all.expected))[1:length(all.expected)]
for (j in seq_along(subcounts)) {
    for (k in 1:ncol(subcounts[[j]])) {
        lord <- order( all.subexpected[["truth"]][[j]][,k] )
        plot( subcounts[[j]][lord,k], xaxt='n', xlab='', main=colnames(subcounts[[j]])[k], log='y' )
        axis(1,at=1:nrow(subcounts[[j]]),labels=rownames(subcounts[[j]])[lord],las=3)
        invisible( lapply(seq_along(all.subexpected),function(y) { lines(all.subexpected[[y]][[j]][lord,k],col=cols[y]) } ) )
        legend("topleft",legend=names(all.subexpected),lty=1,col=cols)
    }
}
dev.off()

# residuals of (shorter) counts 
pdf(file=paste(plotfile,"-shortresids.pdf",sep=''),width=6, height=4, pointsize=10)
layout(matrix(seq_along(subcounts)))
cols <- rainbow(2+length(all.expected))[1:length(all.expected)]
all.subresids <- lapply( all.subexpected, function (x) mapply(function(u,v) (u-v)/sqrt(v),x,subcounts) )
for (j in seq_along(counts)) {
    z <- sapply( lapply( all.subresids[c("truth","ans","cheating")], "[[", j ), as.vector )
    rownames(z) <- paste( rownames(subcounts[[j]])[row(subcounts[[j]])], colnames(subcounts[[j]])[col(subcounts[[j]])], sep="->" )
    plot( 0, type='n', xlab="", ylab="normalized residuals", xlim=c(0,ncol(z)+1), ylim=range(z), xaxt='n' )
    text( jitter(col(z),fac=2), z, labels=rownames(z), col=col(z) )
    axis(1, at=1:ncol(z), labels=colnames(z) )
}
dev.off()

# observed vs expected
pdf(file=paste(plotfile,"-obs-exp.pdf",sep=''),width=6, height=4, pointsize=10)
layout(seq_along(counts))
for (j in seq_along(counts)) {
    plot( as.vector(all.expected[["truth"]][[j]]), as.vector(counts[[j]]), log='xy', xlab="true expected counts", ylab="counts" )
    abline(0,1)
    points(as.vector(all.expected[["truth"]][[j]]), as.vector(all.expected[["cheating"]][[j]]), col='red', pch=20 )
    points(as.vector(all.expected[["truth"]][[j]]), as.vector(all.expected[["ans"]][[j]]), col='green', pch=20, cex=.5)
    legend("topleft",pch=c(1,20,20),col=c('black','red','green'),legend=c('observed','cheating','estimated'))
}
dev.off()

# NOT SO USEFUL since we only do a short mcmc run
# # plot each pairwise posterior marginal density
# pdf(file=paste(plotfile,"-mcmc.pdf",sep=''),width=6, height=4, pointsize=10)
# pairs( rbind( mrun$batch, truth ), col=c(rep(1,nrow(x)),2), pch=c(rep(1,nrow(x)),20), cex=c(rep(1,nrow(x)),2) )
# plot(as.data.frame(mrun$batch))
# layout(matrix(c(1,4,2,3),nrow=2))
# par(mar=c(4,4,0,0)+.1)
# for (k in 1:(ncol(mrun$batch)-1)) { 
#     for (j in (k+1):ncol(mrun$batch)) {
#         plot(mrun$batch[,j],mrun$batch[,k],xlab=colnames(estimates)[j],ylab=colnames(estimates)[k],pch=20,cex=.5,col=adjustcolor('black',.5), xlim=range(c(mrun$batch[,j],estimates[c("truth","ans"),j])), ylim=range(c(mrun$batch[,k],estimates[c("truth","ans"),k])))
#         points( estimates["truth",j], estimates["truth",k], pch=20, col='green' )
#         points( estimates["ans",j], estimates["ans",k], pch=20, col='red' )
#     }
# }
# legend("bottomright",pch=20,col=c('green','red'),legend=c("truth","estimated"))
# matplot( mrun$batch, type='l', col=1:length(truth), xlab="mcmc gens" )
# abline(h=truth, col=adjustcolor(1:length(truth),.5), lwd=2)
# abline(h=estimates["ans",], col=1:length(truth) )
# legend("topright",col=c(1:length(truth),1,adjustcolor(1,.5)), lwd=c(rep(1,length(truth)),1,2),legend=c(colnames(estimates)[1:length(truth)],"point estimate","truth"))
# dev.off()


print(format(Sys.time(),"%Y-%m-%d-%H-%M"))
