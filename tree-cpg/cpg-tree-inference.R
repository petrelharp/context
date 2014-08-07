#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Infer parameters from output of sim-tree-cpg.R .\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file containing simulation." ),
        make_option( c("-w","--shortwin"), type="integer", default=1, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--leftwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rightwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-n","--nbatches"), type="integer", default=20, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-s","--stepscale"), type="numeric", default=1e-4, help="Scale of proposal steps for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-m","--mmean"), type="double", default=1, help="Prior mean on single base mutation rates. [default \"%default\"]" ),
        make_option( c("-c","--cpgmean"), type="double", default=1, help="Prior variance on CpG rate. [default \"%default\"]" ),
        make_option( c("-p","--pprior"), type="double", default=1, help="Parameter for Dirichlet prior on base frequencies. [default \"%default\"]" ),
        make_option( c("-v","--tprior"), type="double", default=.5, help="Parameter for Beta prior on branch length. [default \"%default\"]" ),
        make_option( c("-d","--boundary"), type="character", default="none", help="Boundary conditions for generator matrix. [default \"%default\"]"),
        make_option( c("-y","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases in computing generator matrix. [default \"%default\"]" ),
        make_option( c("-g","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix, or TRUE [default] to look for one. (otherwise, will compute)"),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)
options(error=traceback)

longwin <- leftwin+shortwin+rightwin

if (gmfile=="TRUE") { gmfile <- paste(paste("genmatrices/genmatrix",longwin,boundary,meanboundary,sep="-"),".RData",sep='') }

if (is.null(infile)) { cat("Run\n  cpg-tree-inference.R -h\n for help.") }

if (logfile!="") {
    logfile <- gsub(".RData",".Rout",infile,fixed=TRUE)
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message", split=interactive()) 
    sink(file=logcon, type="output", split=interactive())   # send both to log file
}

scriptdir <- "../"
source(paste(scriptdir,"context-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))

require(mcmc)

load(infile)
basedir <- gsub(".RData","",infile,fixed=TRUE)
if (!file.exists(basedir)) { dir.create(basedir) }
basename <- paste(basedir,"/win-",leftwin,"-",shortwin,"-",rightwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
resultsfile <- paste( basename ,"-results.tsv",sep='')
plotfile <- paste( basename ,"-plot",sep='')

if (file.exists(gmfile)) {
    load(gmfile)
} else {
    if (meanboundary>0) {
        genmatrix <- meangenmatrix( leftwin=1, rightwin=1, patlen=longwin, mutpats=mutpats, selpats=list(), mutrates=mutrates*tlen[1], selcoef=numeric(0), boundary=boundary )
    } else {
        genmatrix <- makegenmatrix( patlen=longwin, mutpats=mutpats, selpats=list(), mutrates=mutrates*tlen[1], selcoef=numeric(0), boundary=boundary )
    }
}
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

counts <- list(
            "1.2"=counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[1]]$finalseq, simseqs[[2]]$finalseq, leftwin=leftwin ),
            "2.1"=counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[2]]$finalseq, simseqs[[1]]$finalseq, leftwin=leftwin )
        )
# want only patterns that overlap little
initcounts <- lapply( counts, rowSums )
nonoverlapping <- lapply( seq_along(counts), function (k) ( leftchanged(rownames(counts[[k]]),colnames(counts[[k]]),leftwin=leftwin,shortwin=shortwin) & (initcounts[[k]]>0) ) )
nov.counts <- lapply(seq_along(counts), function (k) counts[[k]][nonoverlapping[[k]]] )

# move from base frequencies (what we estimate) to pattern frequencies
nmuts <- length(mutpats)
nfreqs <- length(initfreqs)
npats <- nrow(genmatrix)
patcomp <- apply( do.call(rbind, strsplit(rownames(genmatrix),'') ), 2, match, bases )  # which base is at each position in each pattern
likfun <- function (params) {
    # params are: tlen[1]/sum(tlen), sum(tlen)*mutrates, initfreqs
    branchlens <- c(params[1],1-params[1])
    mutrates <- params[1+(1:nmuts)]
    initfreqs <- params[1+nmuts+(1:nfreqs)]
    initfreqs <- initfreqs/sum(initfreqs)
    patfreqs <- initfreqs[patcomp]
    dim(patfreqs) <- dim(patcomp)
    patfreqs <- apply( patfreqs, 1, prod )
    # these are collapsed transition matrix
    updownbranch <- list(  # note "up" branch is from simpler summaries
            getupdowntrans( genmatrix, projmatrix, mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=rev(branchlens) ),
            getupdowntrans( genmatrix, projmatrix, mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=branchlens )
        )
    if (any(sapply(updownbranch,function(x) any(!is.numeric(x))|any(x<0)))) { browser() }
    # return negative log-likelihood plus a penalty to keep initfreqs summing to (almost) 1
    return( 
                (-1) * ( sum( counts[[1]] * log(updownbranch[[1]]) ) + sum( counts[[2]] * log(updownbranch[[2]]) ) ) 
                + 100*(sum(initfreqs)-1)^2
            )
}

# using only nonoverlapping counts, plus priors -- indep't trials. (multinomial)
mmeans <- c( rep(mmean,nmuts-1), cpgmean )
ppriors <- rep( pprior, nfreqs )
tpriors <- rep( tprior, 2 )
lud <- function (params) {
    # params are: tlen[1]/sum(tlen), sum(tlen)*mutrates, initfreqs[-length(initfreqs)]
    branchlens <- c(params[1],1-params[1])
    mutrates <- params[1+(1:nmuts)]
    initfreqs <- params[1+nmuts+(1:(nfreqs-1))]
    initfreqs <- c(initfreqs,1-sum(initfreqs))
    patfreqs <- initfreqs[patcomp]
    dim(patfreqs) <- dim(patcomp)
    patfreqs <- apply( patfreqs, 1, prod )
    if (any(mutrates<0) | any(initfreqs<0)) {
        return( -Inf )
    } else {
        # only do in one direction... ?
        updownbranch <- getupdowntrans( genmatrix, projmatrix, mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), initfreqs=patfreqs, tlens=rev(branchlens) )
        # return (positive) log-posterior
        return( 
                # (-1)*sum(updownbranch[nonoverlapping[[1]]]) +
                sum( nov.counts[[1]] * log(updownbranch[nonoverlapping[[1]]]) ) 
                + sum( (tpriors-1)*log(branchlens) )
                - sum(mmeans*mutrates) 
                + sum( (ppriors-1)*log(initfreqs) )
            )
    }
}

# point estimates
rand.initfreqs <- 5*rexp(length(initfreqs)); rand.initfreqs <- rand.initfreqs/sum(rand.initfreqs)
initpar <- c( runif(1), 2 * runif( nmuts ) * mean(mutrates) * sum(tlen), rand.initfreqs ) # random init
truth <- c( tlen[1]/sum(tlen), mutrates * sum(tlen), initfreqs )  # truth
names(truth) <- c( "rel.branchlen", paste("tmut:", unlist( sapply( sapply( mutpats, lapply, paste, collapse="->" ), paste, collapse=" | " ) ) ), names(initfreqs) )
lbs <- c( 1e-6, rep(0,nmuts), rep(1e-6,length(initfreqs)) )
ubs <- c( 1, rep(20,nmuts), rep(1,length(initfreqs)) )
# # do in parallel:
# cheating.parjob <- mcparallel( optim( par=truth, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(trace=3,fnscale=likfun(truth)) ) ),
# random.parjob <- mcparallel( optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(trace=3,fnscale=likfun(truth)) ) )
# cheating.ans <- mccollect( cheating.parjob, wait=TRUE )
# random.ans <- mccollect( random.parjob, wait=TRUE )
cheating.ans <- optim( par=truth, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(trace=3,fnscale=abs(likfun(truth))) )
random.ans <- optim( par=initpar, fn=likfun, method="L-BFGS-B", lower=lbs, upper=ubs, control=list(trace=3,fnscale=abs(likfun(truth))) )
# renormalize the initial frequencies
cheating.ans.par <- cheating.ans$par; cheating.ans.par[1+nmuts+(1:nfreqs)] <- cheating.ans.par[1+nmuts+(1:nfreqs)] / sum( cheating.ans.par[1+nmuts+(1:nfreqs)] )
random.ans.par <- random.ans$par; random.ans.par[1+nmuts+(1:nfreqs)] <- random.ans.par[1+nmuts+(1:nfreqs)] / sum( random.ans.par[1+nmuts+(1:nfreqs)] )

estimates <- data.frame( rbind(init=initpar, ans=random.ans.par, cheating=cheating.ans.par, truth=truth ) )
colnames(estimates) <- c( "branchlen", paste("muttime",seq_along(mutrates),sep=''), names(initfreqs) )
estimates$likfun <- apply( estimates, 1, likfun )
write.table( estimates, file=resultsfile, quote=FALSE, sep="\t" )

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
                    predicttreecounts( shortwin, leftwin, rightwin, initcounts=rowSums(counts[[1]]), mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), genmatrix=genmatrix, projmatrix=projmatrix, initfreqs=initfreqs, tlens=rev(branchlens) ),
                    predicttreecounts( shortwin, leftwin, rightwin, initcounts=rowSums(counts[[2]]), mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), genmatrix=genmatrix, projmatrix=projmatrix, initfreqs=initfreqs, tlens=branchlens )
                )
    } )
names(all.expected) <- rownames(estimates)

# look at observed/expected counts in smaller windows
cwin <- min(2,shortwin); lrcwin <- min(1,leftwin,rightwin)
subcounts <- lapply( counts, function (x) 
        projectcounts( leftwin=leftwin, countwin=cwin, lcountwin=lrcwin, rcountwin=lrcwin, counts=x ) )
all.subexpected <- lapply( all.expected, lapply, function (x)
        projectcounts( leftwin=leftwin, countwin=cwin, lcountwin=lrcwin, rcountwin=lrcwin, counts=x ) )

save( opt, counts, genmatrix, projmatrix, subtransmatrix, lud, likfun, truth, cheating.ans, random.ans, estimates, initpar, nonoverlapping, nov.counts, mmeans, ppriors, tpriors, all.expected, cwin, subcounts, all.subexpected, mrun, shortwin, leftwin, rightwin, nmuts, nfreqs, npats, patcomp, file=datafile )

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
