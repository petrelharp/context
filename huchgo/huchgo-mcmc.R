#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Run mcmc longer.\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file containing simulation." ),
        make_option( c("-w","--win"), type="integer", default=1, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--lwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-n","--nbatches"), type="integer", default=1000, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-c","--stepscale"), type="character", default="1e-4", help="Scale of proposal steps for mutation*time parameters for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-e","--initscale"), type="character", default="3e-3", help="Scale of proposal steps for initial frequencies for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-f","--treescale"), type="character", default="1e-4", help="Scale of proposal steps for relative tree branch lengths for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-s","--restart"), action="store_true", default=FALSE, help="Start a whole new MCMC run?" ),
        make_option( c("-d","--boundary"), type="character", default="none", help="Boundary conditions for generator matrix. [default \"%default\"]"),
        make_option( c("-y","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases in computing generator matrix. [default \"%default\"]" ),
        make_option( c("-g","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix, or TRUE [default] to look for one. (otherwise, will compute)"),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.character(opt$stepscale)) { opt$stepscale <- eval(parse(text=opt$stepscale)) }
if (is.character(opt$initscale)) { opt$initscale <- eval(parse(text=opt$initscale)) }
if (is.character(opt$treescale)) { opt$treescale <- eval(parse(text=opt$treescale)) }
attach(opt)

if (interactive()) { nbatches <- 100; blen <- 10; restart <- FALSE }

if (!'infile'%in%names(opt)) { stop("Run\n  cpg-inference.R -h\n for help.") }

# options(error=traceback)

scriptdir <- "../"
source(paste(scriptdir,"codon-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))
require(mcmc)

basedir <- gsub(".RData","",infile,fixed=TRUE)
load(infile)
if (!file.exists(basedir)) {
    dir.create(basedir)
}

basename <- paste(basedir,"/win-",lwin,"-",win,"-",rwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
plotfile <- paste( basename ,"-plot",sep='')
mcmcdatafiles <- list.files(path=basedir,pattern="-mcmc.*RData",full.names=TRUE)
mcmcnum <- 1+max(c(0,as.numeric(gsub(".*-mcmc-","",gsub(".RData","",mcmcdatafiles)))),na.rm=TRUE)

if (logfile=="") { logfile <- paste(basename,"-mcmc-run-",mcmcnum,".Rout",sep='') }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="output", split=interactive()) 
    if (!interactive()) { sink(file=logcon, type="message",append=TRUE) }
}

load(datafile)  # has mrun and previous things
if (length(mcmcdatafiles)>0) { load(grep(paste("-mcmc-",mcmcnum-1,".RData",sep=''),mcmcdatafiles,fixed=TRUE,value=TRUE)) }

########

# Inference.
winlen <- lwin+win+rwin
if (gmfile=="TRUE") { gmfile <- paste(paste("genmatrices/genmatrix",winlen,boundary,meanboundary,sep="-"),".RData",sep='') }
if (file.exists(gmfile)) {
    load(gmfile) 
} else {
    if (meanboundary>0) {
        genmatrix <- meangenmatrix( lwin=1, rwin=1, patlen=winlen, mutpats=mutpats, selpats=list(), mutrates=mutrates*tlen, selcoef=numeric(0), boundary=boundary )
    } else {
        genmatrix <- makegenmatrix( patlen=winlen, mutpats=mutpats, selpats=list(), mutrates=mutrates*tlen, selcoef=numeric(0), boundary=boundary )
    }
}
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin )

counts <- list(
            "1.2"=counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[1]]$finalseq, simseqs[[2]]$finalseq, lwin=lwin ),
            "2.1"=counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[2]]$finalseq, simseqs[[1]]$finalseq, lwin=lwin )
        )
# want only patterns that overlap little
initcounts <- lapply( counts, rowSums )
nonoverlapping <- lapply( seq_along(counts), function (k) ( leftchanged(rownames(counts[[k]]),colnames(counts[[k]]),lwin=lwin,win=win) & (initcounts[[k]]>0) ) )
nov.counts <- lapply(seq_along(counts), function (k) counts[[k]][nonoverlapping[[k]]] )

# move from base frequencies (what we estimate) to pattern frequencies
nmuts <- length(mutpats)
nfreqs <- length(initfreqs)
npats <- nrow(genmatrix)
patcomp <- apply( do.call(rbind, strsplit(rownames(genmatrix),'') ), 2, match, bases )  # which base is at each position in each pattern
# using only nonoverlapping counts, plus priors -- indep't poisson.
lud <- function (params) {
    # params are: tlen[1]/sum(tlen), sum(tlen)*mutrates, initfreqs
    branchlens <- c(params[1],1-params[1])
    mutrates <- params[1+(1:nmuts)]
    initfreqs <- params[1+nmuts+(1:nfreqs)]
    patfreqs <- initfreqs[patcomp]
    dim(patfreqs) <- dim(patcomp)
    patfreqs <- apply( patfreqs, 1, prod )
    if (any(mutrates<0) | any(initfreqs<0)) {
        return( -Inf )
    } else {
        # only do in one direction... ?
        updownbranch <- getupdowntrans( genmatrix, projmatrix, 
            mutrates=list(mutrates,mutrates), selcoef=list(numeric(0),numeric(0)), 
            initfreqs=patfreqs, tlens=rev(branchlens) )
        meancounts <- initcounts[[1]] * updownbranch
        # return (positive) log-posterior
        return( 
                (-1)*sum(meancounts[nonoverlapping[[1]]]) 
                + sum( nov.counts[[1]] * log(meancounts[nonoverlapping[[1]]]) ) 
                + sum( (tpriors-1)*log(branchlens) )
                - sum(mmeans*mutrates) 
                + sum( (ppriors-1)*log(initfreqs) )
            )
    }
}

# construct matrix for proposal distribution for mcmc steps:
#  First, need base frequencies to sum to 1.
# we want (root distr'n for x) * (transition rate x -> y) to stay the same,
#  so adjust proposed changes to transition rates
#  by the proposed changes in root distr'n.
# The proposals are, in metrop( ), controlled by scale %*% (standard normal vector).
scalemat <- diag(1+nmuts+nfreqs)
# proposed base freq changes will have zero sum:
scalemat[1+nmuts+(1:nfreqs),1+nmuts+(1:nfreqs)] <- ( diag(nfreqs) - 1/nfreqs )
changes <- mutpatchanges(mutpats)
fac <- 1 / ( mean(truth[1+nmuts+1:nfreqs]) / mean( truth[1+1:nmuts] ) )
for (k in seq_along(mutpats)) {
    frombases <- match( changes[[k]][,1], bases )
    scalemat[1+k,1+nmuts+frombases] <- (-1)*fac
}
scalemat <- scalemat * c( treescale, rep_len( stepscale, nmuts ), rep_len( initscale, nfreqs ) )
rownames(scalemat) <- colnames(scalemat) <- c( 'tree', sapply(lapply(mutpats[1:13],"[[",1),paste,collapse='-'), bases )

if (restart) {
    # mrun <- metrop( lud, initial=random.ans$par[-length(random.ans$par)], nbatch=nbatches, blen=blen, scale=scalemat )
    if (is.finite(lud(random.ans$par[-length(random.ans$par)]))) {
        mrun <- metrop( lud, initial=random.ans$par, nbatch=nbatches, blen=blen, scale=scalemat )
    } else {  # the point estimate broke; cheat to get a good starting point
        mrun <- metrop( lud, initial=truth, nbatch=nbatches, blen=blen, scale=scalemat )
    }
} else {
    mrun <- metrop( mrun, nbatch=nbatches, blen=blen, scale=scalemat )
}
colnames(mrun$batch) <- names(truth[-length(truth)])

save( lwin, win, rwin, lud, mrun, initcounts, nov.counts, nonoverlapping, file=paste(basename,"-mcmc-",mcmcnum,".RData",sep='') )

pdf(file=paste(plotfile,"-mcmc-",mcmcnum,".pdf",sep=''),width=6, height=4, pointsize=10)
subs <- seq.int(1, nrow(mrun$batch), by=max(1,floor(nrow(mrun$batch)/1000)) )
matplot( mrun$batch[subs,], type='l', col=1:length(truth) )
abline(h=truth[-length(truth)], col=adjustcolor(1:length(truth),.5), lwd=2)
abline(h=estimates["ans",], col=1:(length(truth)-1) )
legend("topright", col=c(1:length(truth),1,adjustcolor(1,.5)), lwd=c(rep(1,length(truth)-1),1,2),legend=c(colnames(estimates)[1:(length(truth)-1)],"point estimate","truth"))
dev.off()

pdf(file=paste(plotfile,"-mcmc-",mcmcnum,"-pairwise.pdf",sep=''),width=6,height=6,pointsize=10)
mutlabels <- c("branchlen", paste("mut:", unlist( sapply( sapply( mutpats, lapply, paste, collapse="->" ), paste, collapse=" | " ) ) ), names(initfreqs) )
subs <- seq.int(1, nrow(mrun$batch), by=max(1,floor(nrow(mrun$batch)/1000)) )
pairs( rbind( mrun$batch[subs,], truth[-length(truth)] ), col=c( rep(adjustcolor("black",.1),length(subs)), adjustcolor("red",1)), pch=20, cex=c( rep(.25,length(subs)), 2), labels=mutlabels, gap=.1 )
dev.off()
