#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Run mcmc longer.\
"

option_list <- list(
        make_option( c("-u","--basedir"), type="character", default=NULL, help="Directory to look for input in, and write output files to." ),
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file containing simulation." ),
        make_option( c("-w","--win"), type="integer", default=1, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--lwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-k","--patlen"), type="integer", default=1, help="Include mutation rates for all tuples of this length. [default \"%default\"]" ),
        make_option( c("-n","--nbatches"), type="integer", default=1000, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-c","--stepscale"), type="numeric", default=1e-4, help="Scale of proposal steps for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-s","--restart"), action="store_true", default=FALSE, help="Start a whole new MCMC run?" ),
        make_option( c("-d","--boundary"), type="character", default="none", help="Boundary conditions for generator matrix. [default \"%default\"]"),
        make_option( c("-y","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases in computing generator matrix. [default \"%default\"]" ),
        make_option( c("-g","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix, or TRUE [default] to look for one. (otherwise, will compute)"),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$infile) & is.null(opt$basedir)) { stop("No input file.  Run\n  bcells-inference.R -h\n for help.\n") }
attach(opt)
options(error = quote({dump.frames(to.file = TRUE); q()}))

if (interactive()) { nbatches <- 100; blen <- 10; restart <- FALSE; gmfile <- TRUE }
if (!interactive()) { options(error = quote({dump.frames(to.file = TRUE); q()})) }

winlen <- lwin+win+rwin

if (substr(basedir,nchar(basedir),nchar(basedir)) %in% c("/","\\")) { basedir <- substr(basedir,1,nchar(basedir)-1) }
basename <- paste(basedir,"/win-",lwin,"-",win,"-",rwin,"-",patlen,sep='')
if (is.null(opt$infile)) { infile <- paste( basename ,"-results.RData",sep='') }
if (!file.exists(infile)) { stop("Cannot read file ", infile) }

scriptdir <- "../"
source(paste(scriptdir,"codon-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))
require(mcmc)

load(infile)
if (!file.exists(basedir)) { dir.create(basedir) }

plotfile <- paste( basename ,"-plot",sep='')
mcmcdatafiles <- list.files(path=basedir,pattern="-mcmc.*RData",full.names=TRUE)
mcmcnum <- 1+max(c(0,as.numeric(gsub(".*-mcmc-","",gsub(".RData","",mcmcdatafiles)))),na.rm=TRUE)

if (logfile=="") { logfile <- paste(basename,"-mcmc-run-",mcmcnum,".Rout",sep='') }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    if (!interactive()) { sink(file=logcon, type="message") }
    sink(file=logcon, type="output", split=interactive()) 
}

load(datafile)  # has mrun and previous things
if (length(mcmcdatafiles)>0) { load(grep(paste("-mcmc-",mcmcnum-1,".RData",sep=''),mcmcdatafiles,fixed=TRUE,value=TRUE)) }

########

# Inference.
if (gmfile=="TRUE") { gmfile <- paste(paste("genmatrices/genmatrix",winlen,boundary,meanboundary,sep="-"),".RData",sep='') }
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
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE, time="gamma" )

count.table <- read.table(infile,header=TRUE,stringsAsFactors=FALSE)
counts <- Matrix(0,nrow=nrow(genmatrix),ncol=ncol(projmatrix))
rownames(counts) <- rownames(genmatrix)
colnames(counts) <- colnames(projmatrix)
stopifnot( all( count.table$reference %in% rownames(genmatrix) ) & all(count.table$derived %in% colnames(projmatrix)) ) 
counts[cbind( match(count.table$reference,rownames(genmatrix)), match(count.table$derived,colnames(projmatrix)) )] <- count.table$count
initcounts <- rowSums(counts)

nmuts <- length(mutpats)
nsel <- length(selpats)
stopifnot(nsel==0)

point.estimate <- adhoc.ans$par

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

if (restart) {
    mrun <- metrop( mrun, initial=adhoc.ans$par, nbatch=nbatches, blen=blen, scale=stepscale )
} else {
    mrun <- metrop( mrun, nbatch=nbatches, blen=blen, scale=stepscale )
}


save( lwin, win, rwin, lud, mrun, initcounts, nov.counts, nonoverlapping, file=paste(basename,"-mcmc-",mcmcnum,".RData",sep='') )

pdf(file=paste(plotfile,"-mcmc-",mcmcnum,".pdf",sep=''),width=6, height=4, pointsize=10)
matplot( mrun$batch, type='l', col=1:length(point.estimate) )
abline(h=point.estimate, col=adjustcolor(1:length(point.estimate),.5), lwd=2)
abline(h=estimates["ans",], col=1:length(point.estimate) )
legend("topright", col=c(1:length(point.estimate),1,adjustcolor(1,.5)), lwd=c(rep(1,length(point.estimate)),1,2),legend=c(colnames(estimates)[1:length(point.estimate)],"point estimate","point.estimate"))
dev.off()

mutlabels <- paste("mut:", unlist( sapply( sapply( mutpats, lapply, paste, collapse="->" ), paste, collapse=" | " ) ) )
pdf(file=paste(plotfile,"-mcmc-",mcmcnum,"-pairwise.pdf",sep=''),width=6,height=6,pointsize=10)
pairs( rbind( mrun$batch, point.estimate ), col=c( rep(adjustcolor("black",.1),nrow(mrun$batch)), adjustcolor("red",1)), pch=20, cex=c( rep(.25,nrow(mrun$batch)), 2), labels=mutlabels, gap=.1 )
dev.off()
