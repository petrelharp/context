#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Run mcmc longer.\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file containing results from initial run." ),
        make_option( c("-j","--jobid"), type="character", default=formatC(1e6*runif(1),width=6,format="d",flag="0"), help="Unique job id. [default random]"),
        make_option( c("-n","--nbatches"), type="integer", default=1000, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-c","--mutstepscale"), type="numeric", default=5e-6, help="Scale of proposal steps for mutation rates in the Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-t","--timestepscale"), type="numeric", default=2e-4, help="Scale of proposal steps for time values in the Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-s","--restart"), action="store_true", default=FALSE, help="Start a whole new MCMC run?" ),
        make_option( c("-d","--boundary"), type="character", default="none", help="Boundary conditions for generator matrix. [default \"%default\"]"),
        make_option( c("-y","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases in computing generator matrix. [default \"%default\"]" ),
        make_option( c("-g","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix, or TRUE [default] to look for one. (otherwise, will compute)"),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" ),
        make_option( c("-e","--debug"), type="logical", default=FALSE, help="Debug output, e.g. dump to file on error?" )
    )
mcmcopt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(mcmcopt$infile)) { stop("No input file.  Run\n  bcells-mcmc.R -h\n for help.\n") }
attach(mcmcopt)

if (interactive()) { nbatches <- 100; blen <- 10; restart <- FALSE; gmfile <- TRUE }
if (mcmcopt$debug & !interactive()) { options(error = quote({dump.frames(to.file = TRUE); q()})) }

if (!file.exists(infile)) { stop("Cannot read file ", infile) }
basename <- gsub("-results.RData","",infile)

scriptdir <- "../"
source(paste(scriptdir,"context-inference-fns.R",sep=''),chdir=TRUE)
# source(paste(scriptdir,"sim-context-fns.R",sep=''),chdir=TRUE)
require(mcmc)

load(infile)  # has mrun and previous things (called 'datafile' in -inference.R)

longwin <- leftwin + shortwin + rightwin
if (!exists("patlen")) { patlen <- opt$patlen } # workaround

plotfile <- paste( basename ,"-plot",sep='')
mcmcdatafiles <- list.files(path=".",pattern=paste(basename,"-mcmc.*RData",sep=''),full.names=TRUE)
mcmcnum <- 1+max(c(0,as.numeric(gsub(".*-mcmc-","",gsub(".RData","",mcmcdatafiles)))),na.rm=TRUE)

if (logfile=="") { logfile <- paste(basename,"-mcmc-run-",mcmcnum,".Rout",sep='') }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    if (!interactive()) { sink(file=logcon, type="message") }
    sink(file=logcon, type="output", split=interactive()) 
}
date()
cat("basename: ", basename, "\n")
print(mcmcopt)

if (length(mcmcdatafiles)>0) { load(grep(paste("-mcmc-",mcmcnum-1,".RData",sep=''),mcmcdatafiles,fixed=TRUE,value=TRUE)) }

########

# Inference.
if (gmfile=="TRUE") { gmfile <- paste(paste("genmatrices/genmatrix",longwin,boundary,meanboundary,patlen,sep="-"),".RData",sep='') }
if (file.exists(gmfile)) {
    load(gmfile)
} else {
    stop("Can't find the generator matrix in", gmfile, " . Specify filename directly?")
}
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE, time="gamma" )

# counts loaded from infile
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

if (length(mutstepscale)==1) { mutstepscale <- rep(mutstepscale,nmuts) }
stepscale <- c(mutstepscale,timestepscale)

if (restart) {
    mrun <- metrop( mrun, initial=adhoc.ans$par, nbatch=nbatches, blen=blen, scale=stepscale )
} else {
    mrun <- metrop( mrun, nbatch=nbatches, blen=blen, scale=stepscale )
}


date()
savefile <- paste(basename,"-mcmc-",mcmcnum,".RData",sep='')
cat("saving to: ", savefile, "\n")
save( leftwin, shortwin, rightwin, patlen, lud, mrun, initcounts, mcmcopt, file=savefile )

param.names <- c( sapply(mutpats,function(x){paste(sapply(x,paste,collapse='->'),collapse='|')}), "shape" )

# scale tends to be much larger than mutation rates; scale it.
scaled.batch <- sweep(mrun$batch,2,c(rep(1,nmuts),.01),"*")

pdf(file=paste(plotfile,"-mcmc-",mcmcnum,".pdf",sep=''), width=11, height=8, pointsize=10)
matplot( scaled.batch, type='l', col=1:length(point.estimate) )
abline(h=point.estimate, col=adjustcolor(1:length(point.estimate),.5), lwd=1, lty=2)
legend("topright", col=(1:length(point.estimate)), lty=1, legend=param.names)
dev.off()

pdf(file=paste(plotfile,"-mcmc-",mcmcnum,"-pairwise.pdf",sep=''),width=11,height=8,pointsize=10)
pairs( rbind( scaled.batch, point.estimate ), col=c( rep(adjustcolor("black",.1),nrow(mrun$batch)), adjustcolor("red",1)), pch=20, cex=c( rep(.25,nrow(mrun$batch)), 2), labels=param.names, gap=.1 )
dev.off()

cat("done.\n")
date()
