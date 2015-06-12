
scriptdir <- "../"
source(paste(scriptdir,"context-inference-fns.R",sep=''),chdir=TRUE)
source(paste(scriptdir,"sim-context-fns.R",sep=''),chdir=TRUE)

require(mcmc)

# mean of exponential prior on mutation rate
pmean <- 1
# variance of gaussian prior on selection coefficient
pvar <- 1

####
# available simulated sequences
simdir <- "ising-sims"
simfiles <- list.files(simdir,"*.RData",full.names=TRUE)
siminfo <- do.call(rbind, lapply(simfiles, function (x) {
        load(x)
        y <- do.call( data.frame, c(list(date=now,seqlen=seqlen, tlen=tlen, file=x), as.list(mutrates), as.list(selcoef), stringsAsFactors=FALSE) )
        colnames(y) <- c( colnames(y)[1:4], paste("mutrate",seq_along(mutrates),sep=''), paste("selcoef",seq_along(selcoef),sep='') )
        y
    } ) )
siminfo$meandist <- siminfo$tlen * colSums(siminfo[,grep("mutrate",colnames(siminfo)),drop=FALSE])
siminfo <- siminfo[ order(siminfo$meandist), ]

# pick one
thissim <- siminfo[1,"file"]

load(thissim)

######
# Inference.
shortwin <- 2
leftwin <- rightwin <- 2
longwin <- leftwin+shortwin+rightwin

genmatrix <- makegenmatrix( mutpats, selpats=selpats, patlen=longwin, boundary="wrap")
genmatrix@x <- update(genmatrix,mutrates,selcoef=selcoef,Ne=1)
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

counts <- list( counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[1]]$initseq, simseqs[[1]]$finalseq, leftwin=leftwin ) )
# want only patterns with leftmost possible position changed
nonoverlapping <- leftchanged(rownames(counts[[1]]),colnames(counts[[1]]),leftwin=leftwin,shortwin=shortwin)
ncounts <- counts[[1]][nonoverlapping]

nmuts <- length(mutrates)
nsel <- length(selcoef)
# move from base frequencies (what we estimate) to pattern frequencies
lud <- function (params) {
    if (any(params<0)) {
        return( -Inf )
    } else {
        # params are: mutrates*tlen, selcoef
        # this is collapsed transition matrix
        mutrates <- params[1:nmuts]
        selcoef <- params[nmuts+1:nsel]
        genmatrix@x <- update(genmatrix,mutrates=mutrates,selcoef=selcoef)
        subtransmatrix <- computetransmatrix( genmatrix, projmatrix )
        # return negative log-likelihood 
        # return( sum( counts[[1]] * log(subtransmatrix) ) + sum(pmean*mutrates) + sum(selcoef^2)/pvar )
        return( sum( ncounts * log(subtransmatrix[nonoverlapping]) ) + sum(pmean*mutrates) + sum(selcoef^2)/pvar )
    }
}


truth <- c(mutrates*tlen,selcoef)
system.time( mrun <- metrop( lud, initial=truth, nbatch=1000, blen=10, scale=1e-2 ) )

