require(mcmc)

set.seed(42)

# mean of exponential prior on mutation rate
pmean <- 1

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

win <- 2

# pick one
thissim <- siminfo[1,"file"]

load(thissim)

######
# Inference.
lwin <- rwin <- thiswin
winlen <- lwin+win+rwin

genmatrix <- makegenmatrix( mutpats, selpats=selpats, patlen=winlen, boundary="wrap")
genmatrix@x <- update(genmatrix,mutrates,selcoef=selcoef,Ne=1)
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

counts <- list(
            counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[1]]$initseq, simseqs[[1]]$finalseq, lwin=lwin )
        )
# allow nonmatches with small prob?
eps <- 0
nmuts <- length(mutrates)
nsel <- length(selcoef)
pmeans <- rep(pmean,nmuts+nsel)
# move from base frequencies (what we estimate) to pattern frequencies
lud <- function (params) {
    if (any(params<0)) {
        return( -Inf )
    } else {
        # params are: mutrates, selcoef
        # this is collapsed transition matrix
        mutrates <- params[1:nmuts]
        selcoef <- params[nmuts+1:nsel]
        genmatrix@x <- update(genmatrix,mutrates=mutrates,selcoef=selcoef)
        subtransmatrix <- computetransmatrix( genmatrix, projmatrix )
        # return negative log-likelihood 
        return( sum( counts[[1]] * log(eps+subtransmatrix) ) + sum(pmean*params) )
    }
}


system.time( mrun <- metrop( lud, initial=truth, nbatch=1000, blen=10, scale=1e-2 ) )

