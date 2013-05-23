require(mcmc)

set.seed(42)

# mean of exponential prior on mutation rate
pmean <- 1

# Inference.
lwin <- rwin <- thiswin
winlen <- lwin+win+rwin

genmatrix <- makegenmatrix( mutpats, selpats=c(), patlen=winlen, boundary="wrap")
genmatrix@x <- update(genmatrix,mutrates,selcoef=numeric(0),Ne=1)
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

counts <- list(
            counttrans( rownames(projmatrix), colnames(projmatrix), simseqs[[1]]$initseq, simseqs[[1]]$finalseq, lwin=lwin )
        )
# allow nonmatches with small prob
eps <- 1e-6

lud <- function (mutrates) {
    if (mutrates <= 0) { 
        return( -Inf )
    } else {
        # params are: Ne, branchlens[-1], mutrates, selcoef, basefreqs
        # this is collapsed transition matrix
        genmatrix@x <- update(genmatrix,mutrates,selcoef=numeric(0),Ne=1)
        subtransmatrix <- computetransmatrix( genmatrix, projmatrix )
        # return POSITIVE log-likelihood 
        return( sum( counts[[1]] * log(eps+subtransmatrix) ) + pmean * mutrates )
    }
}

system.time( mrun <- metrop( lud, initial=truth, nbatch=20, blen=10 ) )

