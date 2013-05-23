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
# allow nonmatches with small prob?
eps <- 0

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

nmuts <- length(mutrates)
nsel <- length(selcoef)
# move from base frequencies (what we estimate) to pattern frequencies
lud <- function (params) {
    if (any(params)<0) {
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


system.time( mrun <- metrop( lud, initial=truth, nbatch=20, blen=10 ) )

