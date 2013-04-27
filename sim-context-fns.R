# Simulate a sequence, with context-dependent mutation rates

simseq <- function (seqlen, tlen, genmatrix, bases, count.trans=TRUE, basefreqs=rep(1/length(bases),length(bases)), ... ) {
    # simulate a random sequence and evolve it with genmatrix.
    #  record transition counts if count.trans it TRUE (e.g. for debugging)
    patterns <- rownames(genmatrix)
    patlen <- nchar(patterns[1])
    # max mutation rate
    maxrate <- max( rowSums(genmatrix) )
    ####
    # number and locations of possible changes, ordered by time they occur at
    n.events <- rpois(1,lambda=maxrate*tlen*(seqlen-patlen+1))
    loc.events <- sample(seqlen-patlen+1,n.events,replace=TRUE)
    # count transitions, for debugging
    ntrans <- if (count.trans) { Matrix(0,nrow=length(patterns),ncol=length(patterns)) } else { NULL }
    # initial and final sequences
    finalseq <- initseq <- paste( sample(bases,seqlen,replace=TRUE,prob=basefreqs), collapse="" )
    for (k in loc.events) {
        subseq <- substr(finalseq, k, k+patlen-1 )
        msubseq <- match( subseq, patterns )
        isubseq <- which( (genmatrix@i + 1) == msubseq )
        # indices of possible replacement patterns
        jsubseq <- (genmatrix@j+1)[isubseq]
        # probabilities of choosing these
        psubseq <- genmatrix@x[isubseq]/maxrate
        replind <- sample( c(msubseq,1+seq_along(patterns)[jsubseq]), 1, prob=c(max(0,1-sum(psubseq)),psubseq) )
        if (count.trans) { ntrans[msubseq,jsubseq] <- ntrans[msubseq,jsubseq] + 1 }
        replstr <- c(subseq,patterns)[replind]
        substr( finalseq, k, k+patlen-1 ) <- replace
    }
    return( list( initseq=initseq, finalseq=finalseq, maxrate=maxrate, ntrans=ntrans ) )
}

gettransmatrix <- function (mutpats, mutrates, selpats, selcoef, Ne, tlen, win, lwin=0, rwin=0, ... ) {
    # get reduced transition matrix: given (lwin, win, rwin) context probability of pattern in win
    winlen <- lwin+win+rwin
    ipatterns <- getpatterns(winlen)
    fpatterns <- getpatterns(win)

    fullgenmatrix <- makegenmatrix( mutpats, selpats, ipatterns )
    fullgenmatrix@x <- update(fullgenmatrix,mutrates,selcoef,Ne)

    transmatrix <- expm( tlen * (fullgenmatrix-Diagonal(nrow(fullgenmatrix),rowSums(fullgenmatrix))), method="Higham08" )

    subgenmatrix <- collapsepatmatrix( fullgenmatrix, lwin=lwin, rwin=rwin )
    subtransmatrix <- transmatrix %*% subgenmatrix
    rownames(subtransmatrix) <- ipatterns
    colnames(subtransmatrix) <- fpatterns
    return( subtransmatrix )
}

whichchanged <- function (ipatterns,fpatterns,lwin=0,win=nchar(ipatterns[0])) {
    # return indicator corresponding to entries of output of gettransmatrix that have changed
    return( outer( ipatterns, fpatterns, function (x,y) { substr(x,lwin+1,lwin+win)!=y } ) )
}

counttrans <- function (ipatterns, fpatterns, simseqs) {
    initseq <- simseqs[["initseq"]]
    finalseq <- simseqs[["finalseq"]]
    # count number of times initseq matches ipatterns while finalseq matches fpatterns
    seqlen <- nchar(initseq)
    stopifnot( seqlen == nchar(finalseq) )
    # Ok, count occurrences.  Note needs perl "lookahead" to count overlapping patterns.
    #   (see http://stackoverflow.com/questions/7878992/finding-the-indexes-of-multiple-overlapping-matching-substrings)
    initmatches <- sapply( lapply( ipatterns, function (p) gregexpr(paste("(?=",p,")",sep=''),initseq,perl=TRUE) ), "[[", 1 )
    finalmatches <- sapply( lapply( fpatterns, function (p) gregexpr(paste("(?=",p,")",sep=''),finalseq,perl=TRUE) ), "[[", 1 )
    counts <- Matrix( 0, nrow=length(initmatches), ncol=length(finalmatches) )
    for (x in seq_along(initmatches))  {
        for (y in seq_along(finalmatches)) {
            counts[x,y] <- length(intersect(initmatches[[x]],(-lwin)+finalmatches[[y]]))
        } 
    }
    stopifnot(sum(counts)==seqlen-winlen+1)
    return(counts)
}
