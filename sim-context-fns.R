# Simulate a sequence, with context-dependent mutation rates

simseq <- function (seqlen, tlen, genmatrix, bases, count.trans=TRUE, initseq, basefreqs=rep(1/length(bases),length(bases)), ... ) {
    # simulate a random sequence and evolve it with genmatrix.
    #  record transition counts if count.trans it TRUE (e.g. for debugging)
    patterns <- rownames(genmatrix)
    patlen <- nchar(patterns[1])
    # max mutation rate
    stopifnot( inherits(genmatrix,"dgTMatrix") )  # we access @i, @j, and @x directly...
    stopifnot( all(genmatrix>=0) )  # assume this comes WITHOUT the diagonal
    maxrate <- max( rowSums(genmatrix) )
    diags <- maxrate - rowSums(genmatrix)
    ####
    # number and locations of possible changes, ordered by time they occur at
    n.events <- rpois(1,lambda=maxrate*tlen*(seqlen-patlen+1))
    loc.events <- sample(seqlen-patlen+1,n.events,replace=TRUE)
    # count transitions, for debugging
    ntrans <- if (count.trans) { data.frame( i=rep(NA,n.events), j=rep(NA,n.events), loc=loc.events ) } else { NULL }
    # initial and final sequences
    if (missing(initseq)) { initseq <- paste( sample(bases,seqlen,replace=TRUE,prob=basefreqs), collapse="" ) }
    finalseq <- initseq
    for (k in seq_along(loc.events)) {
        subseq <- substr(finalseq, loc.events[k], loc.events[k]+patlen-1 )  # from string
        msubseq <- match( subseq, patterns )  # which row for this string?
        isubseq <- which( (genmatrix@i + 1) == msubseq ) # transitions are these entries of genmatrix
        jsubseq <- c( msubseq, (genmatrix@j[isubseq]+1) ) # indices of possible replacement patterns (self is first)
        psubseq <- c( diags[msubseq], genmatrix@x[isubseq] )  # probabilities of choosing these
        replind <- sample( jsubseq, 1, prob=psubseq ) # choose this replacement string (index)
        if (count.trans) {  ntrans$i[k] <- msubseq; ntrans$j[k] <- replind } # record this
        replstr <- patterns[replind] # what is the replacement string
        substr( finalseq, loc.events[k], loc.events[k]+patlen-1 ) <- replstr  # put it in
    }
    return( list( initseq=initseq, finalseq=finalseq, maxrate=maxrate, ntrans=ntrans ) )
}


counttrans <- function (ipatterns, fpatterns, simseqs, lwin=0) {
    initseq <- simseqs[["initseq"]]
    finalseq <- simseqs[["finalseq"]]
    # count number of times initseq matches ipatterns while finalseq matches fpatterns
    seqlen <- nchar(initseq)
    stopifnot( seqlen == nchar(finalseq) )
    # Ok, count occurrences.  Note needs perl "lookahead" to count overlapping patterns.
    #   (see http://stackoverflow.com/questions/7878992/finding-the-indexes-of-multiple-overlapping-matching-substrings)
    initmatches <- lapply( ipatterns, function (p) gregexpr(paste("(?=",p,")",sep=''),initseq,perl=TRUE) )
    finalmatches <- lapply( fpatterns, function (p) gregexpr(paste("(?=",p,")",sep=''),finalseq,perl=TRUE) )
    counts <- Matrix( 0, nrow=length(initmatches), ncol=length(finalmatches) )
    for (x in seq_along(initmatches))  {
        for (y in seq_along(finalmatches)) {
            counts[x,y] <- length(intersect(initmatches[[x]][[1]],(-lwin)+finalmatches[[y]][[1]]))
        } 
    }
    return(counts)
}
