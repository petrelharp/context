# Simulate a sequence, with context-dependent mutation rates
require(Biostrings)
require(IRanges)

rinitseq <- function (seqlen,bases,basefreqs=rep(1/length(bases),length(bases))) {
    # return a random sequence.  bases can in fact be longer patterns.
    # note that if bases are longer than 1 character, 
    #  DOES NOT actually produce a sting with frequencies basefreqs.
    patlen <- mean( nchar(bases) )
    rpats <- c()
    while( sum(nchar(rpats))<seqlen ) {
        rpats <- c( rpats, sample( bases, floor(1.05*seqlen/patlen), replace=TRUE, prob=basefreqs ) )
    }
    return( DNAString(substr( paste( rpats, collapse="" ), 1, seqlen )) )
}

simseq <- function (seqlen, tlen, patlen, mutpats, mutrates, selpats, selcoef, Ne, bases=c("A","C","G","T"), count.trans=TRUE, initseq, basefreqs=rep(1/length(bases),length(bases)), ... ) {
    # simulate a random sequence and evolve it with genmatrix.
    #  record transition counts if count.trans it TRUE (e.g. for debugging)
    # First make transition matrix
    #   Note that for a mutation pattern of length less than patlen, we will overcount.
    #    e.g. if the rate of CG -> TG is 1.5, and pattern length is 4,
    #    then we'll match on CG.., .CG., and ..CG, so the transition rates of e.g. CGTT -> TGTT should be 1.5/(4-2+1) = 0.5
    #   So, rescale mutrates:
    genmatrix <- makegenmatrix( mutpats, selpats, patlen=patlen )
    mutpatlens <- lapply( mutpats, function (x) unique(nchar(unlist(x))) )
    if (! all( sapply(mutpatlens,length)==1 ) ) { stop("need each list in mutpats to have patterns of the same length") }
    mutrates <- mutrates / (patlen-unlist(mutpatlens)+1)  # avoid multiple counting
    genmatrix@x <- update(genmatrix,mutrates,selcoef,Ne)
    patterns <- rownames(genmatrix)
    patstrings <- DNAStringSet( patterns )
    # max mutation rate
    stopifnot( inherits(genmatrix,"dgTMatrix") )  # we access @i, @j, and @x directly...
    stopifnot( all(genmatrix>=0) )  # assume this comes WITHOUT the diagonal
    maxrate <- max( rowSums(genmatrix) )
    diags <- maxrate - rowSums(genmatrix)
    ####
    # number and locations of possible changes, ordered by time they occur at
    n.events <- rpois(1,lambda=maxrate*tlen*seqlen)
    loc.events <- sample(seqlen,n.events,replace=TRUE)
    wrap.events <- (loc.events+patlen>seqlen+1) #events wrapping around the end
    # count transitions, for debugging
    ntrans <- if (count.trans) { data.frame( i=factor(rep(NA,n.events),levels=rownames(genmatrix)), j=factor(rep(NA,n.events),levels=colnames(genmatrix)), loc=loc.events ) } else { NULL }
    # initial and final sequences
    if (missing(initseq)) { initseq <- rinitseq(seqlen,bases,basefreqs) }
    finalseq <- initseq
    for (k in seq_along(loc.events)) {
        subseq <- if (wrap.events[k]) { # from string (cyclical)
                xscat( subseq( finalseq, loc.events[k], seqlen), subseq(finalseq, 1, loc.events[k]+patlen-seqlen-1) )
            } else {  
                subseq(finalseq, loc.events[k], loc.events[k]+patlen-1 ) 
            }
        msubseq <- match( subseq, patstrings )  # which row for this string?
        isubseq <- which( (genmatrix@i + 1) == msubseq ) # transitions are these entries of genmatrix
        jsubseq <- c( msubseq, (genmatrix@j[isubseq]+1) ) # indices of possible replacement patterns (self is first)
        psubseq <- c( diags[msubseq], genmatrix@x[isubseq] )  # probabilities of choosing these
        replind <- sample( jsubseq, 1, prob=psubseq ) # choose this replacement string (index)
        replstr <- patstrings[[replind]] # what is the replacement string
        if (count.trans) {  ntrans$i[k] <- patterns[msubseq]; ntrans$j[k] <- patterns[replind] } # record this (is a factor so not storing actual strings)
        if (wrap.events[k]) { # put this back in (cyclical)
            subseq( finalseq, loc.events[k], seqlen) <- subseq(replstr,1,seqlen-loc.events[k]+1)
            subseq(finalseq, 1, loc.events[k]+patlen-seqlen) <- subseq(replstr,seqlen-loc.events[k]+2,patlen)
        } else {
            subseq( finalseq, loc.events[k], loc.events[k]+patlen-1 ) <- replstr
        }
    }
    output <-  list( initseq=initseq, finalseq=finalseq, maxrate=maxrate, ntrans=ntrans ) 
    class(output) <- "simseq"
    return(output)
}

countpats <- function (patterns, simseq ) {
    # count number of times each pattern appears in character string simseq, cyclic-ized.
    patlen <- nchar(patterns[1])
    seqlen <- nchar(simseq)
    simseq <- xscat( simseq, subseq(simseq,1,patlen-1) ) # cyclic-ize
    return( sapply( patterns, countPattern, simseq ) )
}

counttrans <- function (ipatterns, fpatterns, initseq=simseqs[["initseq"]], finalseq=simseqs[["finalseq"]], simseqs, lwin=0) {
    # count number of times initseq matches ipatterns while finalseq matches fpatterns
    # again, cyclical
    patlen <- nchar(ipatterns[1])
    seqlen <- nchar(initseq)
    stopifnot( seqlen == nchar(finalseq) )
    # cyclic-ize
    initseq <- xscat( initseq, subseq(initseq,1,patlen-1) )
    finalseq <- xscat( finalseq, subseq(finalseq,1,patlen-1) )
    # Ok, count occurrences.  
    initmatches <- lapply( ipatterns, matchPattern, initseq )
    finalmatches <- lapply( fpatterns, matchPattern, finalseq )
    counts <- Matrix( 0, nrow=length(initmatches), ncol=length(finalmatches), dimnames=list(ipatterns,fpatterns) )
    for (x in seq_along(initmatches))  {
        for (y in seq_along(finalmatches)) {
            counts[x,y] <- length(intersect(start(initmatches[[x]]),(-lwin)+start(finalmatches[[y]])))
        } 
    }
    return(counts)
}

show.simseq <- function (x,printit=FALSE,maxchar=1e8) { 
    # display the output of simulation by lowercasing or replacing with "." bases that didn't change
    x$initseq <- subseq(x$initseq,1,maxchar)
    fseq <- BString( subseq(x$finalseq,1,maxchar) )
    x$ntrans <- subset(x$ntrans,loc<=maxchar)
    patlen <- nchar( levels( x$ntrans$i )[1] )
    events <- unique(x$ntrans$loc)
    eventlocs <- (seq_len(nchar(fseq)) %in% outer(events,0:(patlen-1),"+"))
    if (length(events)>0) {
        whichchanged <- sapply( events, function (k) { kev <- (x$ntrans$loc==k); any( x$ntrans$i[kev] != x$ntrans$j[kev] ) } )
        changedlocs <- ( seq_len(nchar(fseq)) %in% outer( events[whichchanged], 0:(patlen-1), "+" ) )
        fseq[eventlocs & !changedlocs] <- tolower(fseq[eventlocs & !changedlocs])
    }
    fseq[!eventlocs] <- "."
    if (printit) {
        cat(paste(as.character(x$initseq),"\n",sep=''))
        cat(paste(as.character(fseq),"\n",sep=''))
    }
    return(invisible(fseq))
} 

print.simseq <- function (x,...) {
    invisible( show.simseq(x,printit=TRUE,maxchar=200) )
}
