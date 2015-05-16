# Simulate a sequence, with context-dependent mutation rates
suppressMessages( require(Biostrings) ) # these are NOISY.
suppressMessages( require(IRanges) ) # these are NOISY.

rinitseq <- function (seqlen,bases,basefreqs=rep(1/length(bases),length(bases))) {
    # return a random sequence.  bases can in fact be longer patterns.
    # note that if bases are longer than 1 character,
    #  DOES NOT actually produce a sting with frequencies basefreqs.
    stringfun <- if ( all(bases %in% c("A","C","G","T","-") ) ) { DNAString } else { BString }
    patlen <- mean( nchar(bases) )
    rpats <- c()
    while( sum(nchar(rpats))<seqlen ) {
        rpats <- c( rpats, sample( bases, floor(1.05*seqlen/patlen), replace=TRUE, prob=basefreqs ) )
    }
    return( stringfun(substr( paste( rpats, collapse="" ), 1, seqlen )) )
}

simseq <- function (seqlen, tlen, mutpats, mutrates, selpats=list(), selfactors=lapply(selpats,sapply,function(x)1), selcoef=numeric(0), patlen, bases=c("A","C","G","T"), count.trans=FALSE, initseq, basefreqs=rep(1/length(bases),length(bases)), ... ) {
    # Simulate a random sequence and evolve it with genmatrix, wrapping mutations around as needed.
    #  record transition counts if count.trans it TRUE (e.g. for debugging)
    # First make transition matrix
    #   Note that for a mutation pattern of length less than patlen, we will overcount.
    #    e.g. if the rate of CG -> TG is 1.5, and pattern length is 4,
    #    then we'll match on CG.., .CG., and ..CG, so the transition rates of e.g. CGTT -> TGTT should be 1.5/(4-2+1) = 0.5
    #    To fix this, we rescale mutrates below.
    stringsetfun <- if ( all(bases %in% c("A","C","G","T","-") ) ) { DNAStringSet } else { BStringSet }
    # determine size of padding window
    sellen <- if (length(selpats)>0) { max( sapply( unlist( selpats ), nchar ) ) } else { 0 }
    mutpatlens <- unlist( lapply( mutpats, function (x) unique(nchar(unlist(x))) ) )
    if (! all( sapply(mutpatlens,length)==1 ) ) { stop("need each list in mutpats to have patterns of the same length") }
    mutlen <- max(mutpatlens)
    if (missing(patlen)) { patlen <- mutlen }
    if (patlen<mutlen) { stop("patlen too short") }
    pad.patlen <- patlen+2*max(0,(sellen-1))
    # construct generator matrix for (sellen-1,patlen,sellen-1) but with outer padding not changing
    full.genmatrix <- makegenmatrix( mutpats, selpats, patlen=pad.patlen, boundary="none", bases=bases, selfactors=selfactors, ... )
    mutrates <- mutrates / (patlen-mutpatlens+1)  # avoid overcounting (see above)
    full.genmatrix@x <- update(full.genmatrix,mutrates,selcoef,...)
    patterns <- rownames(full.genmatrix)
    patstrings <- stringsetfun( patterns )
    # max mutation rate
    stopifnot( inherits(full.genmatrix,"dgCMatrix") )  # we access @i, @p, and @x directly...
    full.genmatrix.j <- rep( 0:(ncol(full.genmatrix)-1), times=diff(full.genmatrix@p) )  # column indices corresp to @i
    # remove entries that change in outer wings
    unch <- ( ( substr( patterns[full.genmatrix@i+1L], 1, sellen-1 ) ==  substr( patterns[full.genmatrix.j+1L], 1, sellen-1 ) )
        & ( substr( patterns[full.genmatrix@i+1L], sellen+patlen, patlen+2*(sellen-1) ) ==  substr( patterns[full.genmatrix.j+1L], sellen+patlen, patlen+2*(sellen-1) ) ) )
    genmatrix <- new( "dgCMatrix",
            i=full.genmatrix@i[unch],
            x=full.genmatrix@x[unch],
            p=sapply(0:ncol(full.genmatrix), function(k) sum(full.genmatrix.j[unch]<k)),
            Dim=dim(full.genmatrix), Dimnames=dimnames(full.genmatrix)
        )
    genmatrix.j <- full.genmatrix.j[unch]
    stopifnot( all(genmatrix>=0) )  # assume this comes WITHOUT the diagonal
    maxrate <- max( rowSums(genmatrix) )
    diags <- maxrate - rowSums(genmatrix)
    ####
    # number and locations of possible changes, ordered by time they occur at
    n.events <- rpois(1,lambda=maxrate*tlen*seqlen)
    loc.events <- sample(seqlen,n.events,replace=TRUE)
    wrap.events <- (loc.events+pad.patlen>seqlen+1) #events wrapping around the end
    # count transitions, for debugging
    ntrans <- if (count.trans) { data.frame( i=factor(rep(NA,n.events),levels=rownames(genmatrix)), j=factor(rep(NA,n.events),levels=colnames(genmatrix)), loc=loc.events ) } else { NULL }
    # initial and final sequences
    if (missing(initseq)) { initseq <- rinitseq(seqlen,bases,basefreqs) }
    finalseq <- initseq
    for (k in seq_along(loc.events)) {
        subseq <- if (wrap.events[k]) { # from string (cyclical)
                xscat( subseq( finalseq, loc.events[k], seqlen), subseq(finalseq, 1, loc.events[k]+pad.patlen-seqlen-1) )
            } else {
                subseq(finalseq, loc.events[k], loc.events[k]+pad.patlen-1 )
            }
        msubseq <- match( subseq, patstrings )  # which row for this string?
        isubseq <- which( (genmatrix@i + 1) == msubseq ) # transitions are these entries of genmatrix
        if (length(isubseq)>0) {   # do nothing if this pattern doesn't mutate
            jsubseq <- c( msubseq, (genmatrix.j[isubseq]+1) ) # indices of possible replacement patterns (self is first)
            psubseq <- c( diags[msubseq], genmatrix@x[isubseq] )  # probabilities of choosing these
            replind <- sample( jsubseq, 1, prob=psubseq ) # choose this replacement string (index)
            replstr <- patstrings[[replind]] # what is the replacement string
            if (count.trans) {  ntrans$i[k] <- patterns[msubseq]; ntrans$j[k] <- patterns[replind] } # record this (is a factor so not storing actual strings)
            if (wrap.events[k]) { # put this back in (cyclical)
                subseq( finalseq, loc.events[k], seqlen ) <- subseq( replstr, 1, seqlen-loc.events[k]+1 )
                subseq( finalseq, 1, loc.events[k]+pad.patlen-seqlen-1 ) <- subseq( replstr, seqlen-loc.events[k]+2, pad.patlen )
            } else {
                subseq( finalseq, loc.events[k], loc.events[k]+pad.patlen-1 ) <- replstr
            }
        }
        # sanity check:
        # if (nchar(finalseq) != nchar(initseq)) { browser() }
    }
    output <- list( initseq=initseq, finalseq=finalseq, maxrate=maxrate, ntrans=ntrans, mutpats=mutpats, selpats=selpats, mutrates=mutrates, selcoef=selcoef, tlen=tlen, seqlen=seqlen, bases=bases, fixfn.params=list(...) )
    class(output) <- "simseq"
    return(output)
}

simseq.tree <- function (seqlen,config,...) {
    # return a list of the simulated sequences in the same order as the tips,nodes of the tree
    simseqs <- lapply(nodenames(config$tree),function(e)NULL)
    more.args <- list(...)
    simseqs[[rootname(config$tree)]] <- if ("initseq" %in% names(more.args)) {
        list(finalseq=more.args[["initseq"]])
    }else {
        list(finalseq=rinitseq(seqlen,config$bases,basefreqs=config$initfreqs))
    }
    more.args[["initseq"]] <- NULL # remove it now if it's there
    for (k in 1:nrow(config$tree$edge)) {
        # edge.pair is (from, to)
        edge.pair <- config$tree$edge[k,]
        modconfig <- config[[ nodenames(config$tree)[ edge.pair[2] ] ]]
        # can't help but allow this
        nn <- 1; while (is.character(modconfig)) { modconfig <- config[[ modconfig ]]; nn <- nn+1; if (nn>100){stop("Circular model reference.")} }
        # simulate, with config for corresponding edge
        simseqs[[ edge.pair[2] ]] <- do.call( simseq, c( 
                list( initseq=simseqs[[ edge.pair[1] ]]$finalseq, 
                    tlen=config$tree$edge.length[k],
                    seqlen=seqlen, 
                    mutpats=modconfig$mutpats, mutrates=modconfig$mutrates, 
                    selpats=modconfig$selpats, selfactors=modconfig$selfactors, 
                    selcoef=modconfig$selcoef, bases=config$bases, fixfn=modconfig$fixfn ), 
                more.args,
                modconfig$fixfn.params ) )
    }
    return(simseqs)
}

countpats <- function (patterns, simseq ) {
    # count number of times each pattern appears in character string simseq, cyclic-ized.
    patlen <- nchar(patterns[1])
    seqlen <- nchar(simseq)
    simseq <- xscat( simseq, subseq(simseq,1,patlen-1) ) # cyclic-ize
    return( sapply( patterns, countPattern, simseq ) )
}

counttrans.list <- function (lpatterns, seqlist=lapply(simseqs,"[[","finalseq"), simseqs, leftwin, shift=0, cyclic=FALSE, bases=sort(unique(unlist(strsplit(lpatterns[[1]],""))))) {
    # count number of times each of seqlist matches corresponding lpatterns
    #   optionally, cyclical
    # if shift is nonzero, return a list  with the counts in each of (shift) frames
    patlen <- max( sapply( lapply( lpatterns, "[", 1 ), nchar ) )
    seqlen <- unique( sapply(seqlist, nchar) )
    stopifnot(length(seqlen)==1)
    if (length(leftwin)<length(seqlist)) { shiftwin <- c(0,rep(leftwin,length.out=length(seqlist)-1)) } else { shiftwin <- leftwin }
    # cyclic-ize
    if (cyclic) { seqlist <- lapply( seqlist, function (x) { xscat( x, subseq(x,1,patlen-1) ) } ) }
    # Ok, count occurrences.  uses bioconductor stuff.
    lmatches <- lapply( seq_along(seqlist), function (k) {
            lapply( lpatterns[[k]], matchPattern, seqlist[[k]] )
        } )
    npats <- sapply(lmatches,length)
    counts <- lapply( 1:max(1,shift), function (k) array( 0, dim=npats, dimnames=lpatterns ) )
    ii <- matrix( rep(1,length(npats)), nrow=1 )
    while( ii[1] <= npats[1] ) {
        xycounts <- start(lmatches[[1]][[ii[1]]]) - shiftwin[1]
        for (j in seq_along(ii)[-1]) {
            stopifnot( j %in% seq_along(lmatches) & ii[j] %in% seq_along(lmatches[[j]]) )
            xycounts <- intersect( xycounts, (-shiftwin[j]) + start(lmatches[[j]][[ii[j]]]) )
        }
        if (length(xycounts)>0) {
            for (k in 0:max(0,shift-1)) {
                    counts[[k+1]][ii] <- sum( ( ( xycounts+k ) %% max(shift,1) ) == 0 )
            }
        }
        ii[length(ii)] <- ii[length(ii)]+1
        for (j in rev(seq_along(ii)[-1])) { ii[j-1] <- ii[j-1] + ( ii[j] > npats[j] ) }
        ii[-1] <- 1 + ( (ii[-1]-1) %% npats[-1] )
    }
    colpatterns <- do.call( expand.grid, lpatterns[-1] )
    colnames(colpatterns) <- names(seqlist)[-1]
    counts <- lapply( counts, function (x) {
            dim(x) <- c(dim(x)[1],prod(dim(x)[-1]))
            rownames(x) <- lpatterns[[1]]
            colnames(x) <- apply(colpatterns,1,paste,collapse='.')
            new("tuplecounts", 
                leftwin=leftwin,
                counts=Matrix(x),
                bases=bases,
                rowtaxon=names(seqlist)[1],
                colpatterns=colpatterns
                ) 
        } )
    if (shift==0) { counts <- counts[[1]] }
    return(counts)
}

counttrans <- function (ipatterns, fpatterns, initseq=simseqs[["initseq"]], finalseq=simseqs[["finalseq"]], simseqs, leftwin=0, shift=0, cyclic=FALSE, bases=sort(unique(unlist(strsplit(ipatterns,""))))) {
    # count number of times initseq matches ipatterns while finalseq matches fpatterns
    #   optionally, cyclical
    # if shift is nonzero, return a list  with the counts in each of (shift) frames
    patlen <- nchar(ipatterns[1])
    seqlen <- nchar(initseq)
    stopifnot( seqlen == nchar(finalseq) )
    # cyclic-ize
    if (cyclic) {
        initseq <- xscat( initseq, subseq(initseq,1,patlen-1) )
        finalseq <- xscat( finalseq, subseq(finalseq,1,patlen-1) )
    }
    # Ok, count occurrences.  uses bioconductor stuff.
    initmatches <- lapply( ipatterns, matchPattern, initseq )
    finalmatches <- lapply( fpatterns, matchPattern, finalseq )
    counts <- lapply( 1:max(1,shift), function (k) Matrix( 0, nrow=length(initmatches), ncol=length(finalmatches), dimnames=list(ipatterns,fpatterns) ) )
    for (x in seq_along(initmatches))  {
        for (y in seq_along(finalmatches)) {
            xycounts <- intersect(start(initmatches[[x]]),(-leftwin)+start(finalmatches[[y]]))
            if (length(xycounts)>0) {
                for (k in 0:max(0,shift-1)) {
                    counts[[k+1]][x,y] <- sum( ( ( xycounts+k ) %% max(shift,1) ) == 0 )
                }
            }
        }
    }
    counts <- lapply( counts, function (x) new("tuplecounts",counts=x,leftwin=leftwin,bases=bases,rowtaxon="long",colpatterns=data.frame(short=colnames(x))) )
    if (shift==0) { counts <- counts[[1]] }
    return(counts)
}

changepos <- function (mutpats) {
    # return list of indices corresponding to site the mutation pattern changes
    lapply( mutpats, lapply, function (x) { which(do.call("!=",strsplit(x,"")) ) } )
}

show.simseq <- function (x,printit=FALSE,maxchar=min(nchar(x$initseq),200),latex=FALSE) {
    # display the output of simulation by lowercasing or replacing with "." bases that didn't change
    x$initseq <- subseq(x$initseq,1,maxchar)
    fseq <- BString( subseq(x$finalseq,1,maxchar) )
    if (!is.null(x$ntrans)) {
        x$ntrans <- subset(x$ntrans,loc<=maxchar)
        patlen <- nchar( levels( x$ntrans$i )[1] )
        events <- unique(x$ntrans$loc)
        eventlocs <- (seq_len(nchar(fseq)) %in% outer(events,0:(patlen-1),"+"))
        if (length(events)>0) {
            whichchanged <- sapply( events, function (k) { kev <- (x$ntrans$loc==k); any( x$ntrans$i[kev] != x$ntrans$j[kev] ) } )
            changedlocs <- ( seq_len(nchar(fseq)) %in% outer( events[whichchanged], 0:(patlen-1), "+" ) )
            fseq[eventlocs & !changedlocs] <- BString( tolower(fseq[eventlocs & !changedlocs]) )
        }
        fseq[!eventlocs] <- "."
    }
    outstrings <- c(
            paste(as.character(x$initseq),"\n",sep=''),
            paste(as.character(fseq),"\n",sep='')
        )
    if (latex) {
        outstrings <- gsub("&\n&$", "\\\\\\\\\n", gsub( "^&", "", gsub("\\.","\\$\\\\centerdot\\$", gsub("","&", outstrings ) ) ) )
        outstrings <- paste( c(
                paste( "\\begin{center} \\setlength{\\tabcolsep}{0pt} \\begin{tabular}{",
                    paste(rep('c',nchar(outstrings[[1]])),collapse=''), "}\n", sep='' ),
                outstrings,
                "\\end{tabular} \\end{center} \n"
            ), collapse="" )
    }
    if (printit | latex) { lapply( outstrings, cat ) }
    return(invisible(outstrings))
}

print.simseq <- function (x,...) {
    invisible( show.simseq(x,printit=TRUE,maxchar=min(nchar(x$initseq),150)) )
}
