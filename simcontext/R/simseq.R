
#' Sample a Random Sequence
#'
#' Return a random sequence, to be used for the initial state, for instance.  
#' bases can in fact be longer patterns.
#' Note that if bases are longer than 1 character,
#'  DOES NOT actually produce a sting with frequencies basefreqs.
#'
#' @export
rinitseq <- function (seqlen,bases,basefreqs=rep(1/length(bases),length(bases))) {
    stringfun <- if ( all(bases %in% c("A","C","G","T","-") ) ) { Biostrings::DNAString } else { Biostrings::BString }
    patlen <- mean( nchar(bases) )
    rpats <- c()
    while( sum(nchar(rpats))<seqlen ) {
        rpats <- c( rpats, sample( bases, floor(1.05*seqlen/patlen), replace=TRUE, prob=basefreqs ) )
    }
    return( stringfun(substr( paste( rpats, collapse="" ), 1, seqlen )) )
}

#' Simualte From a Configuration
#'
#' Run simseq but from a config object, basically a named list
#'  note: probably want to run contextual::fill.default.config( ) on config first.
#'
#' @export
simseq.config <- function (seqlen,tlen,config,...) {
    do.call( simseq, c( list( seqlen=seqlen, tlen=tlen), config, list(...) ) )
}

#' Simulate Context-Dependent Evolution of a Sequence
#'
#' Simulate a random sequence and evolve it with genmatrix, wrapping mutations around as needed.
#'  record transition counts if count.trans it TRUE (e.g. for debugging)
#'
#' Computation of the generator matrix may be costly; so it can be stored for later use,
#'  and reloaded, by passing a file name in gmfile.  Note that this differs from other genmatrix'es stored,
#'  because we store a few more things, remove changes in the wings, and divide through rates to avoid overcounting.
#'
#' IF only.num.events is TRUE, then JUST set things up, and estimate how many mutation events will be needed,
#'   and return this (to check before embarking on something ridiculous).
#'
#' @param seqlen Length of the sequence in base pairs.
#' @param tlen Length of time to simulate for (may be a vector).
#' @param mutpats List of list of T-mers describing mutation patterns.
#' @param mutrates Vector of mutation rates.
#' @param selpats List of list of selection patterns.
#' @param selfactors List of vectors of relative selection coefficients.
#' @param selcoef Vector of selection coefficients.
#' @param patlen Maximum pattern length (determines size of genmatrix), safe to let it be automatically determined.
#' @param bases Character string of bases.
#' @param initseq List of initial sequences to start from, of the same length as tlen (if missing will be simulated).
#' @param basefreqs Frequencies to simulate initial sequences from (independently) if missing.
#' @param count.trans Whether to count up the number of transitions.
#' @param only.num.events Just count how many events to expect?
#' @param gmfile File where a precomputed generator matrix and associated other objects are found (NOTE: differs from genmatrix for inference!),
#'          or will be saved to if the file doesn't exist.
#' @param simplify Strip out the list if only simulating one sequence?
#' @param quiet Be quiet?
#' @param ... Additional arguments passed to contextual::makegenmatrix.
#'
#' @export
simseq <- function (seqlen, tlen, mutpats, mutrates, 
        selpats=list(), 
        selfactors=lapply(selpats,sapply,function(x)1), 
        selcoef=numeric(0), 
        patlen, bases=c("A","C","G","T"), 
        initseq, basefreqs=rep(1/length(bases),length(bases)), 
        count.trans=FALSE, only.num.events=FALSE, 
        gmfile=NULL, simplify=TRUE, quiet=FALSE, ... ) {
    # First make transition matrix
    #   Note that for a mutation pattern of length less than patlen, we will overcount.
    #    e.g. if the rate of CG -> TG is 1.5, and pattern length is 4,
    #    then we'll match on CG.., .CG., and ..CG, so the transition rates of e.g. CGTT -> TGTT should be 1.5/(4-2+1) = 0.5
    #    To fix this, we rescale mutrates below.
    stringsetfun <- if ( all(bases %in% c("A","C","G","T","-") ) ) { Biostrings::DNAStringSet } else { Biostrings::BStringSet }
    # determine size of padding window
    sellen <- if (length(selpats)>0) { max( sapply( unlist( selpats ), nchar ) ) } else { 0 }
    mutpatlens <- unlist( lapply( mutpats, function (x) unique(nchar(unlist(x))) ) )
    if (! all( sapply(mutpatlens,length)==1 ) ) { stop("need each list in mutpats to have patterns of the same length") }
    mutlen <- max(mutpatlens)
    if (missing(patlen)) { patlen <- mutlen }
    if (patlen<mutlen) { stop("patlen too short") }
    pad.patlen <- patlen+2*max(0,(sellen-1))
    if (!is.null(gmfile) && file.exists(gmfile)) {  # load "generator matrix" -- note this is NOT a "genmatrix" object.
        gm.obj <- load(gmfile)
        if (!all(c("genmatrix","genmatrix.j","mutrates","patterns","gm.params") %in% gm.obj)) {
            stop(paste("Didn't find expected objects in presaved generator matrix file",gmfile))
        }
        # check to make sure this genmatrix was constructed with the right structures
        if ( !isTRUE(all.equal( gm.params, list(bases=bases,mutpats=mutpats,selpats=selpats,selfactors=selfactors) ) ) ) {
            stop(paste("Parameters in", gmfile, "don't match given ones."))
        }
    } else {
        # construct generator matrix for (sellen-1,patlen,sellen-1) but with outer padding not changing
        full.genmatrix <- contextual::makegenmatrix( mutpats, selpats, patlen=pad.patlen, boundary="none", bases=bases, selfactors=selfactors, ... )
        mutrates <- mutrates / (patlen-mutpatlens+1)  # avoid overcounting (see above)
        full.genmatrix@x <- contextual::update_x(full.genmatrix,mutrates,selcoef,...)
        patterns <- rownames(full.genmatrix)
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
        if (!is.null(gmfile)) {
            cat("Saving to ", gmfile, " -- remember this differs from the model's generator matrix.\n")
            gm.params <- list(bases=bases,mutpats=mutpats,selpats=selpats,selfactors=selfactors) 
            save( genmatrix, genmatrix.j, mutrates, patterns, gm.params, file=gmfile )
        }
    }
    stopifnot( all( genmatrix@x >= 0 ) )  # assume this comes WITHOUT the diagonal
    patstrings <- stringsetfun( patterns )
    maxrate <- max( Matrix::rowSums(genmatrix) )
    diags <- maxrate - Matrix::rowSums(genmatrix)
    ####
    # number and locations of possible changes, ordered by time they occur at
    mean.n.events <- maxrate * tlen * seqlen
    if (!quiet) { cat("simseq: Total number of potential substitutions around ", paste(unique(range(mean.n.events)),collapse=" to "), "\n") }
    if (only.num.events) {  # don't actually simulate them, just say how many
        return(mean.n.events)
    }
    ###
    # if tlen is a vector, return a list of sequences
    seqlen <- rep_len(seqlen,length(tlen))
    if (missing(initseq)) { initseq <- lapply(seqlen, rinitseq, bases, basefreqs) }
    output.list <- lapply( seq_along(tlen), function (k.tlen) {
        n.events <- rpois(1,lambda=maxrate*tlen[k.tlen]*seqlen[k.tlen])
        loc.events <- sample(seqlen[k.tlen],n.events,replace=TRUE)
        wrap.events <- (loc.events+pad.patlen>seqlen[k.tlen]+1) #events wrapping around the end
        # count transitions, for debugging
        ntrans <- if (count.trans) { data.frame( i=factor(rep(NA,n.events),levels=rownames(genmatrix)), j=factor(rep(NA,n.events),levels=colnames(genmatrix)), loc=loc.events ) } else { NULL }
        # initial and final sequences
        finalseq <- initseq[[k.tlen]]
        for (k in seq_along(loc.events)) {
            subseq <- if (wrap.events[k]) { # from string (cyclical)
                    Biostrings::xscat( Biostrings::subseq( finalseq, loc.events[k], seqlen[k.tlen]), Biostrings::subseq(finalseq, 1, loc.events[k]+pad.patlen-seqlen[k.tlen]-1) )
                } else {
                    Biostrings::subseq(finalseq, loc.events[k], loc.events[k]+pad.patlen-1 )
                }
            msubseq <- Biostrings::match( subseq, patstrings )  # which row for this string?
            isubseq <- which( (genmatrix@i + 1) == msubseq ) # transitions are these entries of genmatrix
            if (length(isubseq)>0) {   # do nothing if this pattern doesn't mutate
                jsubseq <- c( msubseq, (genmatrix.j[isubseq]+1) ) # indices of possible replacement patterns (self is first)
                psubseq <- c( diags[msubseq], genmatrix@x[isubseq] )  # probabilities of choosing these
                replind <- sample( jsubseq, 1, prob=psubseq ) # choose this replacement string (index)
                replstr <- patstrings[[replind]] # what is the replacement string
                if (count.trans) {  ntrans$i[k] <- patterns[msubseq]; ntrans$j[k] <- patterns[replind] } # record this (is a factor so not storing actual strings)
                if (wrap.events[k]) { # put this back in (cyclical)
                    Biostrings::subseq( finalseq, loc.events[k], seqlen[k.tlen] ) <- Biostrings::subseq( replstr, 1, seqlen[k.tlen]-loc.events[k]+1 )
                    Biostrings::subseq( finalseq, 1, loc.events[k]+pad.patlen-seqlen[k.tlen]-1 ) <- Biostrings::subseq( replstr, seqlen[k.tlen]-loc.events[k]+2, pad.patlen )
                } else {
                    Biostrings::subseq( finalseq, loc.events[k], loc.events[k]+pad.patlen-1 ) <- replstr
                }
            }
            # sanity check:
            # if (nchar(finalseq) != nchar(initseq)) { browser() }
        }
        output <- list( initseq=initseq[[k.tlen]], finalseq=finalseq, maxrate=maxrate, ntrans=ntrans, mutpats=mutpats, 
                       selpats=selpats, mutrates=mutrates, selcoef=selcoef, tlen=tlen[k.tlen], seqlen=seqlen[k.tlen], bases=bases, fixfn.params=list(...) )
        class(output) <- "simseq"
        output
    })
    # should not do this horrible hack-ey conditional return...
    return( if (length(output.list)==1 && simplify) { output.list[[1]] } else { output.list } )
}

#' Simulate on a Tree
#'
#' @return A list of the simulated sequences in the same order as the tips,nodes of the tree
#'
#' @export
simseq.tree <- function (seqlen,config,...) {
    simseqs <- lapply(contextual::nodenames(config$tree),function(e)NULL)
    more.args <- list(...)
    simseqs[[contextual::rootname(config$tree)]] <- if ("initseq" %in% names(more.args)) {
        list(finalseq=more.args[["initseq"]])
    } else {
        list(finalseq=rinitseq(seqlen,config$bases,basefreqs=config$initfreqs))
    }
    more.args[["initseq"]] <- NULL # remove it now if it's there
    for (k in 1:nrow(config$tree$edge)) {
        # edge.pair is (from, to)
        edge.pair <- config$tree$edge[k,]
        modconfig <- config[[ contextual::nodenames(config$tree)[ edge.pair[2] ] ]]
        # can't help but allow this
        nn <- 1; while (is.character(modconfig)) { modconfig <- config[[ modconfig ]]; nn <- nn+1; if (nn>100){stop("Circular model reference.")} }
        # simulate, with config for corresponding edge
        simseqs[[ edge.pair[2] ]] <- do.call( simseq, c( 
                list( initseq=list(simseqs[[ edge.pair[1] ]]$finalseq), 
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


