#!/usr/bin/R
require(Matrix)
require(expm)
require(stringdist)
# # find what directory this file is in
frame_files <- lapply(sys.frames(), function(x) x$ofile)
frame_files <- Filter(Negate(is.null), frame_files)
.PATH <- dirname(frame_files[[length(frame_files)]])
# source(paste(.PATH,"/expm-simple.R",sep=''))  # expAtv is faster

source(paste(.PATH,"/expAtv.R",sep=''))  # fixed upstream
source(paste(.PATH,"/gammaAtv.R",sep=''))
source(paste(.PATH,"/input-output.R",sep=''))

getpatterns <- function(patlen,bases) {
    # construct a list of all patterns of a given length
    patterns <- do.call( expand.grid, rep( list(bases), patlen ) )
    return( apply(patterns,1,paste,collapse="") )
}

mutpatchanges <- function (mutpats) {
    # return list of single-base changes for each list of patterns in mutpats
    # a `mutpat` is a 2-element string list of the form (from, to)
    lapply( mutpats, function (x) { do.call(rbind, lapply( x, function (y) {
                ysplit <- strsplit(y,'')
                differs <- do.call("!=",ysplit)
                cbind( ysplit[[1]][differs], ysplit[[2]][differs] )
            } ) ) } )
}

getmutpats <- function(bases,patlen,nchanges=1) {
    # YYY suggested name change: gen_mutpatl for generate mutpat list
    # all kmer -> kmer changes where k <= `patlen`
    # that involve at most nchanges changes
    mutpats <- list()
    patterns <- getpatterns(patlen,bases)
    for (k in 1:patlen) {
        kmers <- getpatterns(k,bases)
        mutpats <- c( mutpats,
                apply(combn(kmers,2),2,list), # make lists of rows (2 = apply over columns), giving 2-element lists of kmers
                apply(combn(kmers,2)[2:1,,drop=FALSE],2,list) # and the reverse of the 2-element lists
            )
    }
    obschanges <- sapply(mutpatchanges(mutpats),nrow)
    return( mutpats[obschanges%in%nchanges] )
}

npatterns <- function (patlen,bases) {
    return( length(bases)^patlen )
}

mutnames <- function (mutpats) {
    # Stringify a list of mutpat lists.
    return( unlist( sapply( sapply( mutpats, lapply, paste, collapse="->" ), paste, collapse="|" ) ) )
}

###
# Point estimates

divergence <- function (counts, leftwin) {
    # mean density of nucleotide differences
    require(stringdist) # http://www.inside-r.org/packages/cran/stringdist/docs/stringdist
    patlen <- nchar(colnames(counts)[1])
    nchanges <- stringdistmatrix( substr(rownames(counts),leftwin+1,leftwin+patlen), colnames(counts), method="hamming" )
    return( sum(nchanges*counts)/(sum(counts)*patlen) )
}

countmuts <- function (counts, mutpats, leftwin, ...) {
    # YYY suggest mutpatl rather than mutpats, because that's what we have
    # given a contingency table `counts` of Tmers observed in a data set, a collection of
    # mutation patterns, the leftwin, and a list of arguments to `sum`, this function
    # returns a matrix with columns indexed by the mutpats, the first row of which gives counts of how many of the counted mutations could be produced by each of mutpats
    #   and the second of which gives the number of "from" matches of the mutpats (called the "total possible").
    #
    # in other words, for each mutation pattern a -> b
    #  sum the values of counts[u,v] over choices of u,v such that:
    #   (i) 'a' matches 'u' at some position
    #   (ii) in addition, 'b' matches 'v' at the same position;
    # store this in the first row
    #
    # note that if we estimate rates by
    #    r.est <- countmuts(...)[1,]/countmuts(...)[2,]
    # then something like
    #    sum( r.est ) / 4
    # should be close to the mean density of nucleotide changes
    #    divergence(...)
    counts <- as.matrix(counts)
    # `shortwin` is the length of the inner window
    shortwin <- nchar(colnames(counts)[1])
    stopifnot( length(shortwin)>0 & shortwin>0 ) # length statement catches the case that there are no colnames for counts
    # trim off windows
    xx <- substr(rownames(counts),leftwin+1,leftwin+shortwin)
    yy <- colnames(counts)
    # Note that `observed` counts a change multiply if it can occur in different ways.
    observed <- possible <- numeric(length(mutpats)) # observed and possible are each numeric vectors
    for (j in seq_along(mutpats)) {
        mutpat <- mutpats[[j]]
        for (k in 1:(shortwin-1)) {
            mx <- sapply(mutpat, function (mp) {
                    patlen <- nchar(mp[1])
                    if ( k+patlen-1 <= shortwin ) {
                        sum( counts[ ( substr(xx,k,k+patlen-1) == mp[1] ), ( substr(yy,k,k+patlen-1) == mp[2] ) ], ... )
                    } else {
                        0
                    }
                } )
            observed[j] <- observed[j] + sum( mx, ... )
            px <- sapply(mutpat, function (mp) {
                    patlen <- nchar(mp[1])
                    if ( k+patlen-1 <= shortwin ) {
                        sum( abs(counts)[ ( substr(xx,k,k+patlen-1) == mp[1] ), ], ... )
                    } else {
                        0
                    }
                } )
            possible[j] <- possible[j] + sum( px, ... )
        }
    }
    x <- rbind(observed,possible)
    colnames(x) <- mutnames(mutpats)
    return(x)
}


###
# optimization helper functions

gradest <- function (likfun, params, eps=mean(params)/1000) {
    # estimate gradient, crudely
    gradup <- sapply( seq_along(initpar), function (k) { likfun(initpar+ifelse(seq_along(initpar)==k,eps,0)) } )
    graddown <- sapply( seq_along(initpar), function (k) { likfun(initpar-ifelse(seq_along(initpar)==k,eps,0)) } )
    return((gradup-graddown)/(2*eps))
}

###
# functions to prepare mutation and selection matrices, and build fixation functions
# In this section we start seing mutation pattern list lists. These objects make sense because the outer list
# corresponds to the mutation rates, and the inner list corresponds to a set of mutpats that all have the
# same rate.

getmutmats <- function(mutpats,patterns,boundary=c("none","wrap")) {
    # The purpose of this function is to construct a translation table between
    # mutpats and the entries of the generator matrix that can happen as a
    # consequence of those mutpats.
    # Specifically, returns i and j's indexing the generator matrix in terms of the given patterns.
    # that is, given a list of mutation patterns,
    #   which can be either pairs or lists of pairs,
    # YYY this option is kind of driving me batty. Could we have it just always be a mutpat list list?
    # return a corresponding list of two-column matrices with (1-based) indices of changes corresponding to mutation patterns
    #   i.e. if (i,j) is a row of output[[k]], then patterns[j] can be obtained from patterns[i]
    #   by performing the substitution from mutpats[[k]][1] -> mutpats[[k]][2]
    #   at some location within the string.
    #   YYY I believe this description is a little flawed, in that mutpats is a mutpat list list
    boundary <- match.arg(boundary)
    longwin <- nchar(patterns[1])
    mutpats <- lapply( mutpats, function (x) { if (is.list(x)) { x } else { list(x) } } )
    lapply( mutpats, function (y) {  # y is mutpat list
            do.call( rbind, lapply(y, function (x) {  # x is mutpat (from, to)
                patlen <- nchar(x[1])
                switch( boundary,
                    wrap={ # patterns are circular
                        wpatterns <- paste( patterns, substr(patterns,1,patlen), sep='' )
                        maxshift <- longwin
                        subsfun <- function (pat,topat,k) { wrapsubstr( pat, k, k+patlen-1 ) <- topat; return(pat) }
                    },
                    none={
                        wpatterns <- patterns
                        maxshift <- longwin-patlen+1
                        subsfun <- function (pat,topat,k) { substr( pat, k, k+patlen-1 ) <- topat; return(pat) }
                    }
                )
                do.call( rbind, lapply( 1:maxshift, function (k) {  # k is position of mutpat "from" in long pattern
                        i <- which( substr( wpatterns, k, k+patlen-1 ) == x[1] ) # which long patterns match mutpat "from" at position k
                        replstr <- subsfun( patterns[i], x[2], k ) # substitute in mutpat "to"
                        j <- match( replstr, patterns )  # indices of mutated strings
                        data.frame( i=i, j=j )
                    } ) )
        } ) )
    } )
}

getselmatches <- function (selpats, patterns, boundary=c("none","wrap"), names=FALSE) {
    # selpats can be a vector or a list of vectors,
    #   each element gets one parameter
    # selmatches[i,j] is number of times anything in selpat[[i]] matches pattern[j]
    boundary <- match.arg(boundary)
    substrfun <- switch( boundary, wrap=wrapsubstr, none=substr )
    patlen <- nchar(patterns[1])
    if (!is.list(selpats)) { selpats <- as.list(selpats) }
    selmatches <- do.call( rbind, lapply(selpats, function (y) {
            rowSums( sapply(y, function (x) {
                    maxshift <- patlen - switch( boundary, wrap=1, none=nchar(x) )
                    rowSums( sapply( 0:maxshift, function (k) {
                                x == substrfun( patterns, 1+k, k+nchar(x) )
                                # xx <- paste( c(rep(".",k), x, rep(".", patlen-regexplen(x)-k)), collapse='' )
                                # grepl( xx, patterns )
                        } ) )
                } ) )
        } ) )
    if (names) {
        rownames(selmatches) <- names(selpats)
        colnames(selmatches) <- patterns
    }
    return(selmatches)
}


popgen.fixfn <- function (ds,Ne,...) {
    # total influx of fixation given selection coefficient (s[to] - s[from]) difference ds
    if (length(ds)==0) { 1 } else { ifelse( ds==0, 1, Ne*expm1(-2*ds)/expm1(-2*Ne*ds) ) }
}

ising.fixfn <- function (ds,...) { 1/(1+exp(-ds)) }

null.fixfn <- function (...) { 1 }


###
# genmatrix code
# genmatrix gives the instantaneous rate for going from patterns x -> y
# genmatrix extends a sparse matrix class, carrying along more information.
# it is laid out in classical Markov fashion, with the rows indexing the "from" states
# (headpats) and columns indexing the "to" states (tailpats)
# genmatrix extends the sparse matrix class, carrying along more information.

check.context <- function (cont) {
    return(
        all( rownames(cont@genmatrix) == cont@headpats ) &
        all( rownames(cont@projmatrix) == cont@headpats ) &
        all( colnames(cont@projmatrix) == cont@tailpats ) &
        all( nchar(cont@headpats) == cont@longwin ) &
        all( nchar(cont@tailpats) + cont@leftwin <= cont@longwin ) &
        all( dim(cont@data) == dim(cont@projmatrix) ) &
        ( length(cont@mutrates) == length(cont@genmatrix@mutpats) ) &
        ( length(cont@selcoef) == length(cont@genmatrix@selpats) )
    )
}

# genmatrix extends the sparse matrix class, carrying along more information.
setClass("genmatrix", representation(
                                     muttrans="Matrix",
                                     seltrans="Matrix",
                                     bases="character",
                                     mutpats="list",
                                     selpats="list",
                                     boundary="character",
                                     fixfn="function"),
         contains = "dgCMatrix")

setClass("tuplecounts",representation(leftwin="numeric",counts="Matrix",bases="character"))
# things to make tuplecounts act like the matrix inside of it:
setMethod("dim", signature=c(x="tuplecounts"), definition=function (x) { dim(x@counts) } )
setMethod("dimnames", signature=c(x="tuplecounts"), definition=function (x) { dimnames(x@counts) } )
setMethod("dimnames<-", signature=c(x="tuplecounts",value="ANY"), definition=function (x,value) { dimnames(x@counts)<-value } )
setMethod("as.matrix", signature=c(x="tuplecounts"), definition=function (x) { as.matrix(x@counts) } )
setMethod("as.vector", signature=c(x="tuplecounts"), definition=function (x) { as.vector(x@counts) } )
setMethod("%*%", signature=c(x="tuplecounts",y="ANY"), definition=function (x,y) { x@counts %*% y } )
setMethod("%*%", signature=c(x="ANY",y="tuplecounts"), definition=function (x,y) { x %*% y@counts } )
setMethod("head", signature=c(x="tuplecounts"), definition=function (x) { head(x@counts) } )
setMethod("image", signature=c(x="tuplecounts"), definition=function (x) { image(x@counts) } )


setClass("context",
         representation(data="tuplecounts",
                        genmatrix="genmatrix",
                        mutrates="numeric",
                        selcoef="numeric",
                        params="numeric",
                        projmatrix="Matrix",
                        likfun="function",
                        results="list"
                    )
         )

setClass("contextMCMC", representation(
                                       mutprior="numeric",
                                       selprior="numeric"
                                       ),
         contains="context")

setMethod("dimnames", signature=c(x="context"), definition=function (x) { dimnames(x@data) } )

# We would like the likfun function to automatically have access to the stuff in the context object
#  ... where is 'self'?!?
# EVIL:
setGeneric("likfun", function (x) { standardGeneric("likfun") } )
setMethod("likfun", signature=c(x="context"), definition=function(x) {
          f <- x@likfun
          environment(f) <- list2env( list(genmatrix=x@genmatrix,projmatrix=x@projmatrix,counts=x@data), parent=globalenv())
          return(f) } )

# extract window lengths from these objects:
setGeneric("longwin", function(x) { standardGeneric("longwin") })
setGeneric("shortwin", function(x) { standardGeneric("shortwin") })
setGeneric("leftwin", function(x) { standardGeneric("leftwin") })
setMethod("longwin", signature=c(x="genmatrix"), definition=function(x) { nchar(rownames(x)[1]) } )
setMethod("longwin", signature=c(x="tuplecounts"), definition=function(x) { nchar(rownames(x@counts)[1]) } )
setMethod("longwin", signature=c(x="context"), definition=function(x) { longwin(x@data) } )
setMethod("shortwin", signature=c(x="tuplecounts"), definition=function(x) { nchar(colnames(x@counts)[1]) } )
setMethod("shortwin", signature=c(x="context"), definition=function(x) { shortwin(x@data) } )
setMethod("leftwin", signature=c(x="tuplecounts"), definition=function(x) { x@leftwin } )
setMethod("leftwin", signature=c(x="context"), definition=function(x) { leftwin(x@data) } )
# convenience functions
setGeneric("nmuts", function (x) { standardGeneric("nmuts") } )
setGeneric("nsel", function (x) { standardGeneric("nsel") } )
setMethod("nmuts", signature=c(x="genmatrix"), definition=function (x) { length(x@mutpats) } )
setMethod("nsel", signature=c(x="genmatrix"), definition=function (x) { length(x@selpats) } )
setMethod("nmuts", signature=c(x="context"), definition=function (x) { nmuts(x@genmatrix) } )
setMethod("nsel", signature=c(x="context"), definition=function (x) { nsel(x@genmatrix) } )

# and methods related to model fitting
setGeneric("tuplecounts", function (x) { standardGeneric("tuplecounts") } )
setMethod("tuplecounts", signature=c(x="context"), definition=function (x) {x@data} )
setMethod("coef", signature=c(object="context"), definition=function (object) {
          coef <- c( object@mutrates, object@selcoef, object@params )
          names(coef) <- c( mutnames( object@genmatrix@mutpats ), mutnames( object@genmatrix@selpats ), names(object@params) )
          return(coef) } )
setMethod("rowSums", signature=c(x="tuplecounts"), definition=function (x) { rowSums(x@counts) } )
setMethod("image", signature=c(x="tuplecounts"), definition=function (x) { image(x@counts) } )
setMethod("rowSums", signature=c(x="context"), definition=function (x) { rowSums(x@data@counts) } )
setMethod("fitted", signature=c(object="context"), definition=function (object,...) { predictcounts.context(object,...) } )
setMethod("residuals", signature=c(object="context"), definition=function (object,...) { resid.context(object,...) } )

## This makes everything go all to hell, for some reason:
# setMethod("-",signature=c("tuplecounts","ANY"), definition=function(e1,e2) { new("tuplecounts",leftwin=e1@leftwin,counts=Matrix(e1@counts-e2),bases=e1@bases) } )
# setMethod("-",signature=c("ANY","tuplecounts"), definition=function(e1,e2) { new("tuplecounts",leftwin=e2@leftwin,counts=Matrix(e1-e2@counts),bases=e2@bases) } )


####
#

makegenmatrix <- function (mutpats, selpats=list(), patlen=nchar(patterns[1]), patterns=getpatterns(patlen,bases), bases, mutrates=rep(1,length(mutpats)),selcoef=rep(1,length(selpats)), boundary="none", fixfn=function(...){1}, ...) {
    # YYY suggest mutpatll rather than mutpats
    # Make the generator matrix G on the specified set of patterns,
    # carrying with it the means to quickly update itself.
    #  DON'T do the diagonal, so that the updating is easier.
    if (!is.numeric(patlen)|(missing(patlen)&missing(patterns))) { stop("need patlen or patterns") }
    if ( (length(selpats)>0 && max(sapply(unlist(selpats),nchar))>patlen) | max(sapply(unlist(mutpats),nchar))>patlen ) { stop("some patterns longer than patlen") }
    # mutmats is a list of matrices, with one matrix for each of the mutpatls in mutpatll describing the induced mutation process on the specified patterns (see getmutmats function def'n)
    mutmats <- getmutmats(mutpats,patterns,boundary=boundary)
    allmutmats <- do.call( rbind, mutmats )
    # start converting to dgCMatrix format by getting the order of all of the "from" (i) and "to" (j) indices.
    dgCord <- order( allmutmats$j, allmutmats$i )
    # nmutswitches is a vector with kth component equal to the number of mutations between patterns that can be induced by performing one of the substitutions in the kth mutpat.
    nmutswitches <- sapply(mutmats,NROW)
    # muttrans: build a matrix to translate the specified mutation rates into rates on the specified pattern space.
    # We take these mutation patterns with repetition to be the implicit order on mutation patterns.
    # We put a 1 in every (i,j) such that i indexes a pattern mutation that happens at the jth mutpat.
    muttrans <- dgTtodgC( new( "dgTMatrix", i=1:sum(nmutswitches)-1L, j=rep(seq_along(mutrates),times=nmutswitches)-1L, x=rep(1,sum(nmutswitches)), Dim=c(sum(nmutswitches),length(mutrates)) ) )
    muttrans <- muttrans[ dgCord , , drop=FALSE ]
    # selection?
    if (length(selpats)>0) {
        selmatches <- getselmatches( selpats, patterns, boundary=boundary )
        # transfer selection coefficients to selective differences involved in each mutation
        #    these are ( transitions ) x ( mutpats ) matrix
        #     ... make these sparse?
        ##  the following has (fromsel-tosel); combined to reduce memory usage
        ##  fromsel <- selmatches[,  allmutmats$i, drop=FALSE ]
        ##  tosel <- selmatches[,  allmutmats$j, drop=FALSE ]
        seltrans <- Matrix(t( ( selmatches[,  allmutmats$i, drop=FALSE ] - selmatches[,  allmutmats$j, drop=FALSE ] )[,dgCord,drop=FALSE] ))
    } else {
        seltrans <- Matrix(numeric(0),nrow=nrow(muttrans),ncol=0)
    }
    # there may be duplicated rows (matching multiple patterns); deal with this
    # XXX EM doesn't get this code yet
    dups <- c(FALSE,diff(allmutmats[dgCord,][,1])==0)
    dupproj <- new("dgTMatrix",i=cumsum(!dups)-1L,j=seq_along(dups)-1L,x=rep(1,length(dups)),Dim=c(nrow(allmutmats)-sum(dups),nrow(allmutmats)))
    muttrans <- dupproj %*% muttrans
    seltrans <- dupproj %*% seltrans
    # full instantaneous mutation, and transition matrix
    genmatrix <- with( allmutmats[dgCord,], new( "genmatrix",
            i=(i-1L)[!dups],
            p=sapply(0:length(patterns), function (k) sum(j[!dups]<k+1)),
            x=rep(1,sum(!dups)),
            Dim=c(length(patterns),length(patterns)),
            muttrans=muttrans,
            seltrans=seltrans,
            mutpats=mutpats,
            selpats=selpats,
            bases=bases,
            fixfn=fixfn,
            boundary=boundary
        ) )
    rownames( genmatrix ) <- colnames( genmatrix ) <- patterns
    genmatrix@x <- update( genmatrix, mutrates, selcoef, ... )
    # diag(genmatrix) <- (-1)*rowSums(genmatrix)  # this makes genmatrix a dgCMatrix
    return(genmatrix)
}

update <- function (G, mutrates, selcoef, ...) {
    # use like: genmatrix@x <- update(genmatrix,...)
    fixprob <- if (length(selcoef)>0) { G@fixfn( as.vector(G@seltrans%*%selcoef), ... ) } else { 1 }
    as.vector( G@muttrans %*% mutrates ) * fixprob
}

collapsepatmatrix <- function (ipatterns, leftwin, shortwin=nchar(fpatterns[1]), rightwin=nchar(ipatterns[1])-shortwin-leftwin, fpatterns=getpatterns(nchar(ipatterns[1])-leftwin-rightwin,bases), bases ) {
    # YYY it seems like calling this `getprojmatrix` would make things more consistent with the way the other functions are named.
    # Construct the matrix U described in the tex.
    # ipatterns are the "input" patterns, while fpatterns are the "final" projected patterns
    # returns a (nbases)^k by (nbases)^{k-leftwin-rightwin} matrix projection matrix
    # mapping patterns onto the shorter patterns obtained by deleting leftwin characters at the start and rightwin characters at the end.
    # This function assumes that all input patterns are the same length.
    patlen <- nchar(ipatterns[1])
    shortwin <- patlen - leftwin - rightwin
    stopifnot(shortwin>0)
    matchpats <- match( substr(ipatterns,leftwin+1L,leftwin+shortwin), fpatterns )
    matchmatrix <- new( "dgTMatrix", i=(seq_along(ipatterns)-1L), j=(matchpats-1L), x=rep(1,length(ipatterns)), Dim=c(length(ipatterns),length(fpatterns)) )
    rownames(matchmatrix) <- ipatterns
    colnames(matchmatrix) <- fpatterns
    return( matchmatrix )
}

meangenmatrix <- function (leftwin,rightwin,patlen,...) {
    # create a generator matrix that averages over possible adjacent states
    longpatlen <- patlen+leftwin+rightwin
    genmat <- makegenmatrix(...,patlen=longpatlen)  # this is G
    if (longpatlen == patlen) {  # no need to do anything else...
        return(genmat)
    }
    projmat <- collapsepatmatrix(ipatterns=rownames(genmat),leftwin=leftwin,rightwin=rightwin, bases=genmatrix@bases)  # this is P
    # Divide entries by column sums, then transpose. Makes M, which is now has rows indexed by
    # the short version and columns the usual one.
    meanmat <- t( sweep( projmat, 2, colSums(projmat), "/" ) )
    pgenmat <- meanmat %*% genmat %*% projmat   # this is H = M G P
    ii <- pgenmat@i
    # XXX EM doesn't get this line... what is pgenmat@p?
    jj <- rep(1:ncol(pgenmat),times=diff(pgenmat@p)) - 1L
    nondiag <- ( ii != jj )
    # construct matrix to project from x values in big dgCMatrix to little one:
    #  note that H_ij = M_i. G  P_.j  = sum_kl M_ik G_kl P_lj
    #  ... and we have constructed G and H so we already know which elements are nonzero
    #    and can use this to find the linear transformation
    #    from nonzero elements of G (genmat) to nonzero elements of H (pgenmat)
    ij.H <- 1L + cbind( i=ii[nondiag], j=jj[nondiag] )
    ij.G <- 1L + cbind( i=genmat@i, j=rep(1:ncol(genmat),times=diff(genmat@p))-1L )
    pnonz <- Matrix( 0, nrow=nrow(ij.H), ncol=nrow(ij.G), sparse=TRUE )
    # for-loop to avoid memory hogging
    for (k in 1:nrow(ij.H)) { pnonz[k,] <-  meanmat[ij.H[k,"i"],ij.G[,"i"]] * projmat[ij.G[,"j"],ij.H[k,"j"]] }
    # pnonz <- t( apply( ij.H, 1, function (ij) { meanmat[ij[1],ij.G[,"i"]] * projmat[ij.G[,"j"],ij[2]] } ) )
    pp <- sapply( 0:ncol(pgenmat), function(k) sum(jj[nondiag]<k) )
    meangenmat <- new( "genmatrix", i=ii[nondiag], p=pp, x=pgenmat@x[nondiag], Dim=pgenmat@Dim, Dimnames=pgenmat@Dimnames,
        muttrans = (pnonz %*% genmat@muttrans), seltrans = (pnonz %*% genmat@seltrans),
        bases=pgenmat@bases, mutpats=pgenmat@mutpats, selpats=pgenmat@selpats,
        boundary=pgenmat@boundary, fixfn=pgenmat@fixfn
        )
    args <- list(...)
    if (is.null(args$mutrates)) { args$mutrates <- rep(1,length(genmat@mutpats)) }
    if (is.null(args$selcoef)) { args$selcoef <- rep(1,length(genmat@selpats)) }
    meangenmat@x <- do.call(update,c( list(G=meangenmat), args ) )
    return( meangenmat )
}

computetransmatrix <- function( genmatrix, projmatrix, tlen=1, shape=1, names=FALSE, transpose=FALSE, time="fixed", ... ) {
    # Compute the product of exp(tlen*genmatrix) and projmatrix, either on the left or the right (as transpose is true or false)
    #   either after a fixed time: exp(tlen*genmatrix)
    #   or after a gamma-distributed time
    if (time=="gamma") {
        if (is.null(dim(projmatrix))) { dim(projmatrix) <- c(length(projmatrix),1) }
        subtransmatrix <- gammaAtv( A=( if (transpose) { t(genmatrix) } else {genmatrix} ), scale=tlen, shape=shape, v=projmatrix )
    } else {
        totalrates <- rowSums(genmatrix)
        scale.t <- mean(totalrates)
        A <- (1/scale.t) * ( ( if (transpose) { t(genmatrix) } else {genmatrix} ) - Diagonal(nrow(genmatrix),totalrates) )
        if (is.null(dim(projmatrix))) { dim(projmatrix) <- c(length(projmatrix),1) }
        subtransmatrix <- sapply( 1:ncol(projmatrix), function (k) { expAtv( A=A, t=tlen*scale.t, v=projmatrix[,k] )$eAtv } )
    }
    if (names) {
        rownames(subtransmatrix) <- rownames(genmatrix)
        colnames(subtransmatrix) <- colnames(projmatrix)
    }
    return( subtransmatrix )
}

getupdowntrans <- function ( genmatrix, projmatrix, mutrates, selcoef, initfreqs, tlens=c(1,1), ... ) {
    # arguments are lists of two: first the "up" branch (leading from simpler summaries), second the "down"
    # returns matrix with entry [x,y] the probability of seeing y on the "down branch" given x was seen on the "up" branch.
    otherparams <- list(...)
    op.1 <- lapply( otherparams, "[[", 1 )
    op.2 <- lapply( otherparams, "[[", 2 )
    genmatrix.up <- genmatrix
    genmatrix.up@x <- do.call( update, c( list( G=genmatrix, mutrates=mutrates[[1]]*tlens[1], selcoef=selcoef[[1]] ), op.1 ) )
    genmatrix.down <- genmatrix
    genmatrix.down@x <- do.call( update, c( list( G=genmatrix, mutrates=mutrates[[2]]*tlens[2], selcoef=selcoef[[2]] ), op.2 ) )
    upbranch <- initfreqs * computetransmatrix( genmatrix.up, projmatrix )   #  prob of root, y
    downbranch <- computetransmatrix( genmatrix.down, initfreqs, transpose=TRUE )  # marginal prob of x
    return( computetransmatrix( genmatrix.down, upbranch, transpose=TRUE ) / as.vector(downbranch) )  # conditional prob of y given x
}

upbranch <- function ( genmatrix, tipmatrix, rootfreqs, mutrates, selcoef, tlen, ... ) {
    # "prune" a dangling edge in the "up" direction --
    #    return diag(rootfreqs) * e^( tlen * genmatrix(mut,sel) ) %*% tipmatrix .
    # So, result has
    #   upbranch[i,j] = rootfreqs[i] * E[ tipmatrix[X(tlen),j] | X(0) = i ]
    if (!missing(mutrates)) { genmatrix@x <- update( G=genmatrix, mutrates=mutrates*tlen, selcoef=selcoef, ... ) }
    upbranch <- rootfreqs * computetransmatrix( genmatrix, projmatrix )   #  prob of root, y
    return( upbranch )
}

downbranch <- function ( genmatrix, rootmatrix, mutrates, selcoef, tlen, ... ) {
    # Collapse down branches, rather than up -- like upbranch, but transposed --
    #   return e^( tlen * t(genmatrix(mut,sel)) ) %*% rootmatrix
    if (!missing(mutrates)) { genmatrix@x <- update( G=genmatrix, mutrates=mutrates*tlen, selcoef=selcoef, ... ) }
    downbranch <- computetransmatrix( genmatrix.down, rootmatrix, transpose=TRUE )  # marginal prob of x
    return(downbranch)
}

##
# pruning-ish algorithm
#
# outline:
#   1) begin at the root, with initial frequencies.
#   2) we'll move down, towards the focal (longest-pattern) tip
#   3) before moving, resolve up the the other, non-focal branch, recursive in the same way.

treetrans <- function (  ) {
    # First, compute a good order

}

###
# tree miscellany

get.descendants  <- function (tree) {
    # return the node-by-node matrix with [i,j] TRUE meaning that j is below i
    adjacency <- matrix( 0, nrow=Nnode(tree)+Ntip(tree), ncol=Nnode(tree)+Ntip(tree) )
    adjacency[ tree$edge ] <- 1  # adjacency matrix: [i,j] means that node j is directly below node i
    # descendants is indexed by NODES
    descendants <- apower <- adjacency  # [i,j] means that node j is somewhere below node i
    while ( any(apower>0) ) {
        apower <- apower %*% adjacency
        descendants <- descendants + apower
    }
    stopifnot( all( descendants %in% c(0,1) ) )
    descendants <- ( descendants > 0 )
    diag(descendants) <- TRUE
    return(descendants)
}


###
# genome-ey things

reverse.complement <- function (mutpats) {
    # return index for each mutpat of the reverse-complement mutation pattern (or NA if none)
    # throw error if there are multiple matches
    rmutpats <- lapply( mutpats, lapply, function (x) { chartr(sapply(strsplit(x,''),function(y){paste(rev(y),collapse='')}), old="ACGT", new="TGCA") } )
    mutpats <- lapply( mutpats, sapply, paste, collapse='->' )
    rmutpats <- lapply( rmutpats, sapply, paste, collapse='->' )
    rinds <- rep(NA,length(mutpats))
    for (k in seq_along(mutpats)) {
        matches <- c()
        for (j in seq_along(rmutpats)) {
            if ( all( mutpats[[k]] %in% rmutpats[[j]] ) ) { matches <- c(j,matches) }
        }
        stopifnot(length(matches) <= 1)
        if (length(matches)==1) { rinds[k] <- matches }
    }
    return(rinds)
}


###
#
predictcounts.context <- function (model, longwin=NULL, shortwin=NULL, leftwin=NULL, initcounts=rowSums(model), mutrates=model@mutrates, selcoef=model@selcoef, genmatrix=model@genmatrix, projmatrix=model@projmatrix ) {
    # default values not cooperating with S4 methods:
    if (is.null(longwin)) { longwin <- longwin(model) }
    if (is.null(shortwin)) { shortwin <- shortwin(model) }
    if (is.null(leftwin)) { leftwin <- leftwin(model) }
    if (!missing(genmatrix) && missing(projmatrix)) {
        projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, fpatterns=getpatterns(shortwin,genmatrix@bases) )
    }
    predictcounts(longwin,shortwin,leftwin,initcounts,mutrates,selcoef,genmatrix,projmatrix)
}


predictcounts <- function (longwin, shortwin, leftwin, initcounts, mutrates, selcoef, genmatrix, projmatrix, ... ) {
    # Compute expected counts of paired patterns:
    #  where the actual computation happens
    rightwin <- longwin-shortwin-leftwin
    if (!missing(mutrates)||!missing(selcoef)) { genmatrix@x <- update(genmatrix,mutrates=mutrates,selcoef=selcoef,...) }
    if (missing(projmatrix)) { projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin, bases=genmatrix@bases ) }
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE, ... )
    fullcounts <- initcounts * subtransmatrix
    return( new("tuplecounts", leftwin=leftwin, counts=Matrix(fullcounts), bases=genmatrix@bases ) )
}

projectcounts <- function( counts, lcountwin, countwin, rcountwin ) {
    # using counts for c(leftwin,shortwin,rightwin) compute counts for c(lcountwin,countwin,rcountwin).
    #   valid ranges for parameters are
    #    (l-lc)^+ <= k < (l+w)-(lc+wc)+(r-rc)^-
    leftwin <- leftwin(counts)
    longwin <- longwin(counts)
    shortwin <- shortwin(counts)
    rightwin <- longwin-shortwin-leftwin
    if ( max(0L,leftwin-lcountwin) > (leftwin+shortwin)-(lcountwin+countwin)+min(0L,rightwin-rcountwin) ) {
        stop("unreconcilable windows specified.")
    }
    pcounts <- matrix(0,nrow=npatterns(lcountwin+countwin+rcountwin,counts@bases),ncol=npatterns(countwin,counts@bases))
    for (k in max(0L,leftwin-lcountwin):((leftwin+shortwin)-(lcountwin+countwin)+min(0L,rightwin-rcountwin))) {
        lpmat <- collapsepatmatrix( ipatterns=rownames(counts), leftwin=k, rightwin=longwin-(k+lcountwin+countwin+rcountwin), bases=counts@bases )
        rpmat <- collapsepatmatrix( ipatterns=colnames(counts), leftwin=k+lcountwin-leftwin, rightwin=shortwin-(k+lcountwin-leftwin+countwin), bases=counts@bases )
        pcounts <- pcounts + t(lpmat) %*% counts %*% (rpmat)
    }
    dimnames(pcounts) <- list( colnames(lpmat), colnames(rpmat) )
    return( new("tuplecounts", leftwin=lcountwin, counts=Matrix(pcounts), bases=counts@bases) )
}

predicttreecounts <- function (shortwin, leftwin=0, rightwin=0, initcounts, mutrates, selcoef, mutpats, selpats, tlens, genmatrix, projmatrix, initfreqs, patcomp, ... ) {
    # Compute expected counts of paired patterns:
    longwin <- leftwin+shortwin+rightwin
    if (missing(genmatrix)) { genmatrix <- makegenmatrix(patlen=leftwin+shortwin+rightwin,mutpats=mutpats,selpats=selpats, ...) }
    if (missing(projmatrix)) { projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin ) }
    if (missing(patcomp) & !is.null(names(initfreqs))) {
        patcomp <- apply( do.call(rbind, strsplit(rownames(genmatrix),'') ), 2, match, names(initfreqs) )  # which base is at each position in each pattern
    }
    patfreqs <- initfreqs[patcomp]
    dim(patfreqs) <- dim(patcomp)
    patfreqs <- apply( patfreqs, 1, prod )
    updownbranch <- getupdowntrans( genmatrix, projmatrix, mutrates=mutrates, selcoef=selcoef, initfreqs=patfreqs, tlens=tlens )
    if (missing(initcounts)) { initcounts <- 1 }
    fullcounts <- initcounts * updownbranch
    dimnames(fullcounts) <- dimnames( projmatrix )
    return( fullcounts )
}

###
# stuff for looking at residuals and finding motifs there

resid.context <- function (object, counts=object@data, genmatrix=object@genmatrix, pretty=FALSE) {
    expected <- fitted(object, initcounts=rowSums(counts),
                       longwin=longwin(counts), shortwin=shortwin(counts), leftwin=leftwin(counts),
                       genmatrix=genmatrix )
    stopifnot( all( rownames(expected)==rownames(counts) ) && all( colnames(expected)==colnames(counts) ) )
    resids <- counts
    resids@counts <- ( counts@counts - expected@counts )
    return( resids )
}

listresids <- function (counts, expected, file, trim=20, leftwin=(nchar(rownames(counts)[1])-nchar(colnames(counts)[1]))/2) {
    # make a readable output ordered by z-score
    #  optionally writing results out to 'file'
    #  and trimming to only patterns with z-score above 'trim'
    resids <- data.frame( inpat=rownames(counts)[row(counts)],
                         outpat=colnames(counts)[col(counts)],
                         observed=as.vector(counts),
                         expected=as.vector(expected),
                    resid=as.vector(counts-expected) )
    longwin <- nchar(rownames(counts)[1])
    resids$inpat <- with(resids, paste( tolower(substr(inpat,1,leftwin)), substr(inpat,leftwin+1,longwin-leftwin), tolower(substr(inpat,longwin-leftwin+1,longwin)), sep='' ) )
    resids$outpat <- paste(resids$outpat)
    resids$z <- resids$resid/sqrt(as.vector(expected))
    if (is.numeric(trim)) { resids <- subset( resids, is.numeric(resids$z) & (abs(resids$z) > trim) ) }
    resids <- resids[order(resids$z),]
    if (!missing(file)) {
        write.table(file=pipe(paste("column -t >", file)), x=resids, sep=' ',quote=FALSE, row.names=FALSE )
        return(invisible(resids))
    } else {
        return(resids)
    }
}

clusterresids <- function (resids,npats=300,nclusts=12) {
    # look for motifs in resids (as above)
    invisible( lapply( c(+1,-1), function (sign) {
            resids$z <- resids$z * sign
            trim <- quantile((resids$z[is.finite(resids$z)]),1-npats/sum(is.finite(resids$z)))
            resids <- subset(resids, is.finite(z) & z>trim )
            longwin <- nchar(resids$inpat[1])
            shortwin <- nchar(resids$outpat[1])
            pats <- paste(resids$inpat,resids$outpat,sep="")
            # 10 secs for 3000x3000
            sdists <- stringdistmatrix(pats,pats,method="hamming")
            rownames(sdists) <- colnames(sdists) <- paste(resids$inpat,resids$outpat,sep="|")
            sclust <- cutree( hclust(as.dist(sdists)), k=nclusts )
            motifs <- lapply( 1:nclusts, function (k) print.motif(names(sclust)[sclust==k],print=FALSE,weights=sqrt(resids$z[sclust==k])) )
            cat( c("", paste(do.call( paste, c( motifs, list(sep="   ") ) ),"\n") ) )
            return(motifs)
        } ) )
}

table.weighted <- function(x,weights=rep(1,length(x))) {
    # like table, but with weights
    x <- factor(x)
    tx <- table(x)
    return( tx * tapply(weights,x,sum,na.rm=TRUE) )
}

print.motif <- function (pats,weights=1,n=24,print=TRUE,long=FALSE) {
    pats <- do.call( rbind, strsplit(pats,"") )
    freqs <- lapply( 1:ncol(pats), function(k) {
                    x <- table.weighted(pats[,k],weights); x[order(names(x))]
                } )
    if (long) {
        samps <- do.call( cbind, lapply( freqs, function (x) { names(x)[cut((1:n)/(n+1),breaks=c(0,cumsum(x))/sum(x))] } ) )
        samps <- apply(samps,1,paste,collapse='')
    } else {
        samps <- paste( sapply(freqs, function (x) {
                               themax <- which.max(x)
                               return(
                                   if (max(x)>.75*sum(x)) {
                                       (toupper(names(x)[themax]))
                                   } else if (max(x)>.5*sum(x)) {
                                       (toupper(names(x)[themax]))
                                   } else { "." }
                               )
        } ), collapse='' )
    }
    if (print) { for (x in samps) { cat(x,"\n") } }
    return(invisible(samps))
}

## other stuff

whichchanged <- function (ipatterns,fpatterns,leftwin=0,shortwin=nchar(ipatterns[0])) {
    # return indicator corresponding to entries of output of gettransmatrix that have changed
    if (!is.null(dimnames(ipatterns))) { fpatterns <- colnames(ipatterns); ipatterns <- rownames(ipatterns) }
    return( outer( ipatterns, fpatterns, function (x,y) { substr(x,leftwin+1,leftwin+shortwin)!=y } ) )
}

leftchanged <- function (ipatterns,fpatterns,leftwin=0,shortwin=nchar(ipatterns[0])) {
    # return indicator corresponding to whether average of changed positions is left of middle
    # ... any two such patterns, if they overlap, require an additional change.
    if (!is.null(dimnames(ipatterns))) { fpatterns <- colnames(ipatterns); ipatterns <- rownames(ipatterns) }
    changedchars <- sapply( 1:shortwin, function (k) outer( ipatterns, fpatterns, function (x,y) { ifelse( substr(x,leftwin+k,leftwin+k)!=substr(y,k,k), k, NA ) } ) )
    dim(changedchars) <- c( length(ipatterns), length(fpatterns), shortwin )
    meanpos <- rowMeans(changedchars, na.rm=TRUE, dims=2)
    return( !is.na(meanpos) & meanpos <= (shortwin+1)/2 )
}


# Misc

wrapsubstr <- function (x,start,stop) {
    if (all(nchar(x)==0)) { return("") }
    while( any(nchar(x)<stop) ) {
        x <- paste(x,x,sep='')
    }
    substr(x,start,stop)
}

"wrapsubstr<-" <- function (x,start,stop,value) {
    if (length(value)<length(x)) { value <- rep(value,length(x)) }
    xlen <- nchar(x)
    k <- rep(1,length(x))
    while( TRUE ) {
        stop <- ifelse( xlen>0, stop - xlen*((start-1)%/%xlen), 0 )
        start <- ifelse( xlen>0, (start-1)%%xlen+1, 0 )
        thisstop <- pmin(xlen,stop)
        substr(x,start,thisstop) <- substr(value,k,k+thisstop-start)
        k <- k+thisstop-start+1
        start <- thisstop+1
        if( all(stop<start) ) { break; }
    }
    return(x)
}

dgTtodgC <- function (M) {
    # convert between matrix classes.  For understanding.
    ijx <- data.frame( i=M@i, j=M@j, x=M@x )
    ijx <- ijx[ order( ijx$j, ijx$i ), ]
    with(ijx, new( "dgCMatrix", i=i, p=sapply(0:ncol(M), function(k) sum(j<k)), x=x, Dim=dim(M) ) )
}


