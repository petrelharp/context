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
source(paste(.PATH,"/gammaAtv.R",sep=''))  # fixed upstream
source(paste(.PATH,"/input-output.R",sep='')) 

getpatterns <- function(patlen,bases) {
    # construct a list of all patterns
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

getmutpats <- function(patlen,nchanges=1) {
    # all kmer -> kmer changes
    # that involve at most nchanges changes
    mutpats <- list()
    patterns <- getpatterns(patlen,bases)
    for (k in 1:patlen) {
        kmers <- getpatterns(k,bases)
        mutpats <- c( mutpats,
                apply(combn(kmers,2),2,list), # make lists of rows (2 = apply over columns), giving 2-element lists of kmers
                apply(combn(kmers,2)[2:1,],2,list) # and the reverse of the 2-element lists
            )
    }
    obschanges <- sapply(mutpatchanges(mutpats),nrow)
    return( mutpats[obschanges%in%nchanges] )
}

npatterns <- function (patlen,bases) {
    return( length(bases)^patlen )
}

mutnames <- function (mutpats) {
    return( unlist( sapply( sapply( mutpats, lapply, paste, collapse="->" ), paste, collapse="|" ) ) )
}

###
# Point estimates

divergence <- function (counts, lwin) {
    # mean density of nucleotide differences
    require(stringdist) # http://www.inside-r.org/packages/cran/stringdist/docs/stringdist
    patlen <- nchar(colnames(counts)[1])
    nchanges <- stringdistmatrix( substr(rownames(counts),lwin+1,lwin+patlen), colnames(counts), method="hamming" )
    return( sum(nchanges*counts)/(sum(counts)*patlen) )
}

countmuts <- function (counts, mutpats, lwin, ...) {
    # given a contingency table of kmer changes from data, a collection of
    # mutation patterns, the lwin, and a list of arguments to `sum`, this function
    # returns a matrix, the first row of which gives counts of how many of these could be produced by each of mutpats
    #   and the second of which gives the number of "from" matches of the mutpats (called the "total possible").
    #
    # in other words, for each mutation pattern a -> b
    #  sum the values of counts[u,v] over choices of u,v such that:
    #   (i) 'a' matches 'u' at some position
    #   (ii) in addition, 'b' matches 'v' at the same position
    #
    # note that if we estimate rates by
    #    r.est <- countmuts(...)[1,]/countmuts(...)[2,]
    # then something like
    #    sum( r.est ) / 4
    # should be close to the mean density of nucleotide changes
    #    divergence(...)
    counts <- as.matrix(counts)
    # `win` is the length of the inner window
    win <- nchar(colnames(counts)[1])
    stopifnot( length(win)>0 & win>0 ) # length statement catches the case that there are no colnames for counts
    # trim off windows
    xx <- substr(rownames(counts),lwin+1,lwin+win)
    yy <- colnames(counts)
    # Note that `observed` counts a change multiply if it can occur in different ways.
    observed <- possible <- numeric(length(mutpats)) # observed and possible are each numeric vectors
    for (j in seq_along(mutpats)) {
        mutpat <- mutpats[[j]]
        for (k in 1:(win-1)) {
            mx <- sapply(mutpat, function (mp) {
                    patlen <- nchar(mp[1])
                    if ( k+patlen-1 <= win ) {
                        sum( counts[ ( substr(xx,k,k+patlen-1) == mp[1] ), ( substr(yy,k,k+patlen-1) == mp[2] ) ], ... )
                    } else {
                        0
                    }
                } )
            observed[j] <- observed[j] + sum( mx, ... )
            px <- sapply(mutpat, function (mp) {
                    patlen <- nchar(mp[1])
                    if ( k+patlen-1 <= win ) {
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

getmutmats <- function(mutpats,patterns,boundary=c("none","wrap")) {
    # returns i and j's indexing the generator matrix.
    # that is, given a list of mutation patterns,
    #   which can be either pairs or lists of pairs,
    # return a corresponding list of matrices with (1-based) indices of changes corresponding to mutation patterns
    #   i.e. if (i,j) is a row of output[[k]], then patterns[j] can be obtained from patterns[i]
    #   by performing the substitution from mutpats[[k]][1] -> mutpats[[k]][2]
    #   at some location within the string.
    boundary <- match.arg(boundary)
    winlen <- nchar(patterns[1])
    mutpats <- lapply( mutpats, function (x) { if (is.list(x)) { x } else { list(x) } } )
    lapply( mutpats, function (y) {  # y is list of short pattern pairs
            do.call( rbind, lapply(y, function (x) {  # x is short pattern pair (from, to)
                patlen <- nchar(x[1])
                switch( boundary,
                    wrap={ # patterns are circular
                        wpatterns <- paste( patterns, substr(patterns,1,patlen), sep='' )
                        maxshift <- winlen
                        subsfun <- function (pat,topat,k) { wrapsubstr( pat, k, k+patlen-1 ) <- topat; return(pat) }
                    },
                    none={
                        wpatterns <- patterns
                        maxshift <- winlen-patlen+1
                        subsfun <- function (pat,topat,k) { substr( pat, k, k+patlen-1 ) <- topat; return(pat) }
                    }
                )
                do.call( rbind, lapply( 1:maxshift, function (k) {  # k is position of short pattern in long pattern
                        i <- which( substr( wpatterns, k, k+patlen-1 ) == x[1] ) # which patterns match short from-pattern?
                        replstr <- subsfun( patterns[i], x[2], k ) # substitute in to-pattern
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
# genmatrix extends the sparse matrix class, carrying along more information.

check.context <- function (cont) {
    return(
        all( rownames(cont@genmatrix) == cont@headpats ) &
        all( rownames(cont@projmatrix) == cont@headpats ) &
        all( colnames(cont@projmatrix) == cont@tailpats ) &
        all( nchar(cont@headpats) == cont@winlen ) &
        all( nchar(cont@tailpats) + cont@lwin <= cont@winlen ) &
        all( dim(cont@data) == dim(cont@projmatrix) ) &
        ( length(cont@mutrates) == length(cont@genmatrix@mutpats) ) &
        ( length(cont@selcoef) == length(cont@genmatrix@selpats) )
    )
}

# genmatrix extends the sparse matrix class, carrying along more information.
setClass("genmatrix", representation(
                                     muttrans="Matrix",
                                     seltrans="Matrix",
                                     mutpats="list",
                                     selpats="list",
                                     boundary="character",
                                     fixfn="function"),
         contains = "dgCMatrix")

setClass("tuplecounts",representation(lwin="numeric",counts="Matrix"))
setMethod("dim", signature=(x="tuplecounts"), definition=function (x) { dim(x@counts) } )
setMethod("dimnames", signature=(x="tuplecounts"), definition=function (x) { dimnames(x@counts) } )
setMethod("dimnames<-", signature=c(x="tuplecounts",value="ANY"), definition=function (x,value) { dimnames(x@counts)<-value } )
setMethod("as.matrix", signature=(x="tuplecounts"), definition=function (x) { as.matrix(x@counts) } )
setMethod("as.vector", signature=(x="tuplecounts"), definition=function (x) { as.vector(x@counts) } )

setClass("context",
         representation(data="tuplecounts",
                        genmatrix="genmatrix",
                        mutrates="numeric",
                        selcoef="numeric",
                        params="numeric",
                        projmatrix="Matrix",
                        likfun="function",
                        optim.results="list"
                    )
         )
setMethod("dimnames", signature=c(x="context"), definition=function (x) { dimnames(x@data) } )

# We would like the likfun function to automatically have access to the stuff in the context object
#  ... where is 'self'?!?
# EVIL:
setGeneric("likfun", function (x) { standardGeneric("likfun") } )
setMethod("likfun", signature=c(x="context"), definition=function(x) { 
          f <- x@likfun
          environment(f) <- list2env( list(genmatrix=x@genmatrix,projmatrix=x@projmatrix), parent=globalenv())
          return(f) } )

# extract window lengths from these objects:
setGeneric("winlen", function(x) { standardGeneric("winlen") })
setGeneric("win", function(x) { standardGeneric("win") })
setGeneric("lwin", function(x) { standardGeneric("lwin") })
setMethod("winlen", signature=c(x="genmatrix"), definition=function(x) { nchar(rownames(genmatrix)[1]) } )
setMethod("winlen", signature=c(x="tuplecounts"), definition=function(x) { nchar(rownames(x@counts)[1]) } )
setMethod("winlen", signature=c(x="context"), definition=function(x) { winlen(context@data) } )
setMethod("win", signature=c(x="tuplecounts"), definition=function(x) { nchar(colnames(x@counts)[1]) } )
setMethod("win", signature=c(x="context"), definition=function(x) { win(context@data) } )
setMethod("lwin", signature=c(x="tuplecounts"), definition=function(x) { x@lwin } )
setMethod("lwin", signature=c(x="context"), definition=function(x) { lwin(context@data) } )
# convenience functions
setGeneric("nmuts", function (x) { standardGeneric("nmuts") } )
setGeneric("nsel", function (x) { standardGeneric("nsel") } )
setMethod("nmuts", signature=c(x="genmatrix"), definition=function (x) { length(x@mutpats) } )
setMethod("nsel", signature=c(x="genmatrix"), definition=function (x) { length(x@selpats) } )

makegenmatrix <- function (mutpats, selpats=list(), patlen=nchar(patterns[1]), patterns=getpatterns(patlen,bases), bases, mutrates=rep(1,length(mutpats)),selcoef=rep(1,length(selpats)), boundary="none", fixfn=function(...){1}, ...) {
    # make the generator matrix, carrying with it the means to quickly update itself.
    #  DON'T do the diagonal, so that the updating is easier.
    if (!is.numeric(patlen)|(missing(patlen)&missing(patterns))) { stop("need patlen or patterns") }
    if ( (length(selpats)>0 && max(sapply(unlist(selpats),nchar))>patlen) | max(sapply(unlist(mutpats),nchar))>patlen ) { stop("some patterns longer than patlen") }
    # list of matrices with indices of changes corresponding to mutation patterns
    mutmats <- getmutmats(mutpats,patterns,boundary=boundary)
    allmutmats <- do.call( rbind, mutmats )
    # convert to dgCMatrix format
    dgCord <- order( allmutmats$j, allmutmats$i )
    # use this to transfer these to list of values in mutation matrix
    nmutswitches <- sapply(mutmats,NROW)
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

collapsepatmatrix <- function (ipatterns, lwin, win=nchar(fpatterns[1]), rwin=nchar(ipatterns[1])-win-lwin, fpatterns=getpatterns(nchar(ipatterns[1])-lwin-rwin), bases ) {
    # returns a (nbases)^k x (nbases)^k-m matrix projection matrix
    # map patterns onto the shorter patterns obtained by deleting lwin characters at the start and rwin characters at the end
    patlen <- nchar(ipatterns[1])
    win <- patlen - lwin - rwin
    stopifnot(win>0)
    matchpats <- match( substr(ipatterns,lwin+1L,lwin+win), fpatterns )
    matchmatrix <- new( "dgTMatrix", i=(seq_along(ipatterns)-1L), j=(matchpats-1L), x=rep(1,length(ipatterns)), Dim=c(length(ipatterns),length(fpatterns)) )
    rownames(matchmatrix) <- ipatterns
    colnames(matchmatrix) <- fpatterns
    return( matchmatrix )
}

meangenmatrix <- function (lwin,rwin,patlen,...) {
    # create a generator matrix that averages over possible adjacent states
    longpatlen <- patlen+lwin+rwin
    genmat <- makegenmatrix(...,patlen=longpatlen)  # this is G
    if (longpatlen == patlen) {  # no need to do anything else...
        return(genmat)
    }
    projmat <- collapsepatmatrix(ipatterns=rownames(genmat),lwin=lwin,rwin=rwin)  # this is P
    meanmat <- t( sweep( projmat, 2, colSums(projmat), "/" ) )  # this is M
    pgenmat <- meanmat %*% genmat %*% projmat   # this is H = M G P
    ii <- pgenmat@i
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
        muttrans = (pnonz %*% genmat@muttrans), seltrans = (pnonz %*% genmat@seltrans) )
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

predictcounts <- function (win, lwin=0, rwin=0, initcounts, mutrates, selcoef, mutpats, selpats, genmatrix, projmatrix, ... ) {
    # Compute expected counts of paired patterns:
    winlen <- lwin+win+rwin
    if (missing(genmatrix)) { genmatrix <- makegenmatrix(patlen=lwin+win+rwin,mutpats=mutpats,selpats=selpats, ...) }
    if (!missing(mutrates)|!missing(selcoef)) { genmatrix@x <- update(genmatrix,mutrates=mutrates,selcoef=selcoef,...) }
    if (missing(projmatrix)) { projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin ) }
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE, ... )
    if (missing(initcounts)) { initcounts <- 1 }
    fullcounts <- initcounts * subtransmatrix
    return( fullcounts )
}

projectcounts <- function( lwin, lcountwin, countwin, rcountwin, counts ) {
    # collapse: valid shift ranges is
    #  (l-lc)^+ <= k < (l+w)-(lc+wc)+(r-rc)^-
    #  using counts for (lwin,win,rwin) compute counts for (lcountwin,countwin,rcountwin).
    if ( max(0L,lwin-lcountwin) > (lwin+win)-(lcountwin+countwin)+min(0L,rwin-rcountwin) ) { stop("unreconcilable windows specified.") }
    winlen <- nchar(rownames(counts)[1])
    win <- nchar(colnames(counts)[1])
    rwin <- winlen-win-lwin
    pcounts <- matrix(0,nrow=npatterns(lcountwin+countwin+rcountwin,bases),ncol=npatterns(countwin,bases))
    for (k in max(0L,lwin-lcountwin):((lwin+win)-(lcountwin+countwin)+min(0L,rwin-rcountwin))) {
        lpmat <- collapsepatmatrix( ipatterns=rownames(counts), lwin=k, rwin=winlen-(k+lcountwin+countwin+rcountwin) )
        rpmat <- collapsepatmatrix( ipatterns=colnames(counts), lwin=k+lcountwin-lwin, rwin=win-(k+lcountwin-lwin+countwin) )
        pcounts <- pcounts + t(lpmat) %*% counts %*% (rpmat)
    }
    dimnames(pcounts) <- list( colnames(lpmat), colnames(rpmat) )
    return(pcounts)
}

predicttreecounts <- function (win, lwin=0, rwin=0, initcounts, mutrates, selcoef, mutpats, selpats, tlens, genmatrix, projmatrix, initfreqs, patcomp, ... ) {
    # Compute expected counts of paired patterns:
    winlen <- lwin+win+rwin
    if (missing(genmatrix)) { genmatrix <- makegenmatrix(patlen=lwin+win+rwin,mutpats=mutpats,selpats=selpats, ...) }
    if (missing(projmatrix)) { projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin ) }
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

listresids <- function (counts, expected, file, trim=20, lwin=(nchar(rownames(counts)[1])-nchar(colnames(counts)[1]))/2) {
    # make a readable output ordered by z-score
    #  optionally writing results out to 'file'
    #  and trimming to only patterns with z-score above 'trim'
    resids <- data.frame( inpat=rownames(counts)[row(counts)],
                         outpat=colnames(counts)[col(counts)],
                         observed=as.vector(counts),
                         expected=as.vector(expected),
                    resid=as.vector(counts-expected) )
    winlen <- nchar(rownames(counts)[1])
    resids$inpat <- with(resids, paste( tolower(substr(inpat,1,lwin)), substr(inpat,lwin+1,winlen-lwin), tolower(substr(inpat,winlen-lwin+1,winlen)), sep='' ) )
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
            winlen <- nchar(resids$inpat[1])
            win <- nchar(resids$outpat[1])
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

whichchanged <- function (ipatterns,fpatterns,lwin=0,win=nchar(ipatterns[0])) {
    # return indicator corresponding to entries of output of gettransmatrix that have changed
    if (!is.null(dimnames(ipatterns))) { fpatterns <- colnames(ipatterns); ipatterns <- rownames(ipatterns) }
    return( outer( ipatterns, fpatterns, function (x,y) { substr(x,lwin+1,lwin+win)!=y } ) )
}

leftchanged <- function (ipatterns,fpatterns,lwin=0,win=nchar(ipatterns[0])) {
    # return indicator corresponding to whether average of changed positions is left of middle
    # ... any two such patterns, if they overlap, require an additional change.
    if (!is.null(dimnames(ipatterns))) { fpatterns <- colnames(ipatterns); ipatterns <- rownames(ipatterns) }
    changedchars <- sapply( 1:win, function (k) outer( ipatterns, fpatterns, function (x,y) { ifelse( substr(x,lwin+k,lwin+k)!=substr(y,k,k), k, NA ) } ) )
    dim(changedchars) <- c( length(ipatterns), length(fpatterns), win )
    meanpos <- rowMeans(changedchars, na.rm=TRUE, dims=2)
    return( !is.na(meanpos) & meanpos <= (win+1)/2 )
}

getlikfun <- function (nmuts,nsel,genmatrix,projmatrix,const=0) {
    return( function (params) {
        # params are: mutrates, selcoef, Ne
        mutrates <- params[1:nmuts]
        selcoef <- params[nmuts+(1:nsel)]
        Ne <- params[nmuts+nsel+1]
        # tlen <- params[nmuts+nsel+2]  # confounded.
        # this is collapsed transition matrix
        genmatrix@x <- update(genmatrix,mutrates,selcoef,Ne)
        subtransmatrix <- computetransmatrix( genmatrix, projmatrix )
        # return negative log-likelihood
        (-1) * sum( counts * log(subtransmatrix) ) + const
    } )
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


# Unused?
if (FALSE ){

regexplen <- function (xx) {
    # length of the string matching a regexp that uses only "." and "[...]" (no other special characters!)
    sapply( xx, function (x) {
        y <- diff( c(0,grep("[]\\[]",strsplit(x,"")[[1]],value=FALSE),nchar(x)+1) ) - 1  # lengths of bits in and out of "[]"s
        sum( y[(1 ==  (1:length(y))%%2)] ) + (length(y)-1)/2
    } )
}

gettransmatrix <- function (mutpats, mutrates, selpats, selcoef, Ne, tlen=1, win, lwin=0, rwin=0, expm=expm.poisson, ... ) {
    # get reduced transition matrix: given (lwin, win, rwin) context, return probability of pattern in win
    #   note: alternative is expm=expm::expm(x,method="Higham08")
    winlen <- lwin+win+rwin
    fullgenmatrix <- makegenmatrix( mutpats, selpats, patlen=winlen,...)
    fullgenmatrix@x <- update(fullgenmatrix,mutrates,selcoef,Ne)
    projmatrix <- collapsepatmatrix( ipatterns=rownames(fullgenmatrix), lwin=lwin, rwin=rwin )
    subtransmatrix <- computetransmatrix( genmatrix, tlen, projmatrix, names=TRUE )
    # transmatrix <- expm( tlen * (fullgenmatrix-Diagonal(nrow(fullgenmatrix),rowSums(fullgenmatrix))) )  # exponentiate
    # subtransmatrix <- transmatrix %*% projmatrix        # collapse
    return( subtransmatrix )
}

}
