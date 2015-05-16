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
source(paste(.PATH,"/helper-fns.R",sep=''))
source(paste(.PATH,"/plotting-fns.R",sep=''))

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
    return( unlist( sapply( lapply( mutpats, lapply, paste, collapse="->" ), paste, collapse="|" ) ) )
}

selnames <- function (selpats) {
    # Stringify a list of selpats
    return( sapply(selpats,paste,collapse="|") )
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

countmuts <- function (counts, mutpats, leftwin=leftwin(counts), ...) {
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
    # return a corresponding list of two-column matrices with (1-based) indices of changes corresponding to mutation patterns
    #   i.e. if (i,j) is a row of output[[k]], then patterns[j] can be obtained from patterns[i]
    #   by performing a substitution from mutpats[[k]] at some location within the string.
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

getselmatches <- function (selpats, patterns, selfactors, boundary=c("none","wrap"), names=FALSE) {
    # selpats can be a vector or a list of vectors,
    #   each element gets one parameter
    # selmatches[i,j] is the sum of selfactors[[i]] multiplied by the number of times the corresponding thing in selpat[[i]] matches pattern[j]
    #  ... so if selfactors is all 1's, then :
    #     selmatches[i,j] is number of times anything in selpat[[i]] matches pattern[j]  
    boundary <- match.arg(boundary)
    substrfun <- switch( boundary, wrap=wrapsubstr, none=substr )
    patlen <- nchar(patterns[1])
    if (!is.list(selpats)) { selpats <- as.list(selpats) }
    selmatches <- t( do.call( cbind.dsparseVector, lapply(seq_along(selpats), function (i) {   # loop over selpats
            y <- selpats[[i]]
            Matrix::rowSums( do.call( cbind.dsparseVector, lapply(seq_along(y), function (j) {  # loop over elements of selpats[[y]]
                    x <- y[j]
                    maxshift <- patlen - switch( boundary, wrap=1, none=nchar(x) )
                    nonz.list <- lapply( 0:maxshift, function (k) { which(x == substrfun( patterns, 1+k, k+nchar(x) ) ) } ) 
                    return( selfactors[[i]][j] * Matrix::rowSums( sparseMatrix( i=unlist(nonz.list), p=c(0,cumsum(sapply(nonz.list,length))), dims=c(length(patterns),1+maxshift) ) , sparseResult=TRUE ) )
                } ) ), sparseResult=TRUE )
        } ) ) )
    # selmatches <- do.call( rbind, lapply(selpats, function (y) {
    #         rowSums( sapply(y, function (x) {
    #                 maxshift <- patlen - switch( boundary, wrap=1, none=nchar(x) )
    #                 rowSums( sapply( 0:maxshift, function (k) {
    #                             x == substrfun( patterns, 1+k, k+nchar(x) )
    #                             # xx <- paste( c(rep(".",k), x, rep(".", patlen-regexplen(x)-k)), collapse='' )
    #                             # grepl( xx, patterns )
    #                     } ) )
    #             } ) )
    #     } ) )
    if (names) {
        rownames(selmatches) <- names(selpats)
        colnames(selmatches) <- patterns
    }
    return(selmatches)
}


shape.fixfn <- function (ds,Ne,a,...) {
    # total influx of fixation given selection coefficient difference ds = a * |s[to] - s[from]|
    if (length(ds)==0) { 1 } else { ifelse( ds==0, 1, Ne*expm1(-2*a*abs(ds))/expm1(-2*Ne*a*abs(ds)) ) }
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
        all( dim(cont@counts) == dim(cont@projmatrix) ) &
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

# tuplecounts is a (2D, flattened) contingency table
setClass("tuplecounts",representation(
            leftwin="numeric",
            counts="Matrix",
            bases="character",
            rowtaxon="character",
            colpatterns="data.frame"),
        prototype=list(rowtaxon="long")
    )

# extractor functions
setGeneric("rowtaxon", function(x) { standardGeneric("rowtaxon") })
setMethod("rowtaxon", signature=c(x="tuplecounts"), definition=function (x) { x@rowtaxon } )
setGeneric("coltaxa", function(x) { standardGeneric("coltaxa") })
setMethod("coltaxa", signature=c(x="tuplecounts"), definition=function (x) { colnames(x@colpatterns) } )
setGeneric("colpatterns", function(x) { standardGeneric("colpatterns") })
setMethod("colpatterns", signature=c(x="tuplecounts"), definition=function (x) { x@colpatterns } )
setGeneric("counts", function(x) { standardGeneric("counts") })
setMethod("counts", signature=c(x="tuplecounts"), definition=function (x) { x@counts } )
# this method repackages the counts in a more friendly-looking data frame
setGeneric("countframe", function(x) { standardGeneric("countframe") })
setMethod("countframe", signature=c(x="tuplecounts"), definition=function (x) {
        cf <- cbind(
                data.frame( rep.int(rownames(x@counts),ncol(x@counts)) ),
                colpatterns(x)[ rep(1:nrow(colpatterns(x)),each=nrow(x@counts)), ],
                as.numeric(x@counts)
            )
        colnames(cf) <- c( rowtaxon(x), coltaxa(x), "count" )
        return(cf)
    } )

# things to make tuplecounts act like the matrix inside of it:
setMethod("dim", signature=c(x="tuplecounts"), definition=function (x) { dim(x@counts) } )
setMethod("dimnames", signature=c(x="tuplecounts"), definition=function (x) { dimnames(x@counts) } )
setMethod("dimnames<-", signature=c(x="tuplecounts",value="ANY"), definition=function (x,value) { dimnames(x@counts)<-value } )
setMethod("as.matrix", signature=c(x="tuplecounts"), definition=function (x) { as.matrix(x@counts) } )
setMethod("as.vector", signature=c(x="tuplecounts"), definition=function (x) { as.vector(x@counts) } )
setMethod("head", signature=c(x="tuplecounts"), definition=function (x) { head(x@counts) } )
setMethod("image", signature=c(x="tuplecounts"), definition=function (x) { image(x@counts) } )
setMethod("%*%", signature=c(x="tuplecounts",y="ANY"), definition=function (x,y) { x@counts %*% y } )
setMethod("%*%", signature=c(x="ANY",y="tuplecounts"), definition=function (x,y) { x %*% y@counts } )
# add more fluff to avoid "Note: method with signature ‘ANY#tuplecounts’ chosen for function ..." messages
setMethod("%*%", signature=c(x="tuplecounts",y="tuplecounts"), definition=function (x,y) { x@counts %*% y@counts } )
setMethod("%*%", signature=c(x="dgCMatrix",y="tuplecounts"), definition=function (x,y) { x %*% y@counts } )
setMethod("%*%", signature=c(x="tuplecounts",y="dgCMatrix"), definition=function (x,y) { x@counts %*% y } )
setMethod("%*%", signature=c(x="dgTMatrix",y="tuplecounts"), definition=function (x,y) { x %*% y@counts } )
setMethod("%*%", signature=c(x="tuplecounts",y="dgTMatrix"), definition=function (x,y) { x@counts %*% y } )
setMethod("%*%", signature=c(x="TsparseMatrix",y="tuplecounts"), definition=function (x,y) { x %*% y@counts } )
setMethod("%*%", signature=c(x="tuplecounts",y="TsparseMatrix"), definition=function (x,y) { x@counts %*% y } )

setClass("contextModel",
        # A model of mutation; one for each branch.
         representation(
                       genmatrix="genmatrix",
                       projmatrix="Matrix",
                       mutrates="numeric",
                       selcoef="numeric",
                       params="numeric"
                    )
            )

setClass("context",
         representation(
                       counts="tuplecounts",
                       likfun="function",
                       results="list",
                       invocation="character"
                    ),
            contains="contextModel"
         )

setClass("contextMCMC", representation(
                       mutrates.prior="numeric",
                       selcoef.prior="numeric",
                       fixfn.params.prior="numeric"
                   ),
         contains="context")

setOldClass("phylo")  # phylo in ape is an S3 class
setClass("contextTree", representation(
                       counts="tuplecounts",
                       tree="phylo",  # as in 'ape'
                       initfreqs="numeric",
                       models="list",  # list of 'contextModel's
                       likfun="function",
                       results="list",
                       invocation="character"
                   )
         )


setMethod("dimnames", signature=c(x="context"), definition=function (x) { dimnames(x@counts) } )
setMethod("dimnames", signature=c(x="contextTree"), definition=function (x) { dimnames(x@counts) } )
setMethod("counts", signature=c(x="context"), definition=function (x) { counts(x@counts) } )
setMethod("counts", signature=c(x="contextTree"), definition=function (x) { counts(x@counts) } )
setMethod("countframe", signature=c(x="context"), definition=function (x) countframe(x@counts))
setMethod("countframe", signature=c(x="contextTree"), definition=function (x) countframe(x@counts))

# Currently doing this a different way, assigning to environment(x@likfun) to allow calling it directly.
# # We would like the likfun function to automatically have access to the stuff in the context object
# #  ... where is 'self'?!?
# # SUB-OPTIMAL:
# setGeneric("likfun", function (x) { standardGeneric("likfun") } )
# setMethod("likfun", signature=c(x="context"), definition=function(x) {
#           f <- x@likfun
#           environment(f) <- list2env( lapply( selfname(slotNames(model)), function (n) { slot(model,n)} ), parent=globalenv() )
#           return(f) } )

# extract window lengths from these objects:
setGeneric("longwin", function(x) { standardGeneric("longwin") })
setGeneric("shortwin", function(x) { standardGeneric("shortwin") })
setGeneric("leftwin", function(x) { standardGeneric("leftwin") })
setMethod("longwin", signature=c(x="genmatrix"), definition=function(x) { nchar(rownames(x)[1]) } )
setMethod("longwin", signature=c(x="tuplecounts"), definition=function(x) { nchar(rownames(x@counts)[1]) } )
setMethod("longwin", signature=c(x="context"), definition=function(x) { longwin(x@counts) } )
setMethod("longwin", signature=c(x="contextTree"), definition=function(x) { longwin(x@counts) } )
setMethod("shortwin", signature=c(x="tuplecounts"), definition=function(x) { unique(sapply(lapply(lapply(colpatterns(x),levels),"[",1),nchar)) } )
setMethod("shortwin", signature=c(x="context"), definition=function(x) { shortwin(x@counts) } )
setMethod("shortwin", signature=c(x="contextTree"), definition=function(x) { shortwin(x@counts) } )
setMethod("leftwin", signature=c(x="tuplecounts"), definition=function(x) { x@leftwin } )
setMethod("leftwin", signature=c(x="context"), definition=function(x) { leftwin(x@counts) } )
setMethod("leftwin", signature=c(x="contextTree"), definition=function(x) { leftwin(x@counts) } )
# convenience functions
setGeneric("nmuts", function (x) { standardGeneric("nmuts") } )
setGeneric("nsel", function (x) { standardGeneric("nsel") } )
setGeneric("fixparams", function (x) { standardGeneric("fixparams") } )
setMethod("nmuts", signature=c(x="genmatrix"), definition=function (x) { length(x@mutpats) } )
setMethod("nsel", signature=c(x="genmatrix"), definition=function (x) { length(x@selpats) } )
setMethod("fixparams", signature=c(x="genmatrix"), definition=function (x) { (setdiff(names(as.list(formals(x@fixfn))),"..."))[-1] } )
setMethod("nmuts", signature=c(x="contextModel"), definition=function (x) { nmuts(x@genmatrix) } )
setMethod("nsel", signature=c(x="contextModel"), definition=function (x) { nsel(x@genmatrix) } )
setMethod("fixparams", signature=c(x="contextModel"), definition=function (x) { fixparams(x@genmatrix) } )
setMethod("nmuts", signature=c(x="contextTree"), definition=function (x) { sapply( x@models, nmuts ) } )
setMethod("nsel", signature=c(x="contextTree"), definition=function (x) { sapply( x@models, nsel ) } )
setMethod("fixparams", signature=c(x="contextTree"), definition=function (x) { lapply( x@models, fixparams ) } )

# and methods related to model fitting
setMethod("coef", signature=c(object="contextModel"), definition=function (object) {
          coef <- c( object@mutrates, object@selcoef, object@params )
          names(coef) <- c( mutnames( object@genmatrix@mutpats ), selnames( object@genmatrix@selpats ), names(object@params) )
          return(coef) } )
setMethod("coef", signature=c(object="contextTree"), definition=function (object) { 
            tlens <- object@tree$edge.length
            names(tlens) <- paste("tlen",edge.labels(object@tree),sep='.')
            c( tlens, do.call( c, lapply( object@models, coef ) ) ) 
        } )
setMethod("rowSums", signature=c(x="tuplecounts"), definition=function (x) { rowSums(x@counts) } )
setMethod("image", signature=c(x="tuplecounts"), definition=function (x) { image(x@counts) } )
setMethod("rowSums", signature=c(x="context"), definition=function (x) { rowSums(x@counts@counts) } )
setMethod("rowSums", signature=c(x="contextTree"), definition=function (x) { rowSums(x@counts@counts) } )
setMethod("fitted", signature=c(object="context"), definition=function (object,...) { predictcounts.context(object,...) } )
setMethod("residuals", signature=c(object="context"), definition=function (object,...) { resid.context(object,...) } )

## This makes everything go all to hell, for some reason:
# setMethod("-",signature=c("tuplecounts","tuplecounts"), definition=function(e1,e2) { e1@counts <- (e1@counts-e2@counts); return(e1) } )
# setMethod("-",signature=c("tuplecounts","ANY"), definition=function(e1,e2) { new("tuplecounts",leftwin=e1@leftwin,counts=Matrix(e1@counts-e2),bases=e1@bases) } )
# setMethod("-",signature=c("ANY","tuplecounts"), definition=function(e1,e2) { new("tuplecounts",leftwin=e2@leftwin,counts=Matrix(e1-e2@counts),bases=e2@bases) } )


####
# Here is what happens in makegenmatrix:
#
# Let G be the generator matrix, i.e. G[i,j] is the instantaneous transition rate from patterns[i] to patterns[j].
# Construct this as follows:
#  For each mutation pattern u, let I(u) be the set of index pairs (i,j) such that patterns[i] -> patterns[j] can be obtained by applying u,
#    ordered in some way, so I(u) = ( (i_1,j_1), (i_2,j_2), ..., (i_{n(u)},j_{n(u)}) ).
#  Then let I_k be the concatenation of these lists for all mutation patterns in mutpats[[k]].
#  Also, concatenate the I_k to form I.
#  Note that the same pair (i,j) can occur in some I_k twice, or in different I_k,
#    for instance if { C -> T } occurs at rate mutrates[1], and both { CG -> TG, CCG -> CTG } occur at rate mutrates[2],
#    and patterns[1] == CCG, and patterns[23] == CTG, then the pair (1,23) would occur once in I_1 and twice in I_2.
#  Now,
#    G[i,j] = \sum_k mutrates[k] * ( number of times (i,j) appears in I_k ) .
#  Let (x_1, ..., x_N) be the nonzero entries of G, in some order.
#    and let f_m = (i_m,j_m) be the row,column index of x_m in G.
#  We need to compute P so that x = P %*% mutrates, i.e. the N x length(mutpats) matrix
#    P[m,k] = ( number of times f_m appears in I_k )
#  To do this, also note that if we define the N x length(I) matrix
#    Q[m,l] = 1 if the l-th element of I is equal to f_m, and 0 otherwise
#  and the the length(I) x length(mutpats) matrix
#    J[l,k] = 1 if the l-th element of I came from I_k, and 0 otherwise,
#  then
#    P = Q %*% J .

makegenmatrix <- function (mutpats, selpats=list(), patlen=nchar(patterns[1]), patterns=getpatterns(patlen,bases), bases, fixfn, mutrates=rep(1,length(mutpats)), selcoef=rep(1,length(selpats)), selfactors=lapply(selpats,sapply,function(x)1), boundary="none", ...) {
    # Make the generator matrix G on the specified set of patterns,
    # carrying with it the means to quickly update itself.
    #  DON'T do the diagonal (i.e. it must be left empty), so that the updating is easier.
    #
    # Works with column-oriented sparse matrices (see ?"dgCMatrix-class" and ?"CsparseMatrix-class"),
    # for which the following vectors implicitly index the vector @x that actually stores the entries:
    #   @i ; (0-based) row numbers of nonzero entries
    #   @p : (0-based) indices of the first nonzero entries in each column (so diff(x@p) gives the number of nonzero entries by column)
    if (!is.numeric(patlen)|(missing(patlen)&missing(patterns))) { stop("need patlen or patterns") }
    if ( (length(selpats)>0 && max(sapply(unlist(selpats),nchar))>patlen) | max(sapply(unlist(mutpats),nchar))>patlen ) { stop("some patterns longer than patlen") }
    # mutmats is a list of matrices, with one matrix for each of the mutpatls in mutpatll describing the induced mutation process on the specified patterns (see getmutmats function def'n)
    mutmats <- getmutmats(mutpats,patterns,boundary=boundary)
    # concatenate these matrices into one big matrix, allmutmats;
    #   also compute nmutswitches to keep track of which rows of allmutmats correspond to which of mutpats.
    # the row,column indices of allmutmats are exactly the nonzero entries of the resulting generator matrix,
    #   with the first nmutswitches[1] obtainable by mutpats[[1]], the next nmutswitches[2] obtainable by mutpats[[2]], etcetera.
    #   In other words, nmutswitches[k] gives the number of changes that mutpats[[k]] can induce.
    # note that some rows in allmutmats may be identical:
    #   this means that there is more than one way to apply mutation patterns to change patterns[i] to patterns[j].
    #   We correct for this later.
    allmutmats <- do.call( rbind, mutmats )
    nmutswitches <- sapply(mutmats,NROW)
    # allmutmats is the concatenation of all nonzero entries of the generator matrix,
    #   but they come ordered by which mutpat they correspond to.
    # dgCord is the order in which they will appear in the @x vector of the sparse column-oriented generator matrix,
    #   i.e. ordered by column first, then row. Thus, any identical entries will appear sequentially in dgCord.
    # If there were no duplicated rows in allmutmats, then the nonzero entries of the generator matrix would be
    #    @x == mutrates[ rep(1:length(mutpats),each=nmutswitches) ][ dgCord ]
    dgCord <- order( allmutmats$j, allmutmats$i )
    # Now compute helper matrices:
    #   without selection, each nonzero entry in the generator matrix is a linear combination of mutrates;
    #   muttrans is this linear transformation,
    #   i.e. that @x == muttrans %*% mutrates .
    # This is easily constructed: in the order of allmutmats, this is a nrow(allmutmats) x length(mutpats) matrix,
    #   with 1's along a "fat, irregular, diagonal", the k-th row having nmutswitches[k] 1's, and the rest 0's.
    # We then use dgCord to permute rows.
    # In other words, muttrans[i,j] = 1 if mutpat[[j]] can produce the patterns change of allmutmats[i,].
    # ( Construct via TMatrix, which stores row-and-column indices via @i and @j, then convert to CMatrix. )
    muttrans <- dgTtodgC( new( "dgTMatrix", i=1:sum(nmutswitches)-1L, j=rep(seq_along(mutrates),times=nmutswitches)-1L, x=rep(1,sum(nmutswitches)), Dim=c(sum(nmutswitches),length(mutrates)) ) )
    muttrans <- muttrans[ dgCord , , drop=FALSE ]
    # selection is similar to mutation;
    #   we have to compute differences in numbers of matches of the selection patterns, weighted by selfactors
    if (length(selpats)>0) {
        selmatches <- getselmatches( selpats, patterns, selfactors, boundary=boundary )
        # transfer selection coefficients to selective differences involved in each mutation
        #    these are ( transitions ) x ( selpats ) matrix
        ##  the following has (fromsel-tosel); combined to reduce memory usage
        # Made this sparse in case selpats is large
        ## alternative 1:
        selmatches <- t(selmatches)
        fromsel <- selmatches[ allmutmats$i, , drop=FALSE ]
        tosel <- selmatches[  allmutmats$j, , drop=FALSE ]
        seltrans <- ( tosel - fromsel )
        seltrans <- seltrans[dgCord,,drop=FALSE]
        ## alternative 2:
        # seltrans <- ( selmatches[,  allmutmats$j, drop=FALSE ] - selmatches[,  allmutmats$i, drop=FALSE ] )
        # seltrans <- seltrans[,dgCord,drop=FALSE]
        # seltrans <- Matrix(t( seltrans ))
    } else {
        seltrans <- Matrix(numeric(0),nrow=nrow(muttrans),ncol=0)
    }
    # Now deal with the fact there may be duplicated rows in allmutmats:
    #   dups[k] is TRUE if the k-th row of allmutmats[dgCord,] is equal to the (k-1)-th row (see def'n of dgCord),
    #   meaning that they correspond to the same mutpat.
    # We then translate dups to a projection matrix which projects into the unique-mutpat space.
    dups <- c( FALSE, (diff(allmutmats[dgCord,,drop=FALSE][,1])==0) & (diff(allmutmats[dgCord,,drop=FALSE][,2])==0) )
    dupproj <- new("dgTMatrix",i=cumsum(!dups)-1L,j=seq_along(dups)-1L,x=rep(1,length(dups)),Dim=c(nrow(allmutmats)-sum(dups),nrow(allmutmats)))
    muttrans <- dupproj %*% muttrans
    #  and we want to sum mutation rates but not selection probabilities
    #     ... and note that identical rows in allmutmats have identical seltrans entries
    # re-use dupproj: now takes the mean of duplicated entries
    dupproj@x <- as.vector(1/table(dupproj@i))[dupproj@i+1]
    seltrans <- dupproj %*% seltrans
    # Construct the full instantaneous mutation and transition matrix.
    # The `with` function here opens the allmutmats namespace, defining vectors
    #   `i` and `j`, such that an expression like `i=(i-1L)[!dups]` is not
    #   recursive, but rather delivering the new `i` in terms of the allmutmats' `i`.
    genmatrix <- with( allmutmats[dgCord,,drop=FALSE], new( "genmatrix",
            i=(i-1L)[!dups], # Only include a given index if it's not a duplicate.
            p=sapply(0:length(patterns), function (k) sum(j[!dups]<k+1)),
            x=rep(1,sum(!dups)), # Placeholder to get calculated by update.
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
    return(genmatrix)
}

update <- function (G, mutrates, selcoef, ...) {
    # Calculate the new entries of `x` from muttrans et al.
    # use like: genmatrix@x <- update(genmatrix,...)
    fixprob <- if (length(selcoef)>0) { G@fixfn( as.vector(G@seltrans%*%selcoef), ... ) } else { 1 }
    as.vector( G@muttrans %*% mutrates ) * fixprob
}

collapsepatmatrix <- function (ipatterns, leftwin, shortwin=nchar(fpatterns[1]), rightwin=nchar(ipatterns[1])-shortwin-leftwin, fpatterns=getpatterns(nchar(ipatterns[1])-leftwin-rightwin,bases), bases ) {
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
    projmat <- collapsepatmatrix(ipatterns=rownames(genmat),leftwin=leftwin,rightwin=rightwin, bases=genmat@bases)  # this is P
    # Divide entries by column sums, then transpose. Makes M, which is now has rows indexed by
    # the short version and columns the usual one.
    meanmat <- t( sweep( projmat, 2, colSums(projmat), "/" ) )
    pgenmat <- meanmat %*% genmat %*% projmat   # this is H = M G P
    # ii gives the (zero-based) row indices of nonzero entries of pgenmat:
    ii <- pgenmat@i
    # jj gives the (zero-based) column indices of nonzero entries of pgenmat:
    #   see ?"Csparsematrix-class" for explanation of @p.
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
    meangenmat <- new( "genmatrix",
            i=ii[nondiag],
            p=pp,
            x=pgenmat@x[nondiag],
            Dim=pgenmat@Dim,
            Dimnames=pgenmat@Dimnames,
            muttrans = (pnonz %*% genmat@muttrans),
            seltrans = (pnonz %*% genmat@seltrans),
            bases=genmat@bases,
            mutpats=genmat@mutpats,
            selpats=genmat@selpats,
            boundary=genmat@boundary,
            fixfn=genmat@fixfn
        )
    args <- list(...)
    if (is.null(args$mutrates)) { args$mutrates <- rep(1,length(genmat@mutpats)) }
    if (is.null(args$selcoef)) { args$selcoef <- rep(1,length(genmat@selpats)) }
    meangenmat@x <- do.call(update,c( list(G=meangenmat), args ) )
    return( meangenmat )
}

computetransmatrix <- function( genmatrix, projmatrix, tlen=1, shape=1, names=FALSE, transpose=FALSE, time="fixed",...) {
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
        subtransmatrix <- sapply( 1:ncol(projmatrix), function (k) { expAtv( A=A, t=tlen*scale.t, v=projmatrix[,k], ... )$eAtv } )
    }
    if (names) {
        rownames(subtransmatrix) <- rownames(genmatrix)
        colnames(subtransmatrix) <- colnames(projmatrix)
    }
    return( subtransmatrix )
}

product.index <- function ( longpats, bases ) {
    # which base is at each position in each pattern
    apply( do.call(rbind, strsplit(longpats,'') ), 2, match, bases )  
}

get.root.distrn <- function ( initfreqs, initfreq.index ) {
    # given initial frequencies find the distribution over longpats
    #  (the product measure)
    # with helper matrix computed in product.index()
    patfreqs <- initfreqs[initfreq.index]
    dim(patfreqs) <- dim(initfreq.index)
    return( apply( patfreqs, 1, prod ) )
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

cherry.transmats <- function (m1,m2,do.names=FALSE) {
    mm <- m1[,rep(1:ncol(m1),ncol(m2))] * m2[,rep(1:ncol(m2),each=ncol(m1))]
    if (do.names) {
        stopifnot( all( rownames(m1) == rownames(m2) ) )
        rownames(mm) <- rownames(m1)
        colnames(mm) <- paste( colnames(m1)[rep(1:ncol(m1),ncol(m2))], colnames(m2)[rep(1:ncol(m2),each=ncol(m1))], sep=',' )
    }
    return(mm)
}

peel.transmat <- function (tree, rowtaxon, coltaxa, models, genmatrices, projmatrix, root.distrn,
    tlens=tree$edge.length,
    return.list=TRUE, debug=FALSE) {
    ###
    # Compute probabilities of all combinations, on a tree:
    #    rowtaxon is the name of a single tip, or the root
    #    coltaxa is a vector of names of tips
    #    models is a list of model names, in (tip,node) order
    #    genmatrices is a list of genmatrix'es whose names are model names
    ###
    # indices of nodes:
    setup <- peel.transmat.setup(tree, rowtaxon, coltaxa, models, genmatrices, projmatrix, tlens, debug=debug)
    setup <- peel.transmat.compute( setup, models, genmatrices, root.distrn, tlens, debug=debug )
    if (return.list) {
        # everything
        return( setup )
    } else {
        # just the answer
        return( setup$transmats[[ setup$row.node ]] )
    }
}

peel.transmat.setup <- function (tree, rowtaxon, coltaxa, models, genmatrices, projmatrix, 
    tlens=tree$edge.length,
    return.list=FALSE, debug=FALSE) {
    ###
    # Compute probabilities of all combinations, on a tree:
    #    rowtaxon is the name of a single tip, or the root
    #    coltaxa is a vector of names of tips
    #    models is a list of model names, in (tip,node) order
    #    genmatrices is a list of genmatrix'es whose names are model names
    ###
    # indices of nodes:
    root.node <- get.root(tree)
    row.node <- match( rowtaxon, nodenames(tree) )
    col.nodes <- match( coltaxa, nodenames(tree) )
    # reorder branch lengths to match terminal node
    tlens.ord <- match(seq_along(nodenames(tree)),tree$edge[,2])
    # long x short transition matrices
    transmats <- lapply( nodenames(tree), function (x) NULL )
    for (x in col.nodes) { transmats[[x]] <- projmatrix }
    # find path from root to rowtaxon:
    downpath <- c(row.node)
    while( ! root.node %in% downpath ) { downpath <- c( get.parent(downpath[1],tree), downpath ) }
    # list of active nodes that come off of the path from root to rowtaxon:
    up.twigs <- col.nodes[ get.parent(col.nodes,tree) %in% downpath ]
    # list of other active nodes
    up.active <- setdiff( col.nodes, up.twigs )
    ## these two are for bookkeeping;
    # .resolve.up, .resolve.down, .resolve.root, and .resolve.final each store five things:
    #    ( index of first source in transmats, index of second source in transmats, 
    #           index of first genmatrix in genmatrices, index of second genmatix in genmatrices, 
    #           index of first tlen in tlens, index of second tlen in tlens,
    #           index of destination in transmats )
    # this function finds the indices of generator matrices and the indices of tlens
    add.gm <- function (kk) { c( match( models[ nodenames(tree)[kk] ], names(genmatrices) ), tlens.ord[kk] ) }
    # save the order cherries are resolved in
    .resolve.up <- matrix(nrow=0,ncol=7)
    # list of ordering of taxa for columns of each matrix
    col.order <- lapply( nodenames(tree), function (x) c() )
    for (x in col.nodes) { col.order[[x]] <- c(x) }
    # deal with all the upwards branches
    while (length(up.active)>0) {
        up.cherries <- get.cherries(up.active,tree)
        uppair <- up.cherries[1,]
        # compute and combine transition matrices across each branch
        if (debug) cat( "up: ", uppair, "\n" )
        # transmat1 <- computetransmatrix( genmatrices[[ models[ nodenames(tree)[uppair[1]] ] ]], transmats[[uppair[1]]], tlen=tlens[tlens.ord[uppair[1]]], time="fixed", names=debug)
        # transmat2 <- computetransmatrix( genmatrices[[ models[ nodenames(tree)[uppair[2]] ] ]], transmats[[uppair[2]]], tlen=tlens[tlens.ord[uppair[2]]], time="fixed", names=debug)
        # remove cherry
        up.active <- setdiff( up.active, uppair )
        newly.active <- get.parent( uppair[1], tree )
        # transmats[[newly.active]] <- cherry.transmats( transmat1, transmat2, do.names=debug )
        col.order[[newly.active]] <- c( col.order[[ uppair[1] ]], col.order[[ uppair[2] ]] )
        .resolve.up <- rbind( .resolve.up, c( uppair, add.gm(uppair), newly.active ) )
        if (get.parent(newly.active,tree) %in% downpath) {
            up.twigs <- c( up.twigs, newly.active )
        } else {
            up.active <- c( up.active, newly.active )
        }
    }
    # reorder up.twigs to match downpath
    up.twigs <- up.twigs[ match( get.parent(up.twigs,tree), downpath ) ]
    # now move back down from the root (=downpath[1]) to rowtaxon
    # transmats[[ downpath[1] ]] <- sweep( 
    #         computetransmatrix( genmatrices[[ models[nodenames(tree)[ up.twigs[1] ]] ]], transmats[[ up.twigs[1] ]], tlen=tlens[ tlens.ord[up.twigs[1]] ], time="fixed", names=debug ), 
    #     1, root.distrn, "*" )
    col.order[[ downpath[1] ]] <- col.order[[ up.twigs[1] ]]
    # save ordering things are resolved in
    .resolve.root <- c( NA, up.twigs[1], add.gm(c(NA,up.twigs[1])), downpath[1] )
    .resolve.down <- .resolve.final <- matrix(nrow=0,ncol=7)
    # now move back down
    if (length(downpath)>1) {
        for (k in seq_along(up.twigs)[-1]) {
            if (debug) cat("down: ", downpath[k], up.twigs[k], "\n" )
        #     transmat1 <- computetransmatrix( genmatrices[[ models[ nodenames(tree)[downpath[k]] ] ]], transmats[[downpath[k-1]]], tlen=tlens[ tlens.ord[downpath[k]] ], transpose=TRUE, time="fixed", names=debug )
        #     transmat2 <- computetransmatrix( genmatrices[[ models[ nodenames(tree)[up.twigs[k]] ] ]], transmats[[up.twigs[k]]], tlen=tlens[ tlens.ord[up.twigs[k]] ], time="fixed", names=debug )
        #     transmats[[ downpath[k] ]] <- cherry.transmats( transmat1, transmat2, do.names=debug )
            col.order[[ downpath[k] ]] <- c( col.order[[ downpath[k-1] ]], col.order[[ up.twigs[k] ]] )
            .resolve.down <- rbind( .resolve.down, c( downpath[k-1], up.twigs[k], add.gm(c(downpath[k],up.twigs[k])), downpath[k] ) )
        }
        # and finish off
        if (debug) cat("final: ", row.node, "\n")
        # transmats[[row.node]] <- computetransmatrix( genmatrices[[ models[nodenames(tree)[row.node]] ]],
        #         transmats[[downpath[length(downpath)-1]]], tlen=tlens[ tlens.ord[row.node] ], transpose=TRUE, time="fixed", names=debug )
        col.order[[ row.node ]] <- col.order[[ downpath[length(downpath)-1] ]]
        .resolve.final <- c( downpath[length(downpath)-1], NA, add.gm(c(row.node,NA)), row.node)
        if (debug) cat("order: ", paste( col.order[[ row.node ]], sep=',' ), "\n" )
    }
    return( list( transmats=transmats, 
                    up=.resolve.up, root=.resolve.root, down=.resolve.down, final=.resolve.final, 
                    tlens=tlens, tlens.ord=tlens.ord, col.order=col.order,
                    row.node=row.node  ) 
              )
}

peel.transmat.compute <- function (setup, models, genmatrices, root.distrn, tlens, debug=FALSE, return.list=TRUE ) {
    # do the computation using stuff already set up
    # deal with all the upwards branches
    setup$tlens <- tlens
    for (k.up in seq(1,length.out=nrow(setup$up))) {
        upinds <- setup$up[k.up,]
        transmat1 <- computetransmatrix( genmatrices[[ upinds[3] ]], setup$transmats[[upinds[1]]], 
                                        tlen=setup$tlens[upinds[5]], time="fixed", names=debug)
        transmat2 <- computetransmatrix( genmatrices[[ upinds[4] ]], setup$transmats[[upinds[2]]], 
                                        tlen=setup$tlens[upinds[6]], time="fixed", names=debug)
        setup$transmats[[upinds[7]]] <- cherry.transmats( transmat1, transmat2, do.names=debug )
    }
    # now move back down from the root (=downpath[1]) to rowtaxon
    setup$transmats[[ setup$root[7] ]] <- sweep( 
            computetransmatrix( genmatrices[[ setup$root[4] ]], setup$transmats[[ setup$root[2] ]], 
                               tlen=setup$tlens[ setup$root[6] ], time="fixed", names=debug )
            , 1, root.distrn, "*" )
    # now move back down
    for (k.down in seq(1,length.out=nrow(setup$down))) {
        downinds <- setup$down[k.down,]
        transmat1 <- computetransmatrix( genmatrices[[ downinds[3] ]], setup$transmats[[downinds[1]]], 
                                        tlen=setup$tlens[downinds[5]], transpose=TRUE, time="fixed", names=debug )
        transmat2 <- computetransmatrix( genmatrices[[ downinds[4] ]], setup$transmats[[downinds[2]]], 
                                        tlen=setup$tlens[downinds[6]], time="fixed", names=debug )
        setup$transmats[[ downinds[7] ]] <- cherry.transmats( transmat1, transmat2, do.names=debug )
    }
    # and finish off
    finalinds <- setup$final
    if (length(finalinds)>0) {
        setup$transmats[[ finalinds[7] ]] <- computetransmatrix( genmatrices[[ finalinds[3] ]],
                setup$transmats[[finalinds[1]]], tlen=setup$tlens[finalinds[5]], transpose=TRUE, time="fixed", names=debug )
    }
    if (return.list) {
        return( setup )
    } else {
        return( setup$transmats[[ setup$row.node ]] )
    }
}

reorder.counts <- function (counts, new.ord) {
    # Rearrange the columns of counts to a new ordering of taxa
    # new.ord should be a vector of taxon names, i.e. a reordering of colnames(counts@colpatterns)
    old.ord <- colnames(counts@colpatterns) 
    if (! setequal( new.ord, old.ord ) ) { stop("Cannot reorder to a new set of taxa.") }
    oldpats <- apply( sapply( 1:ncol(counts@colpatterns), function (k) {
                      paste( colnames(counts@colpatterns)[k], levels(counts@colpatterns[[k]])[counts@colpatterns[[k]]], sep=":" )
                } ), 1, paste, collapse=',' )
    colnames(counts@colpatterns) <- new.ord
    newpats <- apply( sapply( 1:ncol(counts@colpatterns), function (k) {
                      paste( colnames(counts@colpatterns)[k], levels(counts@colpatterns[[k]])[counts@colpatterns[[k]]], sep=":" )
                } )[,match(new.ord,old.ord),drop=FALSE], 1, paste, collapse=',' )
    counts@counts <- counts@counts[,match(newpats,oldpats),drop=FALSE]
    colnames(counts@counts) <- apply(counts@colpatterns, 1, paste, collapse='.')
    return( counts )
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
predictcounts.context <- function (model, longwin=NULL, shortwin=NULL, leftwin=NULL, initcounts=rowSums(model), mutrates=model@mutrates, selcoef=model@selcoef, genmatrix=model@genmatrix, projmatrix=model@projmatrix, params=model@params ) {
    # default values not cooperating with S4 methods:
    if (is.null(longwin)) { longwin <- longwin(model) }
    if (is.null(shortwin)) { shortwin <- shortwin(model) }
    if (is.null(leftwin)) { leftwin <- leftwin(model) }
    if (!missing(genmatrix) && missing(projmatrix)) {
        projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, fpatterns=getpatterns(shortwin,genmatrix@bases) )
    }
    predictcounts(longwin,shortwin,leftwin,initcounts,mutrates,selcoef,genmatrix,projmatrix,params)
}


predictcounts <- function (longwin, shortwin, leftwin, initcounts, mutrates, selcoef, genmatrix, projmatrix, params=NULL ) {
    # Compute expected counts of paired patterns:
    #  where the actual computation happens
    rightwin <- longwin-shortwin-leftwin
    if (!missing(mutrates)||!missing(selcoef)||!is.null(params)) { genmatrix@x <- do.call( update, c( list(genmatrix,mutrates=mutrates,selcoef=selcoef), params ) ) }
    if (missing(projmatrix)) { projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin, bases=genmatrix@bases ) }
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE)
    fullcounts <- initcounts * subtransmatrix
    return( new("tuplecounts", 
            leftwin=leftwin, 
            counts=Matrix(fullcounts), 
            bases=genmatrix@bases,
            colpatterns=data.frame(short=colnames(fullcounts))
        ) )
}

projectcounts <- function( counts, new.leftwin, new.shortwin, new.longwin, overlapping=FALSE ) {
    # compute counts for shorter Tmers.
    #   valid ranges for parameters are
    #    (l-lc)^+ <= k < (l+w)-(lc+wc)+(r-rc)^-
    #   where 
    #       l = leftwin, lc = new.leftwin
    #       w = shortwin, wc = new.shortwin
    #       r = rightwin, rc = new.rightwin
    # if the original counts were from overlapping windows, then this will overcount the resulting patterns:
    #    if you slide a window of length L in steps of size 1 then a subwindow of size W
    #      will be seen in ( (L-W) ) big windows;
    #    so we need to divide the counts by the factor 'overcount' below
    #      ... but patterns at the boundary of the sequence will not be overcounted. 
    #     Take the ceiling of the resulting counts to fix these.
    leftwin <- leftwin(counts)
    longwin <- longwin(counts)
    shortwin <- shortwin(counts)
    rightwin <- longwin-shortwin-leftwin
    new.rightwin <- new.longwin-new.shortwin-new.leftwin
    if ( max(0L,leftwin-new.leftwin) > (leftwin+shortwin)-(new.leftwin+new.shortwin)+min(0L,rightwin-new.rightwin) ) {
        stop("unreconcilable windows specified.")
    }
    pcounts <- matrix(0,nrow=npatterns(new.longwin,counts@bases),ncol=npatterns(new.shortwin,counts@bases))
    for (k in max(0L,leftwin-new.leftwin):((leftwin+shortwin)-(new.leftwin+new.shortwin)+min(0L,rightwin-new.rightwin))) {
        lpmat <- collapsepatmatrix( ipatterns=rownames(counts), leftwin=k, rightwin=longwin-(k+new.longwin), bases=counts@bases )
        rpmat <- collapsepatmatrix( ipatterns=colnames(counts), leftwin=k+new.leftwin-leftwin, rightwin=shortwin-(k+new.leftwin-leftwin+new.shortwin), bases=counts@bases )
        pcounts <- pcounts + t(lpmat) %*% counts %*% (rpmat)
    }
    dimnames(pcounts) <- list( colnames(lpmat), colnames(rpmat) )
    if (overlapping) {
        overcount <- sum( 
                ( (0:(longwin-1))+new.longwin <= longwin ) &
                ( (0:(longwin-1))+new.leftwin >= leftwin ) &
                ( (0:(longwin-1))+new.leftwin+new.shortwin <= leftwin+shortwin ) )
        pcounts <- ( pcounts/overcount )
    }
    return( new("tuplecounts", 
            leftwin=new.leftwin, 
            counts=Matrix(pcounts), 
            bases=counts@bases, 
            colpatterns=data.frame(short=colnames(pcounts))
        ) )
}

predicttreecounts <- function (shortwin, leftwin=0, rightwin=0, initcounts, mutrates, selcoef, tlens, genmatrix, projmatrix, initfreqs, patcomp, ... ) {
    # Compute expected counts of paired patterns:
    longwin <- leftwin+shortwin+rightwin
    if (missing(genmatrix)) { genmatrix <- makegenmatrix(patlen=leftwin+shortwin+rightwin,...) }
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

# compute and spit out residuals
computeresids <- function (model, pretty=TRUE, in_longwin=longwin(model), in_shortwin=shortwin(model), in_leftwin=leftwin(model), counts=NULL, genmatrixfile=NULL, overlapping=FALSE) {
    # pull in defaults for things if they are not specified
    # load generator matrix, if needed
    if (!is.null(genmatrixfile)) {
        stopifnot(file.exists(genmatrixfile))
        load(genmatrixfile)  # provides 'genmatrix'
    } else {
        genmatrix <- model@genmatrix
    }

    # get counts
    if (is.null(counts) && ( in_longwin > longwin(model) || in_shortwin > shortwin(model) ) ) {
        stop("If window lengths are longer than fitted model, then need to supply counts.")
    } else if (!is.null(counts)) {
        if ( ( in_longwin > longwin(counts) || in_shortwin > shortwin(counts) ) ) {
            stop("Supplied counts use a window that is too short.")
        }
    } else {
        counts <- model@counts
    }

    expected <- fitted( model, longwin=longwin(model), shortwin=shortwin(model), leftwin=leftwin(model), initcounts=rowSums(counts), genmatrix=genmatrix )

    if ( (in_longwin < longwin(counts)) || (in_shortwin < shortwin(counts)) ) {
        counts <- projectcounts( counts, in_leftwin, in_shortwin, in_longwin, overlapping=overlapping )
        expected <- projectcounts( expected, in_leftwin, in_shortwin, in_longwin, overlapping=overlapping )
    }

    if (pretty) {
        # data frame with columns for long pattern, short pattern, observed, expected, residual, z-score
        residframe <- data.frame( inpat=rownames(counts)[row(counts)],
                            outpat=colnames(counts)[col(counts)],
                            observed=as.vector(counts),
                            expected=as.vector(expected),
                            resid=as.vector(counts)-as.vector(expected),
                            stringsAsFactors=FALSE
                        )
        residframe$z <- residframe$resid/sqrt(as.vector(expected))
        residframe <- residframe[order(residframe$z),]
    } else {
        # concise matrix
        resids <- counts@counts - expected@counts  # residuals( model, counts=counts, genmatrix=genmatrix )
        residframe <- as.vector(resids)
        dim(residframe) <- dim(counts)
        dimnames(residframe) <- dimnames(counts)
    }
    return( invisible(residframe) )
}


resid.context <- function (object, counts=object@counts, genmatrix=object@genmatrix, pretty=FALSE) {
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

cbind.dsparseVector <- function (...) {
    # cbind sparseVectors
    input <- lapply( list(...), as, "dsparseVector" )
    thelength <- unique(sapply(input,length))
    stopifnot( length(thelength)==1 )
    sparseMatrix( 
            x=unlist(lapply(input,slot,"x")), 
            i=unlist(lapply(input,slot,"i")), 
            p=c(0,cumsum(sapply(input,function(x){length(x@x)}))),
            dims=c(thelength,length(input))
        )
}



# number of cores for parallel
getcores <- function (subdir) {
    if ( "parallel" %in% .packages()) {
        cpupipe <- pipe("cat /proc/cpuinfo | grep processor | tail -n 1 | awk '{print $3}'")
        numcores <- 1+as.numeric(scan(cpupipe))
        close(cpupipe)
    } else {
        numcores <- 1
    }
    if ( !missing(subdir) && ( as.numeric(gsub("x","",subdir)) < 50 ) ) {
        numcores <- 1
    }
    return(numcores)
}

