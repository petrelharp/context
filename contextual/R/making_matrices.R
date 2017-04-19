###
# functions to prepare mutation and selection matrices, and build fixation functions
# In this section we start seing mutation pattern list lists. These objects make sense because the outer list
# corresponds to the mutation rates, and the inner list corresponds to a set of mutpats that all have the
# same rate.

#' cbind sparseVectors
cbind.dsparseVector <- function (...) {
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


#' Set Up a List of Mutation Matrices
#'
#' The purpose of this function is to construct a translation table between
#' mutpats and the entries of the generator matrix that can happen as a
#' consequence of those mutpats.
#' Specifically, returns i and j's indexing the generator matrix in terms of the given patterns.
#' that is, given a list of mutation patterns,
#'   which can be either pairs or lists of pairs,
#' return a corresponding list of two-column matrices with (1-based) indices of changes corresponding to mutation patterns
#'   i.e. if (i,j) is a row of output[[k]], then patterns[j] can be obtained from patterns[i]
#'   by performing a substitution from mutpats[[k]] at some location within the string.
#'
getmutmats <- function(mutpats, patterns, boundary=c("none","wrap"),
        do.parallel=("parallel" %in% .packages(all.available=TRUE)),
        numcores=if (do.parallel) { parallel::detectCores() } else { 1 }
        ) {
    boundary <- match.arg(boundary)
    longwin <- nchar(patterns[1])
    mutpats <- lapply( mutpats, function (x) { if (is.list(x)) { x } else { list(x) } } )
    this.lapply <- if ( do.parallel && numcores>1 ) { function (...) { parallel::mclapply( ..., mc.cores=numcores ) } } else { lapply }
    lapply( mutpats, function (y) {  # y is mutpat list
            do.call( rbind, this.lapply(y, function (x) {  # x is mutpat (from, to)
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

#' Find Matches for Selection Patterns
#'
#'   each element gets one parameter
#' selmatches[i,j] is the sum of selfactors[[i]] multiplied by the number of times the corresponding thing in selpat[[i]] matches pattern[j]
#'  ... so if selfactors is all 1's, then :
#'     selmatches[i,j] is number of times anything in selpat[[i]] matches pattern[j]  
#'
getselmatches <- function (selpats, patterns, selfactors, 
        boundary=c("none","wrap"), names=FALSE,
        do.parallel=("parallel" %in% .packages(all.available=TRUE)),
        numcores=if (do.parallel) { parallel::detectCores() } else { 1 }
    ) {
    # selpats can be a vector or a list of vectors,
    boundary <- match.arg(boundary)
    substrfun <- switch( boundary, wrap=wrapsubstr, none=substr )
    patlen <- nchar(patterns[1])
    if (!is.list(selpats)) { selpats <- as.list(selpats) }
    this.lapply <- if ( do.parallel && numcores>1 ) { function (...) { parallel::mclapply( ..., mc.cores=numcores ) } } else { lapply }
    selmatches <- t( do.call( cbind.dsparseVector, lapply(seq_along(selpats), function (i) {   # loop over selpats
            y <- selpats[[i]]
            Matrix::rowSums( do.call( cbind.dsparseVector, this.lapply(seq_along(y), function (j) {  # loop over elements of selpats[[y]]
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



