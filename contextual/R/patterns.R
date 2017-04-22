
#' Construct a List of All Patterns of a Given Length
#' @export
getpatterns <- function(patlen,bases) {
    patterns <- do.call( expand.grid, rep( list(bases), patlen ) )
    return( apply(patterns,1,paste,collapse="") )
}

#' Find all possible mutational changes
#'
#' return list of single-base changes for each list of patterns in mutpats
#' a `mutpat` is a 2-element string list of the form (from, to)
#'
#' @export
mutpatchanges <- function (mutpats) {
    lapply( mutpats, function (x) { do.call(rbind, lapply( x, function (y) {
                ysplit <- strsplit(y,'')
                differs <- do.call("!=",ysplit)
                cbind( ysplit[[1]][differs], ysplit[[2]][differs] )
            } ) ) } )
}

#' Get all k-mer changes
#'
#' Finds all kmer -> kmer changes where k <= `patlen`
#' that involve at most nchanges changes.
#'
#' @export
getmutpats <- function(bases,patlen,nchanges=1) {
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

#' Number of Possible Changes
#'
#' @return Integer number of possible patterns.
npatterns <- function (patlen,bases) {
    return( length(bases)^patlen )
}

#' Print a Mutation Pattern List in Readable Form
#'
#' Stringify a list of mutpat lists,
#'  if pretty=TRUE, then using the names if they are available
#'
#' @export
mutnames <- function (mutpats,pretty=TRUE) {
    stringified <- unlist( sapply( lapply( mutpats, lapply, paste, collapse="->" ), paste, collapse="|" ) ) 
    if (pretty && !is.null(names(mutpats))) {
        stringified[nchar(names(mutpats))>0] <- names(mutpats)[nchar(names(mutpats))>0]
    }
    return( stringified )
}

#' Print a Selection Pattern List in Readable Form
#'
#' Stringify a list of selpats,
#'  but use the names if they are available
#'
#' @export
selnames <- function (selpats,pretty=TRUE) {
    stringified <- sapply(selpats,paste,collapse="|")
    if (pretty && !is.null(names(selpats))) {
        stringified[nchar(names(selpats))>0] <- names(selpats)[nchar(names(selpats))>0]
    }
    return( stringified )
}


#' Construct a Pattern-Substring Matrix
#'
#' Construct the matrix U described in the tex.
#' ipatterns are the "input" patterns, while fpatterns are the "final" projected patterns
#' This function assumes that all input patterns are the same length.
#'
#' @param ipatterns Long patterns index the rows of the result.
#' @param leftwin Left offset of short patterns from long patterns.
#' @param shortwin Length of short patterns.
#' @param rightwin Right offset of short patterns from long patterns.
#' @param fpatterns Short patterns, which index the columns of the result.
#' @param bases Alphabet of bases.
#' @return A (nbases)^k by (nbases)^{k-leftwin-rightwin} sparse matrix
#' projection matrix that maps patterns onto the shorter patterns obtained by
#' deleting leftwin characters at the start and rightwin characters at the end.  
#'
#' @examples
#' bases <- c("X","O")
#' # columns are the right-hand character
#' collapsepatmatrix(ipatterns=getpatterns(2,bases), leftwin=1, fpatterns=bases)
#' # columns are the middle two characters
#' collapsepatmatrix(ipatterns=getpatterns(4,bases), leftwin=1, rightwin=1, bases=bases)
#'
#' @export
collapsepatmatrix <- function (ipatterns, 
                               leftwin, 
                               shortwin=nchar(fpatterns[1]), 
                               rightwin=nchar(ipatterns[1])-shortwin-leftwin, 
                               fpatterns=getpatterns(nchar(ipatterns[1])-leftwin-rightwin,bases), 
                               bases ) {
    patlen <- nchar(ipatterns[1])
    shortwin <- patlen - leftwin - rightwin
    stopifnot(shortwin>0)
    matchpats <- match( substr(ipatterns,leftwin+1L,leftwin+shortwin), fpatterns )
    matchmatrix <- new( "dgTMatrix", i=(seq_along(ipatterns)-1L), j=(matchpats-1L), x=rep(1,length(ipatterns)), Dim=c(length(ipatterns),length(fpatterns)) )
    rownames(matchmatrix) <- ipatterns
    colnames(matchmatrix) <- fpatterns
    return( matchmatrix )
}

#' Which base is at each position in each pattern
#'
#' @export
product.index <- function ( longpats, bases ) {
    apply( do.call(rbind, strsplit(longpats,'') ), 2, match, bases )  
}

#' Product Distribution for the Root of a Tree
#'
#' Given initial frequencies find the distribution over longpats
#' (the product measure) with helper matrix computed in product.index()
#'
#' @export
get.root.distrn <- function ( initfreqs, initfreq.index ) {
    patfreqs <- initfreqs[initfreq.index]
    dim(patfreqs) <- dim(initfreq.index)
    return( apply( patfreqs, 1, prod ) )
}


