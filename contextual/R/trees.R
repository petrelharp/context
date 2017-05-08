#' Reorder Counts by Taxon Reordering
#'
#' Rearrange the columns of counts to a new ordering of taxa
#'
#' @param counts A tuplecounts object.
#' @param new.ord A vector of taxon names, a reordering of coltaxa(counts)
#'
#' @return A new tuplecounts object.
#'
#' @export reorder.counts
reorder.counts <- function (counts, new.ord) {
    old.ord <- colnames(counts@colpatterns) 
    if (! setequal( new.ord, old.ord ) ) { stop("Cannot reorder to a new set of taxa.") }
    oldpats <- apply( sapply( 1:ncol(counts@colpatterns), function (k) {
                      paste( colnames(counts@colpatterns)[k], as.character(counts@colpatterns[[k]]), sep=":" )
                } ), 1, paste, collapse=',' )
    colnames(counts@colpatterns) <- new.ord
    newpats <- apply( sapply( 1:ncol(counts@colpatterns), function (k) {
                      paste( colnames(counts@colpatterns)[k], as.character(counts@colpatterns[[k]]), sep=":" )
                } )[,match(new.ord,old.ord),drop=FALSE], 1, paste, collapse=',' )
    counts@counts <- counts@counts[,match(newpats,oldpats),drop=FALSE]
    colnames(counts@counts) <- apply(counts@colpatterns, 1, paste, collapse='.')
    return( counts )
}


