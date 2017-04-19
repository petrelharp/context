#' Reorder Counts by Taxon Reordering
#'
#' Rearrange the columns of counts to a new ordering of taxa
#' new.ord should be a vector of taxon names, i.e. a reordering of colnames(counts@colpatterns)
#'
#' @export
reorder.counts <- function (counts, new.ord) {
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


