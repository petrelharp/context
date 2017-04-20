
###
# stuff for looking at residuals and finding motifs there

#' Compute Residuals
#'
#' @export
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

    expected <- fitted( model, longwin=longwin(model), shortwin=shortwin(model), leftwin=leftwin(model), initcounts=Matrix::rowSums(counts), genmatrix=genmatrix )

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


#' used in setMethod
resid.context <- function (object, counts=object@counts, genmatrix=object@genmatrix, pretty=FALSE) {
    expected <- fitted(object, initcounts=Matrix::rowSums(counts),
                       longwin=longwin(counts), shortwin=shortwin(counts), leftwin=leftwin(counts),
                       genmatrix=genmatrix )
    stopifnot( all( rownames(expected)==rownames(counts) ) && all( colnames(expected)==colnames(counts) ) )
    resids <- counts
    resids@counts <- ( counts@counts - expected@counts )
    return( resids )
}


