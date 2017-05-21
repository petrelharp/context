
###
# stuff for looking at residuals and finding motifs there

#' Compute Residuals
#'
#' @param model An object inheriting from `context` or `contextTree` having a `fitted` method.
#' @param pretty Whether to return results in a data frame sorted by z-score (otherwise, returns a matrix of the same form as `counts`).
#' @param longwin The `long` window length of Tmers to compute residuals for.
#' @param shortwin The `short` window length of Tmers to compute residuals for.
#' @param leftwin The left overhang of Tmers to compute residuals for.
#' @param counts A tuplecounts object - can be supplied to compute residuals for longer Tmers than the model was fit with.
#' @param ... Additional arguments to be passed to `fitted()`.
#'
#' @return If `pretty` is TRUE, returns a data frame with columns `inpat` (the head of the Tmer), `outpat` (the tail of the Tmer), 
#' `observed`, `expected`, `resid` (observed minus expected), and `z` (residual divided by square root of expected).
#' Otherwise, returns a matrix of the same dimensions as `counts` with rows and columns labeled by the head and tail Tmer patterns.
#'
#' @export
computeresids <- function (model, pretty=TRUE, 
                           longwin=longwin(model), 
                           shortwin=shortwin(model), 
                           leftwin=leftwin(model), 
                           counts=NULL, 
                           ...) {
    # get counts
    if (is.null(counts) && ( longwin > longwin(model) || shortwin > shortwin(model) ) ) {
        stop("If window lengths are longer than fitted model, then need to supply counts and genmatrix.")
    } else if (!is.null(counts)) {
        if ( ( longwin > longwin(counts) || shortwin > shortwin(counts) ) ) {
            stop("Supplied counts use a window that is too short.")
        }
    } else {
        counts <- model@counts
    }

    expected <- fitted( model, longwin=longwin, shortwin=shortwin, leftwin=leftwin, ... )

    if ( (longwin < longwin(counts)) || (shortwin < shortwin(counts)) ) {
        counts <- projectcounts( counts, leftwin, shortwin, longwin )
        # expected <- projectcounts( expected, leftwin, shortwin, longwin )
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
        residframe$z <- residframe$resid/sqrt(as.vector(expected)*longwin)
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


