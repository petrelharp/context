
#' Mean Density of Nucleotide Differences
#'
#' Uses stringdist::stringdistmatrix.
#'
#' @param counts A tuplecounts object or a matrix of counts with rows and columns labeled.
#' @param leftwin The length of the left overhang of the Tmers of the counts object.
#'
#' @return The mean density.
#'
#' @export
divergence <- function (counts, leftwin) {
    if (missing(leftwin)) { leftwin <- leftwin(counts) }
    require(stringdist)
    patlen <- nchar(colnames(counts)[1])
    nchanges <- stringdist::stringdistmatrix( substr(rownames(counts),leftwin+1,leftwin+patlen), colnames(counts), method="hamming" )
    if (inherits(counts, "tuplecounts")) {
        counts <- counts@counts
    }
    return( sum(nchanges*counts)/(sum(counts)*patlen) )
}

#' Count Numbers of Mutations
#'
#' Given a contingency table `counts` of Tmers observed in a data set, a collection of
#' mutation patterns, the leftwin, and a list of arguments to `sum`, this function
#' returns a matrix with columns indexed by the mutpats, the first row of which gives counts of how many of the counted mutations could be produced by each of mutpats
#'   and the second of which gives the number of "from" matches of the mutpats (called the "total possible").
#'
#' in other words, for each mutation pattern a -> b
#'  sum the values of counts[u,v] over choices of u,v such that:
#'   (i) 'a' matches 'u' at some position
#'   (ii) in addition, 'b' matches 'v' at the same position;
#' store this in the first row
#'
#' note that if we estimate rates by
#'    r.est <- countmuts(...)[1,]/countmuts(...)[2,]
#' then something like
#'    sum( r.est ) / 4
#' should be close to the mean density of nucleotide changes
#'    divergence(...)
#'
#' @export
countmuts <- function (counts, mutpats, leftwin=leftwin(counts), ...) {
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


