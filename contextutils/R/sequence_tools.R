
###
# genome-ey things

#' Reverse Complement
#'
#' return index for each mutpat of the reverse-complement mutation pattern (or NA if none)
#' throw error if there are multiple matches
#'
#' @export
reverse.complement <- function (mutpats) {
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


#' Count number of times each pattern appears in character string simseq, cyclic-ized.
#'
#' @export
countpats <- function (patterns, simseq ) {
    patlen <- nchar(patterns[1])
    seqlen <- nchar(simseq)
    simseq <- xscat( simseq, subseq(simseq,1,patlen-1) ) # cyclic-ize
    return( sapply( patterns, countPattern, simseq ) )
}

#' Get list of indices corresponding to site the mutation pattern changes
#'
#' @export
changepos <- function (mutpats) {
    lapply( mutpats, lapply, function (x) { which(do.call("!=",strsplit(x,"")) ) } )
}


