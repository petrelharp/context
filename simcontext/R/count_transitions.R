

#' Count Transitions
#'
#' count number of times initseq matches ipatterns while finalseq matches fpatterns
#'   optionally, cyclical
#' if shift is nonzero, return a list  with the counts in each of (shift) frames
#'
#' @return A tuplecounts object.
#'
#' @export
counttrans <- function (ipatterns, fpatterns, 
                        initseq=simseqs[["initseq"]], finalseq=simseqs[["finalseq"]],
                        simseqs, leftwin=0, shift=0, cyclic=FALSE, bases=sort(unique(unlist(strsplit(ipatterns,""))))) {
    patlen <- nchar(ipatterns[1])
    seqlen <- nchar(initseq)
    stopifnot( seqlen == nchar(finalseq) )
    # cyclic-ize
    if (cyclic) {
        initseq <- xscat( initseq, Biostrings::subseq(initseq,1,patlen-1) )
        finalseq <- xscat( finalseq, Biostrings::subseq(finalseq,1,patlen-1) )
    }
    # Ok, count occurrences.  uses bioconductor stuff.
    initmatches <- lapply( ipatterns, matchPattern, initseq )
    finalmatches <- lapply( fpatterns, matchPattern, finalseq )
    counts <- lapply( 1:max(1,shift), function (k) Matrix( 0, nrow=length(initmatches), ncol=length(finalmatches), dimnames=list(ipatterns,fpatterns) ) )
    for (x in seq_along(initmatches))  {
        for (y in seq_along(finalmatches)) {
            xycounts <- intersect(start(initmatches[[x]]),(-leftwin)+start(finalmatches[[y]]))
            if (length(xycounts)>0) {
                for (k in 0:max(0,shift-1)) {
                    counts[[k+1]][x,y] <- sum( ( ( xycounts+k ) %% max(shift,1) ) == 0 )
                }
            }
        }
    }
    counts <- lapply( counts, function (x) new("tuplecounts",counts=x,leftwin=leftwin,bases=bases,rowtaxon="long",colpatterns=data.frame(short=colnames(x))) )
    if (shift==0) { counts <- counts[[1]] }
    return(counts)
}


#' Count transitions with Shift
#'
#' count number of times each of seqlist matches corresponding lpatterns
#'   optionally, cyclical
#' if shift is nonzero, return a list  with the counts in each of (shift) frames
#'
#' @export
counttrans.list <- function (lpatterns, seqlist=lapply(simseqs,"[[","finalseq"), 
                             simseqs, leftwin, shift=0, cyclic=FALSE, 
                             bases=sort(unique(unlist(strsplit(lpatterns[[1]],""))))) {
    patlen <- max( sapply( lapply( lpatterns, "[", 1 ), nchar ) )
    seqlen <- unique( sapply(seqlist, nchar) )
    stopifnot(length(seqlen)==1)
    if (length(leftwin)<length(seqlist)) { shiftwin <- c(0,rep(leftwin,length.out=length(seqlist)-1)) } else { shiftwin <- leftwin }
    # cyclic-ize
    if (cyclic) { seqlist <- lapply( seqlist, function (x) { xscat( x, Biostrings::subseq(x,1,patlen-1) ) } ) }
    # Ok, count occurrences.  uses bioconductor stuff.
    lmatches <- lapply( seq_along(seqlist), function (k) {
            lapply( lpatterns[[k]], matchPattern, seqlist[[k]] )
        } )
    npats <- sapply(lmatches,length)
    counts <- lapply( 1:max(1,shift), function (k) array( 0, dim=npats, dimnames=lpatterns ) )
    ii <- matrix( rep(1,length(npats)), nrow=1 )
    while( ii[1] <= npats[1] ) {
        xycounts <- start(lmatches[[1]][[ii[1]]]) - shiftwin[1]
        for (j in seq_along(ii)[-1]) {
            stopifnot( j %in% seq_along(lmatches) & ii[j] %in% seq_along(lmatches[[j]]) )
            xycounts <- intersect( xycounts, (-shiftwin[j]) + start(lmatches[[j]][[ii[j]]]) )
        }
        if (length(xycounts)>0) {
            for (k in 0:max(0,shift-1)) {
                    counts[[k+1]][ii] <- sum( ( ( xycounts+k ) %% max(shift,1) ) == 0 )
            }
        }
        ii[length(ii)] <- ii[length(ii)]+1
        for (j in rev(seq_along(ii)[-1])) { ii[j-1] <- ii[j-1] + ( ii[j] > npats[j] ) }
        ii[-1] <- 1 + ( (ii[-1]-1) %% npats[-1] )
    }
    colpatterns <- do.call( expand.grid, lpatterns[-1] )
    colnames(colpatterns) <- names(seqlist)[-1]
    counts <- lapply( counts, function (x) {
            dim(x) <- c(dim(x)[1],prod(dim(x)[-1]))
            rownames(x) <- lpatterns[[1]]
            colnames(x) <- apply(colpatterns,1,paste,collapse='.')
            new("tuplecounts", 
                leftwin=leftwin,
                counts=Matrix(x),
                bases=bases,
                rowtaxon=names(seqlist)[1],
                colpatterns=colpatterns
                ) 
        } )
    if (shift==0) { counts <- counts[[1]] }
    return(counts)
}



