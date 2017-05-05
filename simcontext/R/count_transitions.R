

#' Count Transitions (Tmers)
#'
#' `counttrans` counts the number of times `initseq` matches each of
#' `ipatterns` while `finalseq` matches each of `fpatterns` (in a Tmer). 
#' The result has rows indexed by patterns in `initseq` and columns indexed by
#' patterns in `finalseq`.
#' 
#' `counttrans.list` counts the number of times that each pattern in
#' `lpatterns[[k]]` appears in `seqlist[[k]]` (for each combination).
#' The result has rows indexed by patterns in the first sequence, and columns
#' indexed by combinations of patterns in the remaining sequences.
#' 
#' Both count optionally cyclically. If `shift` is
#' nonzero, return a list  with the counts in each of (shift) frames: for
#' instance, if `shift=2`, then counts[[1]] will contain counts of Tmers
#' starting at odd positions, and counts[[2]] will contain counts starting at
#' even positions.
#'
#' @param ipatterns Patterns to count in `initseq`, which will be put on the rows of the resulting matrix.
#' @param fpatterns Patterns to count in `finalseq`, which will be put on the columns of the resulting matrix.
#' @param lpatterns A list of patterns to count, of the same length as `seqlist`.
#' @param seqlist A list of sequences to count transitions between.
#' @param initseq The sequence we count "long" patterns in, for `counttrans`.
#' @param finalseq The sequence we count "short" patterns in, for `counttrans`.
#' @param simseqs As output by simseq() - will obtain initseq and finalseq from here if they are missing.
#' @param leftwin The left-hand offset of `fpatterns` under `ipatterns` to make the Tmers.
#' @param shift If this is positive, will return a list of tuplecounts of length `shift`, one anchored at each offset from 1 to shift.
#' @param cyclic Whether to treat the sequences as cyclic.
#' @param bases The vector of possible bases (only used to pass on to the resulting tuplecounts).
#'
#' @return A tuplecounts object, or a list of them if `shift>0`.
#'
#' @examples
#' counts <- counttrans(ipatterns=getpatterns(3,bases=c("X","O")), fpatterns=getpatterns(2,bases=c("X","O")), 
#'                         initseq="XOOXXXOXOXOXOOOXXO", finalseq="XOXOXOXOXXOXOOOXXO", leftwin=1)
#' counts(counts)
#'
#' count.list <- counttrans(ipatterns=getpatterns(3,bases=c("X","O")), fpatterns=getpatterns(2,bases=c("X","O")), 
#'                         initseq="XOOXXXOXOXOXOOOXXO", finalseq="XOXOXOXOXXOXOOOXXO", leftwin=1, shift=2)
#' counts(count.list[[1]])
#' counts(count.list[[2]])
#'
#' @export
counttrans <- function (ipatterns, fpatterns, 
                        initseq=simseqs[["initseq"]], 
                        finalseq=simseqs[["finalseq"]],
                        simseqs, 
                        leftwin=0, 
                        shift=0, 
                        cyclic=FALSE, 
                        bases=sort(unique(unlist(strsplit(ipatterns,""))))
                    ) {
    patlen <- nchar(ipatterns[1])
    seqlen <- nchar(initseq)
    stopifnot( seqlen == nchar(finalseq) )
    # cyclic-ize
    if (cyclic) {
        initseq <- Biostrings::xscat( initseq, Biostrings::subseq(initseq,1,patlen-1) )
        finalseq <- Biostrings::xscat( finalseq, Biostrings::subseq(finalseq,1,patlen-1) )
    }
    # Ok, count occurrences.  uses bioconductor stuff.
    initmatches <- lapply( ipatterns, Biostrings::matchPattern, initseq )
    finalmatches <- lapply( fpatterns, Biostrings::matchPattern, finalseq )
    counts <- lapply( 1:max(1,shift), function (k) Matrix::Matrix( 0, nrow=length(initmatches), ncol=length(finalmatches), dimnames=list(ipatterns,fpatterns) ) )
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
    counts <- lapply( counts, function (x) new("tuplecounts", counts=x, leftwin=leftwin, bases=bases, 
                                               rowtaxon="long", colpatterns=data.frame(short=colnames(x))) )
    if (shift==0) { counts <- counts[[1]] }
    return(counts)
}

#' Count transitions as the result of simulation on a tree.
#'
#' @param rowpatterns The patterns that index rows.
#' @param rowtaxon The taxon counted in rows.
#' @param colpatterns A data frame with column names equal to the tips counted in columns, and one row of patterns per row of the output.
#' @param simseqs A named list of (initseq, finalseq) pairs; 'finalseq' is used.
#'
#' @seealso counttrans, tuplecounts
#'
#' @return A tuplecounts object.
#'
#' @export
counttrans.tree <- function (
                        rowpatterns, 
                        rowtaxon, 
                        colpatterns,
                        simseqs, 
                        leftwin=0, 
                        shift=0, 
                        cyclic=FALSE, 
                        bases=sort(unique(unlist(strsplit(rowpatterns,""))))
                    ) {
    stopifnot(rowtaxon %in% names(simseqs))
    stopifnot(all(colnames(colpatterns) %in% names(simseqs)))
    patlen <- nchar(rowpatterns[1])
    rowseq <- simseqs[[rowtaxon]]$finalseq
    colseqs <- lapply( colnames(colpatterns), function (x) simseqs[[x]]$finalseq )
    seqlen <- nchar(rowseq)
    stopifnot( seqlen == nchar(rowseq) )
    # cyclic-ize
    if (cyclic) {
        rowseq <- Biostrings::xscat( rowseq, Biostrings::subseq(rowseq,1,patlen-1) )
        for (k in seq_along(colseqs)) {
            colseqs[[k]] <- Biostrings::xscat( colseqs[[k]], Biostrings::subseq(colseqs[[k]],1,patlen-1) )
        }
    }
    # Ok, count occurrences.  uses bioconductor stuff.
    initmatches <- lapply( rowpatterns, Biostrings::matchPattern, rowseq )
    finalmatches <- lapply( seq_along(colseqs), function (k) {
                        lapply( colpatterns[,k], Biostrings::matchPattern, colseqs[[k]] )
                    } )
    cns <- apply(colpatterns, 1, paste, collapse='.')
    counts <- lapply( 1:max(1,shift), function (k) Matrix::Matrix( 0, nrow=length(initmatches), ncol=nrow(colpatterns), dimnames=list(rowpatterns,cns) ) )
    for (x in seq_along(initmatches))  {
        for (y in seq_len(nrow(colpatterns))) {
            xycounts <- Biostrings::start(initmatches[[x]])
            for (k in 1:ncol(colpatterns)) {
                xycounts <- intersect(xycounts, (-leftwin)+Biostrings::start(finalmatches[[k]][[y]]))
            }
            if (length(xycounts)>0) {
                for (k in 0:max(0,shift-1)) {
                    counts[[k+1]][x,y] <- sum( ( ( xycounts+k ) %% max(shift,1) ) == 0 )
                }
            }
        }
    }
    counts <- lapply( counts, function (x) new("tuplecounts", 
                                               counts=x, 
                                               leftwin=leftwin, 
                                               bases=bases, 
                                               rowtaxon=rowtaxon, 
                                               colpatterns=colpatterns) )
    if (shift==0) { counts <- counts[[1]] }
    return(counts)
}


#' @describeIn counttrans Count Transitions in a List of Sequences
#' @export counttrans.list
counttrans.list <- function (lpatterns, 
                             seqlist=lapply(simseqs,"[[","finalseq"), 
                             simseqs, 
                             leftwin, 
                             shift=0, 
                             cyclic=FALSE, 
                             bases=sort(unique(unlist(strsplit(lpatterns[[1]],""))))
                             ) {
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
                counts=Matrix::Matrix(x),
                bases=bases,
                rowtaxon=names(seqlist)[1],
                colpatterns=colpatterns
                ) 
        } )
    if (shift==0) { counts <- counts[[1]] }
    return(counts)
}



