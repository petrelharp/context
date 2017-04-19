
#' Show Simulated Sequence
#'
#' Display the output of simulation by lowercasing or replacing with "." bases that didn't change
#'
#' @export
show.simseq <- function (x,printit=FALSE,maxchar=min(nchar(x$initseq),200),latex=FALSE) {
    x$initseq <- subseq(x$initseq,1,maxchar)
    fseq <- BString( subseq(x$finalseq,1,maxchar) )
    if (!is.null(x$ntrans)) {
        x$ntrans <- subset(x$ntrans,loc<=maxchar)
        patlen <- nchar( levels( x$ntrans$i )[1] )
        events <- unique(x$ntrans$loc)
        eventlocs <- (seq_len(nchar(fseq)) %in% outer(events,0:(patlen-1),"+"))
        if (length(events)>0) {
            whichchanged <- sapply( events, function (k) { kev <- (x$ntrans$loc==k); any( x$ntrans$i[kev] != x$ntrans$j[kev] ) } )
            changedlocs <- ( seq_len(nchar(fseq)) %in% outer( events[whichchanged], 0:(patlen-1), "+" ) )
            fseq[eventlocs & !changedlocs] <- BString( tolower(fseq[eventlocs & !changedlocs]) )
        }
        fseq[!eventlocs] <- "."
    }
    outstrings <- c(
            paste(as.character(x$initseq),"\n",sep=''),
            paste(as.character(fseq),"\n",sep='')
        )
    if (latex) {
        outstrings <- gsub("&\n&$", "\\\\\\\\\n", gsub( "^&", "", gsub("\\.","\\$\\\\centerdot\\$", gsub("","&", outstrings ) ) ) )
        outstrings <- paste( c(
                paste( "\\begin{center} \\setlength{\\tabcolsep}{0pt} \\begin{tabular}{",
                    paste(rep('c',nchar(outstrings[[1]])),collapse=''), "}\n", sep='' ),
                outstrings,
                "\\end{tabular} \\end{center} \n"
            ), collapse="" )
    }
    if (printit | latex) { lapply( outstrings, cat ) }
    return(invisible(outstrings))
}

#' Print Simulated Sequence
#'
#' @export
print.simseq <- function (x,...) {
    invisible( show.simseq(x,printit=TRUE,maxchar=min(nchar(x$initseq),150)) )
}

