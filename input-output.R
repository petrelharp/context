###
# input and output functions
#
# TODO: move some functions over from the other files
###

read.counts <- function (infile,leftwin,bases,longpats,shortpats) {
    # read in a file of counts of the following form:
    #     reference derived count
    #   1      AAAA      AA     1
    #   2      CAAA      AA     2
    #   3      GAAA      AA     2
    #   4      TAAA      AA     3
    # ... and convert it to a 'tuplecounts' object
    # optionally passing in the orderings of the rows and columns
    count.table <- read.table(infile,header=TRUE,stringsAsFactors=FALSE)
    longwin <- nchar( count.table$reference[1] )
    shortwin <- nchar( count.table$derived[1] )
    if ( missing(bases) ) { bases <- sort( unique( unlist( strsplit( count.table$reference, "" ) ) ) ) }
    if ( missing(longpats) ) { longpats <- getpatterns(longwin,bases) }
    if ( missing(shortpats) ) { shortpats <- getpatterns(shortwin,bases) }
    counts <- Matrix(0,nrow=length(longpats),ncol=length(shortpats))
    rownames(counts) <- longpats
    colnames(counts) <- shortpats
    stopifnot( all( count.table$reference %in% rownames(counts) ) & all(count.table$derived %in% colnames(counts)) ) 
    counts[cbind( match(count.table$reference,rownames(counts)), match(count.table$derived,colnames(counts)) )] <- count.table$count
    return( new("tuplecounts", counts=counts, leftwin=leftwin, bases=bases ) )
}

