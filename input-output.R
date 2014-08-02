###
# input and output functions
#
# TODO: move some functions over from the other files
###

read.counts <- function (infile,lwin) {
    count.table <- read.table(infile,header=TRUE,stringsAsFactors=FALSE)
    bases <- sort( unique( unlist( strsplit( count.table$reference, "" ) ) ) )
    winlen <- nchar( count.table$reference[1] )
    win <- nchar( count.table$derived[1] )
    longpats <- getpatterns(winlen,bases)
    shortpats <- getpatterns(win,bases)
    counts <- Matrix(0,nrow=length(longpats),ncol=length(shortpats))
    rownames(counts) <- longpats
    colnames(counts) <- shortpats
    stopifnot( all( count.table$reference %in% rownames(counts) ) & all(count.table$derived %in% colnames(counts)) ) 
    counts[cbind( match(count.table$reference,rownames(counts)), match(count.table$derived,colnames(counts)) )] <- count.table$count
    return( new("tuplecounts", counts=counts, lwin=lwin ) )
}
