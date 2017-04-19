
#' @export
dgTtodgC <- function (M) {
    # convert between matrix classes.  For understanding.
    ijx <- data.frame( i=M@i, j=M@j, x=M@x )
    ijx <- ijx[ order( ijx$j, ijx$i ), ]
    with(ijx, new( "dgCMatrix", i=i, p=sapply(0:ncol(M), function(k) sum(j<k)), x=x, Dim=dim(M) ) )
}


