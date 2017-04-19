
#' Circular Substring
#'
#' Finds the substring treating x as a circular string
#'
#' @export
wrapsubstr <- function (x,start,stop) {
    if (all(nchar(x)==0)) { return("") }
    while( any(nchar(x)<stop) ) {
        x <- paste(x,x,sep='')
    }
    substr(x,start,stop)
}

#' @describeIn wrapsubstr Assign to circular substring
#' @export
"wrapsubstr<-" <- function (x,start,stop,value) {
    if (length(value)<length(x)) { value <- rep(value,length(x)) }
    xlen <- nchar(x)
    k <- rep(1,length(x))
    while( TRUE ) {
        stop <- ifelse( xlen>0, stop - xlen*((start-1)%/%xlen), 0 )
        start <- ifelse( xlen>0, (start-1)%%xlen+1, 0 )
        thisstop <- pmin(xlen,stop)
        substr(x,start,thisstop) <- substr(value,k,k+thisstop-start)
        k <- k+thisstop-start+1
        start <- thisstop+1
        if( all(stop<start) ) { break; }
    }
    return(x)
}


