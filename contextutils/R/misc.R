
#' Name a list or vector with itself.
#' @param x A vector.
#' @return Itself, named by itself.
#' @export
selfname <- function (x) { names(x) <- x; return(x) }


