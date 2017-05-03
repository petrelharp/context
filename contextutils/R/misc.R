
#' Name a list or vector with itself.
#' @param x A vector.
#' @return Itself, named by itself.
#' @export
selfname <- function (x) { names(x) <- x; return(x) }


#' Open for reading or writing
#'
#' Use to open stdin/stdout or process substitution things correctly
#'   from  http://stackoverflow.com/questions/15784373/process-substitution
#'
#' @param arg A text string: one of "-", "[/dev/]stdin", "[/dev/]stdout", or a file name (including a file descriptor, e.g. "/dev/fd3").
#'
#' @return A connection, from one of `stdout()`, `stdin()`, `fifo()`, or `file()`.
#'
#' @export
openread <- function(arg) {
    if (arg %in% c("-", "/dev/stdin","stdin")) {
       stdin()
    } else if (grepl("^/dev/fd/", arg)) {
       fifo(arg, open = "r")
    } else {
       file(arg, open = "r")
    }
}
#' @describeIn openread Open for writing.
#' @export
openwrite <- function(arg) {
    if (arg %in% c("-", "/dev/stdout","stdout")) {
       stdout()
    } else if (grepl("^/dev/fd/", arg)) {
       fifo(arg, open = "w")
    } else {
       file(arg, open = "w")
    }
}

