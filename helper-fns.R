# miscellaneous utility functions

selfname <- function (x) { names(x) <- x; return(x) }

debug.optparse <- function (command,option_list) {
    # use for interactive debugging of optparse scripts
    # In place of e.g.  in the script this-script.R
    #    opt <- parse_args(OptionParser(option_list=option_list,description=usage))
    # use
    #    opt <- debug.optparse("Rscript this-script.R -a 1 -b 2 --blah", option_list)
    args <- strsplit(command," ")[[1]]
    if (args[1]=="Rscript") { args <- args[-(1:2)] }
    return( parse_args(OptionParser(option_list=option_list,description=usage),args=args) )
}

# load a file, but not into the global environment, rather, into a list.
load.to.list <- function (file) { e <- environment(); n <- load(file,envir=e); names(n) <- n; lapply( n, get, envir=e ) }

# use to open stdin/stdout or process substitution things correctly
#   from  http://stackoverflow.com/questions/15784373/process-substitution
openread <- function(arg) {
    if (arg %in% c("-", "/dev/stdin","stdin")) {
       file("stdin", open = "r")
    } else if (grepl("^/dev/fd/", arg)) {
       fifo(arg, open = "r")
    } else {
       file(arg, open = "r")
    }
}
openwrite <- function(arg) {
    if (arg %in% c("-", "/dev/stdout","stdout")) {
       file("stdout", open = "w")
    } else if (grepl("^/dev/fd/", arg)) {
       fifo(arg, open = "w")
    } else {
       file(arg, open = "w")
    }
}
