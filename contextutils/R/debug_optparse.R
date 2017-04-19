#' @export
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

###
# use for debugging noninteractive stuff
# options(error=print.and.dump)

#' @export
print.and.dump <- function () {
    cat(paste("Error in \"", paste(commandArgs(),collapse=' '), "\": dumping frames.\n")); dump.frames(to.file = TRUE); q(status=1)
} 

