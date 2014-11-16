#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "Process lots of summary .Rmd files."

option_list <- list(
        make_option( c("-o","--outfile"), type="character", help="File to save things in.")
        )
cmd <- parse_args(OptionParser(option_list=option_list,description=usage), positional_arguments = TRUE)

source("../context-inference-fns.R")

digest <- function(path_complete) {
    load(path_complete)
    c(seed = simseq.opt$seed,
      configfile = simseq.opt$configfile,
      cov = cov(mr_compare$simulated, mr_compare$fit))
}

write.table(
    do.call(rbind,lapply(cmd$args, digest)),
    cmd$options$outfile)
