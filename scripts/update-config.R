#!/usr/bin/env Rscript
library(optparse)

invocation <- commandArgs()

usage <- "\
Update parameters in a config JSON file from a fitted model.
"

option_list <- list(
    # input/output
        make_option( c("-m","--modelfile"), type="character", help="RData containing model." ),
        make_option( c("-c","--configfile"), type="character", help="JSON containing previous configuration."),
        make_option( c("-o","--outfile"), type="character", help="Where to write output JSON configuration to.")
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$modelfile) || is.null(opt$configfile) || is.null(opt$outfile)) { stop("No input files.  Run\n  update-config.R -h\n for help.\n") }
if (!file.exists(opt$modelfile)) { stop("Could not find model file.") }
if (!file.exists(opt$configfile) ) { stop("Could not find config file `", opt$configfile, "`.") }

library(contextual)
library(jsonlite)
library(ape)

loaded <- load(opt$modelfile)  # provides 'model'
stopifnot("model" %in% loaded)
stopifnot(inherits(model, "contextTree"))

dirlist <- function (a) { da <- dirname(a); if (da == a) { if (a=="/") { NULL } else { a } } else { c( dirlist(da), basename(a) ) } }
relpath <- function (a,b) {
    # returns relative path from location of a to b
    if (file.exists(a) && !dir.exists(a)) { a <- dirname(a) }
    aa <- dirlist(normalizePath(a))
    bb <- dirlist(normalizePath(b))
    k <- 0
    while (k < length(aa) && k < length(bb) && aa[k+1] == bb[k+1]) { k = k + 1 }
    do.call(file.path, c(list("..")[rep(1,length(aa)-k)], as.list(bb[-(1:k)])))
}

# read in config file
raw.config <- fromJSON(opt$configfile, 
                       simplifyMatrix = FALSE, 
                       simplifyDataFrame = FALSE)
config.pointers <- config.dereference(raw.config, names(model@models))
config <- treeify.config(read.config(opt$configfile))

stopifnot(all(config$tree$tip.label == model@tree$tip.label))
stopifnot(all(config$tree$node.label == model@tree$node.label))
raw.config$tree <- write.tree(model@tree)

for (k in seq_along(model@models)) {
    modelname <- names(model@models)[k]
    put_here <- config.pointers[modelname]
    # mut rates
    stopifnot(all(mutnames(raw.config[[put_here]]$mutpats) == names(model@models[[1]]@mutrates)))
    if (length(raw.config[[put_here]]$mutpats)>0) {
        raw.config[[put_here]]$mutrates <- model@models[[1]]@mutrates
    }
    stopifnot(all(selnames(raw.config[[put_here]]$selpats) == names(model@models[[1]]@selcoef)))
    if (length(raw.config[[put_here]]$selcoef)>0) {
        raw.config[[put_here]]$selcoef <- model@models[[1]]@selcoef
    }
    if (length(raw.config[[put_here]]$fixfn.params)>0) {
        raw.config[[put_here]]$fixfn.params <- model@models[[1]]@params
    }
    # would be better to do this by relative paths
    raw.config[[put_here]]$genmatrix <- file.path(relpath(opt$outfile, normalizePath(file.path(dirname(opt$configfile),dirname(raw.config[[put_here]]$genmatrix)))), basename(raw.config[[put_here]]$genmatrix))

}

cat(toJSON(raw.config, pretty=TRUE), file=opt$outfile)
