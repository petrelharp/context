#!/usr/bin/env Rscript
library(optparse)
library(jsonlite)

usage <- '\
Pre-compute a generator matrix and associated structures,  \
reading list of patterns from a configuration file, \
formatted in json as follows (see singlebase.mut): \
 \
{  \
    "bases" : [ "A", "C", "G", "T" ], \
    "mutpats" : [ \
    [ [ "A", "T" ], [ "T", "A" ] ], \
    [ [ "C", "G" ], [ "G", "C" ] ], \
    [ [ "A", "C" ], [ "T", "G" ] ], \
    [ [ "A", "G" ], [ "T", "C" ] ], \
    [ [ "C", "A" ], [ "G", "T" ] ], \
    [ [ "C", "T" ], [ "G", "A" ] ] \
    ], \
    "selpats" : [ \
    ], \
    "fixfn" : "null.fixfn"
} \
\
Further entries will be ignored.\
\
Entries in the same sub-list have the same associated parameter (reverse-complement, above). \
\
To make this for e.g. all patterns up to 3-way, with one change each: \
\
library(jsonlite) \
bases <- c( "A", "C", "G", "T") \
pats <- getmutpats(3,bases) \
# check this is idempotent \
stopifnot( identical( pats, fromJSON(toJSON(pats,pretty=TRUE),simplifyMatrix=FALSE) ) ) \
sink("triplebase.mut") \
toJSON(list(bases=bases,mutation=pats,selection=NULL)) \
sink(NULL) \
\
'

option_list <- list(
        make_option( c("-c","--configfile"), type="character", help="Input file of mutation and selection motifs." ),
        make_option( c("-o","--outfile"), type="character", help="Save resulting matrix in this file.  [default: genmatrix-(longwin)-(configfile)-(modelname).RData]" ),
        make_option( c("-w","--longwin"), type="integer", help="Size of matching window." ),
        make_option( c("-b","--boundary"), type="character", default="none", help="Boundary conditions, either 'none' or 'wrap'. [default \"%default\"]" ),
        make_option( c("-m","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases. [default \"%default\"]" ),
        make_option( c("-l","--logfile"), type="character", help="Direct logging output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$configfile) || !file.exists(opt$configfile)) { stop("Need a config file.") }
if (is.null(opt$longwin)) { stop("Need a window length.") }

if (is.null(opt$outfile)) { 
    opt$outfile <- paste(dirname(opt$configfile),"/",paste(c(paste("genmatrix",opt$longwin,gsub(".json$","",basename(opt$configfile)),sep="-"),opt$modelname),collapse="-"),".RData",sep='') 
}
basename <- gsub(".RData",'',opt$outfile)
if (!file.exists(dirname(opt$outfile))) { dir.create(dirname(opt$outfile)) }

library(contextual)
library(contextutils)

config <- read.config(opt$configfile)
if (!is.null(config$tree)) {
    config <- parse.models( treeify.config( config ) )
    modelnames <- config$.models
} else {
    config <- fill.default.config( config )
    config$fixfn <- parse.fixfn( config$fixfn, config$fixfn.params )
    config$genmatrix <- opt$outfile
    config <- list( default=config )
    modelnames <- c("default")
}


for (mm in modelnames) {
    if (opt$meanboundary==0) {
        genmatrix <- do.call( makegenmatrix, c( list( 
                                patlen=opt$longwin, 
                                mutpats=config[[mm]]$mutpats, 
                                selpats=config[[mm]]$selpats, 
                                selfactors=config[[mm]]$selfactors, 
                                boundary=opt$boundary, 
                                bases=config[[mm]]$bases, 
                                fixfn=config[[mm]]$fixfn ), 
                        config[[mm]]$fixfn.params ) )
    } else {
        genmatrix <- do.call( meangenmatrix, c( list( 
                                leftwin=opt$meanboundary, rightwin=opt$meanboundary, patlen=opt$longwin, 
                                mutpats=config[[mm]]$mutpats, 
                                selpats=config[[mm]]$selpats, 
                                selfactors=config[[mm]]$selfactors, 
                                boundary=opt$boundary, 
                                bases=config[[mm]]$bases, 
                                fixfn=config[[mm]]$fixfn ), 
                        config[[mm]]$fixfn.params ) )
    }
    if (!is.null(config[[mm]]$mutrates)) {
        genmatrix@x <- do.call( update_x, c( list(G=genmatrix, mutrates=config[[mm]]$mutrates, selcoef=config[[mm]]$selcoef), config[[mm]]$fixfn.params ) )
    }
    outfile <- gsub("%",opt$longwin,config[[mm]]$genmatrix,fixed=TRUE)
    cat("Writing out to ", outfile, "\n")
    save( genmatrix, file=outfile )
}
