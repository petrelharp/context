#!/usr/bin/Rscript --vanilla
require(optparse)
require(jsonlite)

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
    ] \
} \
 \
Entries in the same sub-list have the same associated parameter (reverse-complement, above). \
\
To make this for e.g. all patterns up to 3-way, with one change each: \
\
require(jsonlite) \
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
        make_option( c("-c","--configfile"), type="character", default="-", help="Input file of mutation and selection motifs." ),
        make_option( c("-s","--outfile"), type="character", default="", help="Save resulting matrix in this file.  [default: genmatrix-(winlen)-(configfile).RData]" ),
        make_option( c("-w","--winlen"), type="integer", help="Size of matching window." ),
        make_option( c("-b","--boundary"), type="character", default="none", help="Boundary conditions. [default \"%default\"]" ),
        make_option( c("-m","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases. [default \"%default\"]" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct logging output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$configfile) | !file.exists(opt$configfile)) { stop("Need a config file.") }
if (is.null(opt$winlen)) { stop("Need a window length.") }
options(error=traceback)

config <- fromJSON(opt$configfile,simplifyMatrix=FALSE)
attach(c(opt,config))
# rename (config file is more user-friendly?)
mutpats <- mutation
selpats <- selection

if (outfile=="") { 
    outfile <- paste(dirname(configfile),"/",paste("genmatrix",winlen,gsub(".mut$","",basename(configfile)),sep="-"),".RData",sep='') 
}
basename <- gsub(".RData",'',outfile)
if (logfile=="" & !interactive()) { logfile <- paste(basename,"-flex-genmat.Rout",sep='') }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message") 
    sink(file=logcon, type="output") 
}

if (interactive()) { win <- 2; boundary <- "none"; meanboundary <- 0 }

source("../sim-context-fns.R")
source("../context-inference-fns.R")

fixfn <- function (...) { 1 }

if (meanboundary==0) {
    genmatrix <- makegenmatrix( patlen=winlen, mutpats=mutpats, selpats=selpats, boundary=boundary, Ne=1e-4 )
} else {
    genmatrix <- meangenmatrix( lwin=meanboundary, rwin=meanboundary, patlen=winlen, mutpats=mutpats, selpats=selpats, boundary=boundary, Ne=1e-4 )
}

save( boundary, meanboundary, bases, mutpats, selpats, fixfn, genmatrix, file=outfile )
