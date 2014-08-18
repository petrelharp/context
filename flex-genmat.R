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
    ], \
    "fixfn" : "null.fixfn"
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
        make_option( c("-c","--configfile"), type="character", help="Input file of mutation and selection motifs." ),
        make_option( c("-s","--outfile"), type="character", default="", help="Save resulting matrix in this file.  [default: genmatrix-(longwin)-(configfile).RData]" ),
        make_option( c("-w","--longwin"), type="integer", help="Size of matching window." ),
        make_option( c("-b","--boundary"), type="character", default="none", help="Boundary conditions. [default \"%default\"]" ),
        make_option( c("-m","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases. [default \"%default\"]" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct logging output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$configfile) | !file.exists(opt$configfile)) { stop("Need a config file.") }
if (is.null(opt$longwin)) { stop("Need a window length.") }
options(error=traceback)

config <- fromJSON(opt$configfile,simplifyMatrix=FALSE)
attach(c(opt,config))

if (outfile=="") { 
    outfile <- paste(dirname(configfile),"/",paste("genmatrix",longwin,gsub(".json$","",basename(configfile)),sep="-"),".RData",sep='') 
}
basename <- gsub(".RData",'',outfile)
if (logfile=="" & !interactive()) { logfile <- paste(basename,"-flex-genmat.Rout",sep='') }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message") 
    sink(file=logcon, type="output") 
}

if (interactive()) { shortwin <- 2; boundary <- "none"; meanboundary <- 0 }

source("../sim-context-fns.R")
source("../context-inference-fns.R")

# turn fixfn into an actual function
# either by looking it up as a name
# or parsing it directly
if (!exists("selpats")) { selpats <- list(); selcoef <- numeric(0); fixfn <- null.fixfn; fixfn.params=list() }
if (!exists("fixfn")) { fixfn <- null.fixfn; fixfn.params=list() }
if (is.character(fixfn)) {
    if (exists(fixfn,mode="function")) {
        fixfn <- get(fixfn,mode="function")
    } else {
        fixfn <- eval(parse(text=fixfn))
    }
}

if (meanboundary==0) {
    genmatrix <- makegenmatrix( patlen=opt$longwin, mutpats=mutpats, selpats=selpats, boundary=boundary, bases=bases, fixfn=fixfn )
} else {
    genmatrix <- meangenmatrix( leftwin=meanboundary, rightwin=meanboundary, patlen=opt$longwin, mutpats=mutpats, selpats=selpats, boundary=boundary, bases=bases, fixfn=fixfn )
}

save( boundary, meanboundary, bases, mutpats, selpats, fixfn, genmatrix, file=outfile )
