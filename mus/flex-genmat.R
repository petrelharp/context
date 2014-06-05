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
'

option_list <- list(
        make_option( c("-c","--configfile"), type="character", default="-", help="Input file of mutation and selection motifs." ),
        make_option( c("-s","--outfile"), type="character", default="", help="Save resulting matrix in this file.  [default: genmatrices/genmatrix-(configfile).RData]" ),
        make_option( c("-w","--winlen"), type="integer", default=2, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-b","--boundary"), type="character", default="none", help="Boundary conditions. [default \"%default\"]" ),
        make_option( c("-m","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases. [default \"%default\"]" ),
        make_option( c("-k","--patlen"), type="integer", default=1, help="Include mutation rates for all tuples of this length. [default \"%default\"]" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct logging output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$configfile)) { stop("Need a config file.") }
options(error=traceback)

config <- fromJSON(configfile)
attach(c(opt,config))

if (outfile=="") { 
    if (!file.exists("genmatrices")) { dir.create("genmatrices") }; 
    outfile <- paste(paste("genmatrices/genmatrix",gsub(".mut$","",configfile),sep="-"),".RData",sep='') 
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
attr(genmatrix,"bases") <- bases
attr(genmatrix,"mutpats") <- mutpats
attr(genmatrix,"selpats") <- selpats
attr(genmatrix,"boundary") <- boundary
attr(genmatrix,"meanboundary") <- meanboundary

genmat.opt <- opt
save( genmat.opt, winlen, boundary, meanboundary, bases, mutpats, selpats, fixfn, genmatrix, file=outfile )
