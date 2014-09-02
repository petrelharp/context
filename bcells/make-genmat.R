#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Pre-compute a generator matrix and associated structures, including terms for all k-tuples.\
"

option_list <- list(
        make_option( c("-s","--outfile"), type="character", default="", help="Save resulting matrix in this file.  [default: genmatrix-longwin-boundary-meanboundary.RData]" ),
        make_option( c("-w","--longwin"), type="integer", default=2, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-b","--boundary"), type="character", default="none", help="Boundary conditions. [default \"%default\"]" ),
        make_option( c("-m","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases. [default \"%default\"]" ),
        make_option( c("-k","--patlen"), type="integer", default=1, help="Include mutation rates for all tuples of this length. [default \"%default\"]" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct logging output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)
options(error=traceback)

if (outfile=="") { 
    if (!file.exists("genmatrices")) { dir.create("genmatrices") }; 
    outfile <- paste(paste("genmatrices/genmatrix",longwin,boundary,meanboundary,patlen,sep="-"),".RData",sep='') 
}
basename <- gsub(".RData",'',outfile)
if (logfile=="" & !interactive()) { logfile <- paste(basename,"-make-genmat.Rout",sep='') }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message") 
    sink(file=logcon, type="output") 
}

if (interactive()) { shortwin <- 2; boundary <- "none"; meanboundary <- 0 }

bases <- c("A","T","C","G")

source("../sim-context-fns.R")
source("../context-inference-fns.R")

mutpats <- getmutpats(patlen)
# # DO NOT allow GC-bias
# selpats <- list(
#         c("A","T")
#     )
# fixfn <- popgen.fixfn  # takes (dx,Ne) as arguments
selpats <- list( )
fixfn <- function (...) { 1 }
if (meanboundary==0) {
    genmatrix <- makegenmatrix( patlen=longwin, mutpats=mutpats, selpats=selpats, boundary=boundary, Ne=1e-4 )
} else {
    genmatrix <- meangenmatrix( leftwin=meanboundary, rightwin=meanboundary, patlen=longwin, mutpats=mutpats, selpats=selpats, boundary=boundary, Ne=1e-4 )
}

genmat.opt <- opt
save( genmat.opt, longwin, boundary, meanboundary, bases, mutpats, selpats, fixfn, genmatrix, file=outfile )
