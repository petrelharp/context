#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Pre-compute a generator matrix and associated structures, including terms for all k-tuples.\
"

option_list <- list(
        make_option( c("-s","--outfile"), type="character", default="", help="Save resulting matrix in this file.  [default: genmatrix-winlen-boundary-meanboundary.RData]" ),
        make_option( c("-w","--winlen"), type="integer", default=2, help="Size of matching window. [default \"%default\"]" ),
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
    outfile <- paste(paste("genmatrices/genmatrix",winlen,boundary,meanboundary,patlen,sep="-"),".RData",sep='') 
}
basename <- gsub(".RData",'',outfile)
if (logfile=="" & !interactive()) { logfile <- paste(basename,"-make-genmat.Rout",sep='') }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message") 
    sink(file=logcon, type="output") 
}

if (interactive()) { win <- 2; boundary <- "none"; meanboundary <- 0 }

bases <- c("A","T","C","G")

source("../sim-context-fns.R")
source("../codon-inference-fns.R")

mutpats <- list()
# all kmer -> kmer changes
for (k in 1:patlen) {
    kmers <- getpatterns(k)
    mutpats <- c( mutpats,
            apply(combn(kmers,2),2,list),  
            apply(combn(kmers,2)[2:1,],2,list)
        )
}
# restrict to those only changing one base
nchanges <- sapply(mutpatchanges(mutpats),nrow)
mutpats <- mutpats[nchanges==1]
# # DO NOT allow GC-bias
# selpats <- list(
#         c("A","T")
#     )
# fixfn <- popgen.fixfn  # takes (dx,Ne) as arguments
selpats <- list( )
fixfn <- function (...) { 1 }
if (meanboundary==0) {
    genmatrix <- makegenmatrix( patlen=winlen, mutpats=mutpats, selpats=selpats, boundary=boundary, Ne=1e-4 )
} else {
    genmatrix <- meangenmatrix( lwin=meanboundary, rwin=meanboundary, patlen=winlen, mutpats=mutpats, selpats=selpats, boundary=boundary, Ne=1e-4 )
}

genmat.opt <- opt
save( genmat.opt, winlen, boundary, meanboundary, bases, mutpats, selpats, fixfn, genmatrix, file=outfile )
