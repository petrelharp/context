#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Pre-compute a generator matrix for the process with codon -> codon transitions.\
"

option_list <- list(
        make_option( c("-s","--outfile"), type="character", default="", help="Save resulting matrix in this file.  [default: genmatrix-winlen-boundary-meanboundary.RData]" ),
        make_option( c("-w","--winlen"), type="integer", default=3, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-b","--boundary"), type="character", default="none", help="Boundary conditions. [default \"%default\"]" ),
        make_option( c("-m","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases. [default \"%default\"]" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct logging output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)
options(error=traceback)

if (outfile=="") { 
    if (!file.exists("genmatrices")) { dir.create("genmatrices") }; 
    outfile <- paste(paste("genmatrices/genmatrix",winlen,boundary,meanboundary,sep="-"),".RData",sep='') 
}
basename <- gsub(".RData",'',outfile)
if (logfile=="" & !interactive()) { logfile <- paste(basename,"-simrun.Rout",sep='') }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message") 
    sink(file=logcon, type="output") 
}

if (interactive()) { winlen <- 3; boundary <- "none"; meanboundary <- 0 }

bases <- c("A","T","C","G")

source("../sim-context-fns.R")
source("../codon-inference-fns.R")
source("../codons.R")
codonstrings <- as.character(codons$codon)

# all pairwise trinucleotide
mutpats <- c(
        apply(combn(codonstrings,2),2,list),  # single-base rates
        apply(combn(codonstrings,2)[2:1,],2,list),  # single-base rates
        list( list( c("CG","TG"), c("CG","CA") ) )  # CpG rate
    ) 
# restrict to those with only one transition
nchanges <- sapply( lapply( lapply( changepos(mutpats), lapply, length ), unlist ), max )  # there would be a less elegant way to do this
mutpats <- mutpats[ nchanges==1 ]

selpats <- as.list(codonstrings)
fixfn <- function (...) { 1 }

if (meanboundary>0) {
    genmatrix <- meangenmatrix( lwin=meanboundary, rwin=meanboundary, patlen=winlen, mutpats=mutpats, selpats=selpats, boundary=boundary )
} else {
    genmatrix <- makegenmatrix( patlen=winlen, mutpats=mutpats, selpats=selpats, boundary=boundary )
}

save( opt, winlen, boundary, meanboundary, bases, mutpats, selpats, fixfn, genmatrix, file=outfile )
