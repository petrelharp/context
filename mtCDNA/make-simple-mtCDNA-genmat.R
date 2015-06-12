#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Pre-compute a generator matrix for the process with codon -> codon transitions.\
"

option_list <- list(
        make_option( c("-s","--outfile"), type="character", default="", help="Save resulting matrix in this file.  [default: simple-genmatrices/genmatrix-longwin-boundary-meanboundary.RData]" ),
        make_option( c("-w","--longwin"), type="integer", default=3, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-b","--boundary"), type="character", default="none", help="Boundary conditions. [default \"%default\"]" ),
        make_option( c("-m","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases. [default \"%default\"]" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct logging output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)
options(error=traceback)

if (outfile=="") { 
    if (!file.exists("simple-genmatrices")) { dir.create("simple-genmatrices") }; 
    outfile <- paste(paste("simple-genmatrices/genmatrix",longwin,boundary,meanboundary,sep="-"),".RData",sep='') 
}
basename <- gsub(".RData",'',outfile)
if (logfile=="" & !interactive()) { logfile <- paste(basename,"-simrun.Rout",sep='') }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message") 
    sink(file=logcon, type="output") 
}

if (interactive()) { longwin <- 3; boundary <- "none"; meanboundary <- 0 }

bases <- c("A","T","C","G")

source("../sim-context-fns.R",chdir=TRUE)
source("../context-inference-fns.R",chdir=TRUE)
source("../context.R",chdir=TRUE)
codonstrings <- as.character(codons$codon)

# all pairwise bases
mutpats <- c(
        apply(combn(bases,2),2,list),  # single-base rates
        apply(combn(bases,2)[2:1,],2,list),  # single-base rates
        list( list( c("CG","TG"), c("CG","CA") ) ),  # CpG rate
        list( list( c("GC","AA"), c("GC","TT") ), list( c("GA","TT"), c("TC","AA") ) ) # Kelley's dinucleotides
    ) 

selpats <- list()
fixfn <- function (...) { 1 }

if (meanboundary>0) {
    genmatrix <- meangenmatrix( leftwin=meanboundary, rightwin=meanboundary, patlen=longwin, mutpats=mutpats, selpats=selpats, boundary=boundary )
} else {
    genmatrix <- makegenmatrix( patlen=longwin, mutpats=mutpats, selpats=selpats, boundary=boundary )
}

save( opt, longwin, boundary, meanboundary, bases, mutpats, selpats, fixfn, genmatrix, file=outfile )
