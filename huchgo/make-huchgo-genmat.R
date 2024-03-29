#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Pre-compute a generator matrix for the CpG + GC biased process.\
"

option_list <- list(
        make_option( c("-s","--outfile"), type="character", default="", help="Save resulting matrix in this file.  [default: genmatrix-longwin-boundary-meanboundary.RData]" ),
        make_option( c("-w","--longwin"), type="integer", default=2, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-b","--boundary"), type="character", default="none", help="Boundary conditions. [default \"%default\"]" ),
        make_option( c("-m","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases. [default \"%default\"]" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct logging output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)
options(error=traceback)

if (outfile=="") { 
    if (!file.exists("genmatrices")) { dir.create("genmatrices") }; 
    outfile <- paste(paste("genmatrices/genmatrix",longwin,boundary,meanboundary,sep="-"),".RData",sep='') 
}
basename <- gsub(".RData",'',outfile)
if (logfile=="" & !interactive()) { logfile <- paste(basename,"-simrun.Rout",sep='') }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message") 
    sink(file=logcon, type="output") 
}

if (interactive()) { shortwin <- 2; boundary <- "none"; meanboundary <- 0 }

bases <- c("A","T","C","G")


library(contextual)
library(contextutils)
library(simcontext)

# maximum size of pattern (for simulation)
mutpats <- c(
        apply(combn(bases,2),2,list),  # single-base rates
        apply(combn(bases,2)[2:1,],2,list),  # single-base rates
        list( list( c("CG","TG"), c("CG","CA") ) )  # CpG rate
    ) 
selpats <- list(
        c("A","T")
    )
fixfn <- popgen.fixfn  # takes (dx,Ne) as arguments
if (meanboundary==0) {
    genmatrix <- makegenmatrix( patlen=longwin, mutpats=mutpats, selpats=selpats, boundary=boundary, Ne=1e-4 )
} else {
    genmatrix <- meangenmatrix( leftwin=meanboundary, rightwin=meanboundary, patlen=longwin, mutpats=mutpats, selpats=selpats, boundary=boundary, Ne=1e-4 )
}

save( opt, longwin, boundary, meanboundary, bases, mutpats, selpats, fixfn, genmatrix, file=outfile )
