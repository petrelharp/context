#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Infer parameters from output of sim-tasep.R .\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", default=NULL, help=".RData file containing simulation." ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default parse-mcmc.log]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)
options(error=traceback)

if (is.null(infile) | is.null(nbatches)) { cat("Run\n parse-mcmc.R -h\n for help.") }

if (logfile!="" & !interactive()) { 
    basedir <- gsub(".RData","/",infile,fixed=TRUE)
    logcon <- if (logfile=="-") { stdout() } else { file(paste(basedir,logfile,sep=''),open="wt") }
    sink(file=logcon, type="message") 
    sink(file=logcon, type="output")   # send both to log file
}



