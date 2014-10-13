#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Simulate from the process whose parameters are given in the config file,\
wrapping mutation patterns around to the beginning as needed.\
\
The config file is a JSON file defining: \
    bases : character vector of allowable bases \
    initfreqs : numeric vector of frequencies of the bases at the root \
    mutpats : list of lists of length-two character vectors of mutation patterns \
    mutrates : numeric vector of same length as mutpats giving mutation rates per generation \
    selpats : list of character vectors of patterns under selection (default: empty, in which case don't need the subsequent things either) \
    selcoef : numeric vector of same length as selpats giving selection coefficients \
    fixfn : character string giving name of the fixation function  \
    fixfn.params: named list of additional parameters to pass to fixfn \
\
This extends the information in config files for making a generator matrix. \
\
Example:\
    Rscript ../sim-seq.R -c cpg-model.json -t .01 -s 1000 -o testseq.RData \
"

option_list <- list(
        make_option( c("-c","--configfile"), type="character", help="Configuration file."),
        make_option( c("-t","--tlen"), type="character", help="Time to simulate for." ),
        make_option( c("-s","--seqlen"), type="numeric", help="Number of bases to simulate." ),
        make_option( c("-o","--outfile"), type="character", help="Direct output to this file."),
        make_option( c("-d","--outdir"), type="character", default=".", help="Direct output to this directory with default name if outfile is not specified."),
        make_option( c("-j","--jobid"), type="numeric", default=formatC( floor(runif(1)*1e6) , digits=6,flag='0'), help="Unique identifier to append to output name if not specified."),
        make_option( c("-l","--logfile"), type="character", help="Direct logging output to this file. [default appends -simrun.Rout]" ),
        make_option( c("-z","--seed"), type="integer", help="Seed for pseudorandom number generator; an integer. [default: does not meddle]")
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if ( is.null(opt$configfile) | is.null(opt$seqlen) )  { stop("Rscript sim-seq.R -h for help.") }
if ( !file.exists(opt$configfile) ) { stop("Could not find config file `", opt$configfile, "`.") }

if ( !is.null(opt$seed) ) { set.seed(opt$seed) }

source("../context-inference-fns.R")
source("../sim-context-fns.R")

# identifiers
if (is.null(opt$outfile)) {
    now <- Sys.time()
    opt$outfile <- paste(opt$outdir,"/","simseq-",format(now,"%Y-%m-%d-%H-%M"),"-",opt$jobid,".RData",sep='')
}

if (!is.null(opt$outdir)) {
    if (!file.exists(opt$outdir)) { dir.create(opt$outdir) }
    if (dirname(opt$outfile) != opt$outdir) { opt$outfile <- paste(opt$outdir,opt$outfile,sep='/') }
}
basename <- gsub(".RData","",opt$outfile)
if (is.null(opt$logfile) & !interactive()) { opt$logfile <- paste(basename,".Rout",sep='') }
if (!is.null(opt$logfile)) {
    logcon <- if (opt$logfile=="-") { stdout() } else { file(opt$logfile,open="wt") }
    sink(file=logcon, type="message")
    sink(file=logcon, type="output")
}

config <- read.config(opt$configfile)
# treeify, add defaults, etc
config <- treeify.config(config,tlen=opt$tlen)
config <- parse.models(config)
# error checks
stopifnot( ( length(config$bases) == length(config$initfreqs) ) )

# return a list of the simulated sequences in the same order as the tips,nodes of the tree
simseqs <- simseq.tree(opt$seqlen,config)

simseq.opt <- opt
simseq.config <- config

save( simseq.opt, simseq.config, simseqs, file=opt$outfile )


