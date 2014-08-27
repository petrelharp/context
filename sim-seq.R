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
        make_option( c("-l","--logfile"), type="character", help="Direct logging output to this file. [default appends -simrun.Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if ( is.null(opt$configfile) | is.null(opt$seqlen)) ) { stop("Rscript sim-seq.R -h for help.") }
if ( !file.exists(opt$configfile) ) { stop("Could not find config file `", opt$configfile, "`.") }

source("../context-inference-fns.R")
source("../sim-context-fns.R")

config <- read.config(opt$configfile)

# parse tree, and edge-associated models
if (is.null(config$tree)) {
    if (is.null(opt$tlen)) { stop("Must specify timelengths on the tree.") }
    config <- list( tree="(tip);", bases=config$bases, tip=config, initfreqs=config$initfreqs )
}
config$tree <- ape::read.tree(text=config$tree)
if (is.null(config$tree$edge.length)) { 
    # if tree has no edge lengths, bring these in from tlen
    config$tree$edge.length <- eval(parse(text=opt$tlen))
} else if (!is.null(opt$tlen)) { 
    warning("Branch lengths specified in config file and on command line; ignoring the command line.") 
}

if (is.null(config$tree$tip.label) | is.null(config$tree$node.label)) { stop("Please label tips and nodes on the tree.") }

# edges are labeled by the node/tip below them:
rootname <- c( config$tree$tip.label, config$tree$node.label )[ get.root(config$tree) ]
edges <- selfname( setdiff( c( config$tree$tip.label, config$tree$node.label ), rootname ) )
if (!all(edges %in% names(config)) ) { stop("Must specify named models for each edge in the tree.") }

stopifnot( ( length(config$bases) == length(config$initfreqs) ) )

for (edge in edges) {
    if (!is.character(config[[edge]])) {
        # turn fixfn into a function and check we have the right parameters
        config[[edge]]$fixfn <- parse.fixfn(config$fixfn,config[[edge]]$fixfn.params)
        # put in defaults if no selection
        if (is.null(config[[edge]]$selpats)) { config[[edge]]$selpats <- list(); config[[edge]]$selcoef <- numeric(0); config[[edge]]$fixfn <- null.fixfn; config[[edge]]$fixfn.params=list() }
        stopifnot( with( config[[edge]], ( length(mutpats) == length(mutrates) ) && ( length(selpats) == length(selcoef) )))
    }
}

# config$fixfn <- parse.fixfn(config$fixfn,config$fixfn.params)
# # defaults:
# if (is.null(config$selpats)) { config$selpats <- list(); config$selcoef <- numeric(0); config$fixfn <- null.fixfn; config$fixfn.params=list() }
# stopifnot( with( config, ( length(bases) == length(initfreqs) ) && ( length(mutpats) == length(mutrates) ) && ( length(selpats) == length(selcoef) )))


# identifiers
if (is.null(opt$outfile)) {
    now <- Sys.time()
    basename <- paste(opt$outdir,"/","simseq-",format(now,"%Y-%m-%d-%H-%M"),"-",opt$jobid,sep='')
    opt$outfile <- paste(basename,".RData",sep='')
} else {
    basename <- gsub(".RData","",opt$outfile)
}
if (is.null(opt$logfile) & !interactive()) { opt$logfile <- paste(basename,".Rout",sep='') }
if (!is.null(opt$logfile)) {
    logcon <- if (opt$logfile=="-") { stdout() } else { file(opt$logfile,open="wt") }
    sink(file=logcon, type="message")
    sink(file=logcon, type="output")
}


initseq <- rinitseq(opt$seqlen,config$bases,basefreqs=config$initfreqs)
system.time(
        simseqs <- list(
                do.call( simseq, c( list( seqlen=opt$seqlen, tlen=opt$tlen, mutpats=config$mutpats, mutrates=config$mutrates, selpats=config$selpats, selcoef=config$selcoef, initseq=initseq, bases=config$bases, fixfn=config$fixfn ), config$fixfn.params ) )
            )
    )

simseq.opt <- opt

save( simseq.opt, config, simseqs, file=opt$outfile )


