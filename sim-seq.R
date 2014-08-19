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
    fixparams: named list of additional parameters to pass to fixfn \
\
This extends the information in config files for making a generator matrix. \
"

option_list <- list(
        make_option( c("-c","--configfile"), type="character", help="Configuration file."),
        make_option( c("-t","--tlen"), type="numeric", help="Time to simulate for." ),
        make_option( c("-s","--seqlen"), type="numeric", help="Number of bases to simulate." ),
        make_option( c("-o","--outfile"), type="character", default="", help="Direct output to this file."),
        make_option( c("-d","--outdir"), type="character", default="", help="Direct output to this directory with default name if outfile is not specified."),
        make_option( c("-j","--jobid"), type="numeric", default=formatC( floor(runif(1)*1e6) , digits=6,flag='0'), help="Unique identifier to append to output name if not specified."),
        make_option( c("-l","--logfile"), type="character", default="", help="Direct logging output to this file. [default appends -simrun.Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if ( is.null(opt$configfile) | (is.null(opt$tlen) | is.null(opt$seqlen)) ) { stop("Rscript sim-seq.R -h for help.") }
if ( !file.exists(opt$configfile) ) { stop("Could not find config file `", opt$configfile, "`.") }

source("../context-inference-fns.R")
source("../sim-context-fns.R")

require(jsonlite)
config <- fromJSON(opt$configfile,simplifyMatrix=FALSE)
attach(config)

# turn fixfn into an actual function
# either by looking it up as a name
# or parsing it directly
if (exists("fixfn") && is.character(fixfn)) {
    if (exists(fixfn,mode="function")) {
        fixfn <- get(fixfn,mode="function")
    } else {
        fixfn <- eval(parse(text=fixfn))
    }
}

# defaults:
if (!exists("selpats")) { selpats <- list(); selcoef <- numeric(0); fixfn <- null.fixfn; fixfn.params=list() }
stopifnot(
          ( length(bases) == length(initfreqs) ) &&
           ( length(mutpats) == length(mutrates) ) &&
           ( length(selpats) == length(selcoef) )
        )


# identifiers
if (opt$outfile == "") {
    now <- Sys.time()
    basename <- paste(outdir,"/","simseq-",format(now,"%Y-%m-%d-%H-%M"),"-",opt$jobid,sep='')
    opt$outfile <- paste(basename,".RData",sep='')
} else {
    basename <- gsub(".RData","",opt$outfile)
}
if (opt$logfile=="" & !interactive()) { opt$logfile <- paste(basename,".Rout",sep='') }
if (!is.null(opt$logfile)) {
    logcon <- if (opt$logfile=="-") { stdout() } else { file(opt$logfile,open="wt") }
    sink(file=logcon, type="message")
    sink(file=logcon, type="output")
}


initseq <- rinitseq(opt$seqlen,bases,basefreqs=initfreqs)
system.time(
        simseqs <- list(
                simseq( seqlen=opt$seqlen, tlen=opt$tlen, mutpats=mutpats, mutrates=mutrates, selpats=selpats, selcoef=selcoef, initseq=initseq, bases=bases, fixfn=fixfn )
            )
    )

simseq.opt <- opt

save( simseq.opt, bases, mutpats, mutrates, selpats, selcoef, fixfn, initfreqs, simseqs, file=opt$outfile )


