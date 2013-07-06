#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Simulate from the process.\
"

option_list <- list(
        make_option( c("-t","--tlen"), type="numeric", default=NULL, help="Time to simulate for." ),
        make_option( c("-s","--seqlen"), type="numeric", default=NULL, help="Number of bases to simulate." ),
        make_option( c("-b","--interaction"), type="numeric", default=NULL, help="Interaction (beta)."),
        make_option( c("-g","--field"), type="numeric", default=0.0, help="External field (gamma)."),
        make_option( c("-f","--initfreq"), type="numeric", default=.5, help="Initial frequency of X. [default \"%default\"]"),
        make_option( c("-d","--simdir"), type="character", default="ising-sims/", help="Write simulated sequence in this directory [default \"%default\"]" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct logging output to this file. [default appends -simrun.Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (interactive()) { opt$tlen <- .1; opt$seqlen <- 150 }
if ( (is.null(opt$tlen) | is.null(opt$seqlen)) ) { stop("Rscript sim-ising.R -h for help.") }
attach(opt)

if (interactive()) { seqlen <- 1e3; tlen <- .1; interaction <- 1; field <- .5; Xfreq <- 0.5 }

if (!file.exists(simdir)) { dir.create(simdir,recursive=TRUE) }

# identifiers
thisone <- formatC( floor(runif(1)*1e6) , digits=6,flag='0')
now <- Sys.time()

basename <- paste(simdir,"selsims-",format(now,"%Y-%m-%d-%H-%M"),"-",thisone,sep='')
outfile <- paste(basename,".RData",sep='')
if (logfile=="" & !interactive()) { logfile <- paste(basename,"-simrun.Rout",sep='') }
if (logfile != "" & !interactive()) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message") 
    sink(file=logcon, type="output") 
}

bases <- c("X","O")

source("../sim-context-fns.R")
source("../codon-inference-fns.R")

# maximum size of pattern (for simulation)
mutpats <- list( 
    list( c("O","X"), c("X","O") )
    ) 
mutrates <- c(1)
selpats <- list(
        c("OX","XO"),
        c("X")
    )
selcoef <- c(interaction,field)

fixfn <- ising.fixfn

initfreqs <- c(Xfreq,1-Xfreq)
initseq <- rinitseq(seqlen,bases,basefreqs=initfreqs)
system.time( 
        simseqs <- list(
                simseq( seqlen=seqlen, tlen=tlen, mutpats=mutpats, mutrates=mutrates, selpats=selpats, selcoef=selcoef, initseq=initseq, bases=bases, count.trans=FALSE )
            )
    )

save( thisone, now, bases, patlen, mutpats, mutrates, selpats, selcoef, fixfn, seqlen, tlen, initfreqs, simseqs, file=outfile )

