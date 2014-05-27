#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Simulate from the process.\
\
"

usage <- "\
Simulate from the process.\
"

option_list <- list(
        make_option( c("-s","--seqlen"), type="numeric", default=NULL, help="Number of bases to simulate." ),
        make_option( c("-t","--tlen"), type="numeric", default=NULL, help="Time to simulate for." ),
        make_option( c("-f","--Xfreq"), type="numeric", default=0.5, help="Initial frequency of 'X'. [default \"%default\"]"),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (interactive()) { opt$tlen <- 1; opt$seqlen <- 150 }
if ( (is.null(opt$tlen) | is.null(opt$seqlen)) ) { stop("Rscript sim-tasep -h for help.") }
attach(opt)


# identifiers
thisone <- formatC( floor(runif(1)*1e6) , digits=6,flag='0')
now <- Sys.time()

simdir <- "tasep-sims/"
if (!file.exists(simdir)) { dir.create(simdir,recursive=TRUE) }

bases <- c("X","O")

source("../sim-context-fns.R")
source("../context-inference-fns.R")

# maximum size of pattern (for simulation)
mutpats <- list( 
    list( c("XO","OX") )
    ) 
mutrates <- 1
selpats <- list()
selcoef <- numeric(0)

fixfn <- function (...) { 1 }

initfreqs <- c(Xfreq,1-Xfreq)
initseq <- rinitseq(seqlen,bases,basefreqs=initfreqs)
system.time( 
        simseqs <- list(
                simseq( seqlen=seqlen, tlen=tlen, mutpats=mutpats, mutrates=mutrates, selpats=selpats, selcoef=selcoef, initseq=initseq, bases=bases )
            )
    )

save( thisone, now, bases, patlen, mutpats, mutrates, selpats, selcoef, fixfn, seqlen, tlen, initfreqs, simseqs, file=paste(simdir,"selsims-",format(now,"%Y-%m-%d-%H-%M"),"-",thisone,".RData",sep='') )


