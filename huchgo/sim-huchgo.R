#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Simulate from the process.\
"

option_list <- list(
        make_option( c("-t","--tlen"), type="character", default=NULL, help="Branch length(s), in the order (Root->Go), (Root->HuCh), (HuCh->Ch), (HuCh->Hu).  If one is provided, this is multiplied by (1,.25,.75,.75)."),
        make_option( c("-s","--seqlen"), type="numeric", default=NULL, help="Number of bases to simulate." ),
        make_option( c("-m","--baserates"), type="character", default=NULL, help="Single base transition rate (or list of twelve rates for A->T A->C A->G T->C T->G C->G T->A C->A G->A C->T G->T G->C). [default indep't Uniform[0,1]]" ),
        make_option( c("-c","--cpgrate"), type="numeric", default=NULL, help="Additional CpG rate. [default 2*Uniform[0,1]]"),
        make_option( c("-g","--gcbias"), type="numeric", default=NULL, help="Strength of GC bias. [default 2e-4*Uniform[0,1]]"),
        make_option( c("-n","--Ne"), type="character", default="1e4", help="Effective population sizes, in same order as branch lengths (for computing prob of fixation). [default \"%default\"]"),
        make_option( c("-f","--initfreqs"), type="character", default="c(.25,.25,.25,.25)", help="Initial base frequencies. [default \"%default\"]"),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends -simrun.Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (interactive()) { opt$tlen <- ".1"; opt$seqlen <- 150 }
if ( (is.null(opt$tlen) | is.null(opt$seqlen)) ) { stop("Rscript sim-huchgo.R -h for help.") }
if (is.null(opt$baserates)) { opt$baserates <- runif(12); names(opt$baserates) <- c("A->T", "A->C", "A->G", "T->C", "T->G", "C->G", "T->A", "C->A", "G->A", "C->T", "G->T", "G->C") }
if (is.character(opt$baserates)) { opt$baserates <- eval(parse(text=opt$baserates)) }
if (length(opt$baserates)==1) { opt$baserates <- rep(opt$baserates,12) }
if (is.null(opt$cpgrate)) { opt$cpgrate <- 2*runif(1) }
if (is.null(opt$gcbias)) { opt$gcbias <- 2e-4*runif(1) }
opt$initfreqs <- eval(parse(text=opt$initfreqs))
opt$tlen <- eval(parse(text=opt$tlen))
if (length(opt$tlen)==1) { opt$tlen <- opt$tlen * c(1,.25,.75,.75) }
opt$Ne <- eval(parse(text=opt$Ne))
if (length(opt$Ne)==1) { opt$Ne <- rep( opt$Ne, 4 ) }
attach(opt)

simdir <- "huchgo-sims/"
if (!file.exists(simdir)) { dir.create(simdir,recursive=TRUE) }

# identifiers
thisone <- formatC( floor(runif(1)*1e6) , digits=6,flag='0')
now <- Sys.time()

basename <- paste(simdir,"selsims-",format(now,"%Y-%m-%d-%H-%M"),"-",thisone,sep='')
outfile <- paste(basename,".RData",sep='')
if (logfile=="" & !interactive()) { logfile <- paste(basename,"-simrun.Rout",sep='') }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message") 
    sink(file=logcon, type="output") 
}

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
mutrates <- c(baserates,cpgrate)
selpats <- list(
        c("A","T")
    )
selcoef <- c( gcbias )

if (is.null(names(initfreqs))) { names(initfreqs) <- bases }

fixfn <- popgen.fixfn  # takes (dx,Ne) as arguments

rootseq <- rinitseq(seqlen,bases,basefreqs=initfreqs)
huchseq <- simseq( seqlen=seqlen, tlen=tlen[2], mutpats=mutpats, mutrates=mutrates, selpats=selpats, selcoef=selcoef, initseq=rootseq, bases=bases, Ne=Ne[2] )
system.time( 
        simseqs <- list(
                Go=simseq( seqlen=seqlen, tlen=tlen[1], mutpats=mutpats, mutrates=mutrates, selpats=selpats, selcoef=selcoef, initseq=rootseq, bases=bases, count.trans=TRUE, Ne=Ne[1] ),
                Ch=simseq( seqlen=seqlen, tlen=tlen[3], mutpats=mutpats, mutrates=mutrates, selpats=selpats, selcoef=selcoef, initseq=huchseq$finalseq, bases=bases, count.trans=TRUE, Ne=Ne[3] ),
                Hu=simseq( seqlen=seqlen, tlen=tlen[4], mutpats=mutpats, mutrates=mutrates, selpats=selpats, selcoef=selcoef, initseq=huchseq$finalseq, bases=bases, count.trans=TRUE, Ne=Ne[4] )
            )
    )

save( thisone, now, opt, bases, mutpats, mutrates, selpats, selcoef, fixfn, seqlen, tlen, Ne, initfreqs, simseqs, file=outfile )


