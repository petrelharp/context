#!/usr/bin/Rscript --vanilla
require(optparse)
setwd('/Users/Jessica/Documents/USC/context/sims-selection')
usage <- "\
Simulate from the process.\
"

option_list <- list(
        make_option( c("-t","--tlen"), type="character", default=NULL, help="Branch length(s). [if one provided, both are equal]"),
        make_option( c("-s","--seqlen"), type="numeric", default=NULL, help="Number of bases to simulate." ),
        make_option( c("-m","--baserates"), type="character", default=NULL, help="Single base transition rate (or list of twelve rates for A->T A->C A->G T->C T->G C->G T->A C->A G->A C->T G->T G->C). [default indep't Uniform[0,1]]" ),
        make_option( c("-c","--cpgrate"), type="character", default=NULL, help="Additional CpG rate. [default 2*Uniform[0,1]]"),
        make_option( c("-g","--gcbias"), type="character", default=NULL, help="Strength of GC bias. [default 2e-4*Uniform[0,1]]"),
        make_option( c("-n","--Ne"), type="character", default="1e4", help="Effective population sizes, in same order as branch lengths (for computing prob of fixation). [default \"%default\"]"),
	make_option( c("-f","--initfreqs"), type="character", default="c(.25,.25,.25,.25)", help="Initial base frequencies. [default \"%default\"]"),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends -simrun.Rout]" )
    )

# added options for, Ne, gc bias, dual mutation rates (one set for each branch)
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if ( (is.null(opt$tlen) | is.null(opt$seqlen)) ) { stop("Rscript sim-tree-cpg.R -h for help.") }
if (is.null(opt$baserates)) { opt$baserates <- c(runif(12), runif(12)); names(opt$baserates) <- rep(c("A->T", "A->C", "A->G", "T->C", "T->G", "C->G", "T->A", "C->A", "G->A", "C->T", "G->T", "G->C"), 2) }
if (is.character(opt$baserates)) { opt$baserates <- eval(parse(text=opt$baserates)) }
if (length(opt$baserates)==1) { opt$baserates <- rep(opt$baserates,12) }
if (is.null(opt$cpgrate)) { opt$cpgrate <- 2*runif(2) }
if (is.null(opt$gcbias)) { opt$gcbias <- 2e-4*runif(2) }
opt$Ne <- eval(parse(text=opt$Ne))
if (length(opt$Ne)==1) { opt$Ne <- rep( opt$Ne, 2 ) }
opt$initfreqs <- eval(parse(text=opt$initfreqs))
opt$tlen <- eval(parse(text=opt$tlen))
if (length(opt$tlen)==1) { opt$tlen <- rep(opt$tlen,2) }
attach(opt)

# identifiers
thisone <- formatC( floor(runif(1)*1e6) , digits=6,flag='0')
now <- Sys.time()

basename <- paste("selsims_seqlen_", opt$seqlen, "-simrun" ,sep='')
outfile <- paste(basename,".RData",sep='')
if (logfile=="" & !interactive()) { logfile <- paste(basename,".Rout",sep='') }
if (!is.null(logfile)) { 
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message") 
    sink(file=logcon, type="output") 
}

bases <- c("A","T","C","G")

source("../sim-context-fns-dual.R",chdir=TRUE)
source("../context-inference-fns.R",chdir=TRUE)

# maximum size of pattern (for simulation)
mutpats <- c(
        apply(combn(bases,2),2,list),  # single-base rates
        apply(combn(bases,2)[2:1,],2,list),  # single-base rates
        list( list( c("CG","TG"), c("CG","CA") ) )  # CpG rate
    ) 
mutrates <- list(c(baserates[1:12],cpgrate[1]), c(baserates[13:24], cpgrate[2]))
selpats <- list(
	c("A", "T")
    )
selcoef <- c(gcbias)

if (is.null(names(initfreqs))) { names(initfreqs) <- bases }

fixfn <- popgen.fixfn  # takes (dx,Ne) as arguments

rootseq <- rinitseq(seqlen,bases,basefreqs=initfreqs)
system.time( 
        simseqs <- list(
                simseq( seqlen=seqlen, tlen=tlen[1], mutpats=mutpats, mutrates=mutrates[[1]], selpats=selpats, selcoef=selcoef[1], initseq=rootseq, bases=bases, Ne=Ne[1] ),
                simseq( seqlen=seqlen, tlen=tlen[2], mutpats=mutpats, mutrates=mutrates[[2]], selpats=selpats, selcoef=selcoef[2], initseq=rootseq, bases=bases, Ne=Ne[2] )
            )
    )

save( thisone, now, opt, bases, mutpats, mutrates, selpats, selcoef, fixfn, seqlen, tlen, Ne, initfreqs, simseqs, file=outfile )
