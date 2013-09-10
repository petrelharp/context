#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Infer parameters from output of sim-tree-cpg.R .\
"

option_list <- list(
        make_option( c("-c","--infile"), type="character", default="mtCDNApri-sub.nuc", help=".RData file containing simulation." ),
        make_option( c("-w","--win"), type="integer", default=1, help="Size of matching window. [default \"%default\"]" ),
        make_option( c("-l","--lwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-r","--rwin"), type="integer", default=2, help="Size of left-hand context. [default \"%default\"]" ),
        make_option( c("-n","--nbatches"), type="integer", default=20, help="Number of MCMC batches. [default \"%default\"]" ),
        make_option( c("-b","--blen"), type="integer", default=10, help="Length of each MCMC batch. [default \"%default\"]" ),
        make_option( c("-s","--stepscale"), type="numeric", default=1e-4, help="Scale of proposal steps for Metropolis algorithm. [default \"%default\"]" ),
        make_option( c("-m","--mmean"), type="double", default=1, help="Prior mean on single base mutation rates. [default \"%default\"]" ),
        make_option( c("-c","--cpgmean"), type="double", default=1, help="Prior variance on CpG rate. [default \"%default\"]" ),
        make_option( c("-p","--pprior"), type="double", default=1, help="Parameter for Dirichlet prior on base frequencies. [default \"%default\"]" ),
        make_option( c("-v","--tprior"), type="double", default=.5, help="Parameter for Beta prior on branch length. [default \"%default\"]" ),
        make_option( c("-d","--boundary"), type="character", default="none", help="Boundary conditions for generator matrix. [default \"%default\"]"),
        make_option( c("-y","--meanboundary"), type="integer", default=0, help="Average over this many neighboring bases in computing generator matrix. [default \"%default\"]" ),
        make_option( c("-g","--gmfile"), type="character", default="TRUE", help="File with precomputed generator matrix, or TRUE [default] to look for one. (otherwise, will compute)"),
        make_option( c("-e","--countfile"), type="character", default="", help="Record tuple counts in this file. [default: countdata/counts-(lwin)-(win)-(rwin).RData]" ),
        make_option( c("-o","--logfile"), type="character", default="", help="Direct output to this file. [default appends .Rout]" )
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
attach(opt)
options(error=traceback)

winlen <- lwin+win+rwin

if (gmfile=="TRUE") { gmfile <- paste(paste("genmatrices/genmatrix",winlen,boundary,meanboundary,sep="-"),".RData",sep='') }

if (is.null(infile)) { cat("Run\n  cpg-tree-inference.R -h\n for help.") }

if (logfile!="") {
    logfile <- gsub(".RData",".Rout",infile,fixed=TRUE)
    logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
    sink(file=logcon, type="message", split=interactive()) 
    sink(file=logcon, type="output", split=interactive())   # send both to log file
}

if (countfile=="") { countfile <- paste(paste("countdata/counts",lwin,win,rwin,sep="-"),".RData",sep="") }


# took mtCDNApri.nuc from paml/examples/mtCDNA
#  and: removed extra spaces on first line (only one space between numbers)
#       made each sequence be on a single line
#       removed trailing text

require(Biostrings)

scriptdir <- "../"
source(paste(scriptdir,"codon-inference-fns.R",sep=''))
source(paste(scriptdir,"sim-context-fns.R",sep=''))

mtCDNA <- as( readDNAMultipleAlignment(infile, format="phylip" ), "DNAStringSet" )

basedir <- gsub(".nuc","",infile,fixed=TRUE)
if (!file.exists(basedir)) { dir.create(basedir) }
basename <- paste(basedir,"/win-",lwin,"-",win,"-",rwin,sep='')
datafile <- paste( basename ,"-results.RData",sep='')
resultsfile <- paste( basename ,"-results.tsv",sep='')
plotfile <- paste( basename ,"-plot",sep='')

lwin <- rwin <- 1
win <- 2
winlen <- lwin+win+rwin
boundary <- "none"; meanboundary <- 0
gmfile <- paste(paste("genmatrices/genmatrix",winlen,boundary,meanboundary,sep="-"),".RData",sep='')
load(gmfile)

projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=lwin, rwin=rwin )
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE )

# all pairwise counts
counts <- lapply( names( mtCDNA ), function (x) {
        lapply( names( mtCDNA ), function (y) {
                counttrans( rownames(projmatrix), colnames(projmatrix), mtCDNA[[x]],  mtCDNA[[y]],  lwin=lwin ) 
            } ) } )

initcounts <- lapply( counts, rowSums )

save( lwin, rwin, win, winlen, boundary, meanboundary, gmfile, projmatrix, subtransmatrix, counts, initcounts, file=countfile )

