#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Write out the table of residuals for a fitted model, in the same format as the count data.
"

option_list <- list(
    # input/output
        make_option( c("-i","--infile"), type="character", help="Input file, .RData with fitted model object."),
        make_option( c("-o","--outfile"), type="character", help="File to save results to.  [default: like infile, but with -resids.(longwin).(shortwin).l(leftwin).tsv appended]"),
        make_option( c("-w","--longwin"), type="integer", help="Size of long window." ),
        make_option( c("-s","--shortwin"), type="integer", help="Size of short window." ),
        make_option( c("-l","--leftwin"), type="integer", help="Size of offset of short window from the left."),
        make_option( c("-m","--gmfile"), type="character", help="File with precomputed generator matrix."),
        make_option( c("-c","--countfile"), type="character", help="Input file with tuple counts, tab-separated, with header 'reference', 'derived', 'count'. (only necessary if longwin is longer than used to fit model)")
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$infile)) { stop("No input file.  Run\n  compute-resids.R -h\n for help.\n") }
if (is.null(opt$basedir)) { opt$basedir <- dirname(opt$infile) }
if (is.null(opt$outfile)) { opt$outfile <- paste( opt$basedir, "/", gsub("\\.[^.]*","",basename(opt$infile) ), "-resids", ".", longwin, ".", shortwin, ".l", leftwin, ".tsv", sep='' ) }
print(opt) # this will go in the pbs log
options(error = quote({dump.frames(to.file = TRUE); q()}))

if (interactive()) {
    opt$infile <- "sim-cpg-123456-modelfit-54321.RData"
    opt$countfile <- "sim-cpg-123456.5.3.l1.counts"
    opt$leftwin <- 2
    opt$longwin <- 5
    opt$shortwin <- 1
    opt$gmfile <- "genmatrices/genmatrix-5-singlebase.RData"
}

attach(opt)

source("../context-inference-fns.R")

# load generator matrix
stopifnot(file.exists(gmfile))
load(gmfile)  # provides 'genmatrix'

# and fitted model
stopifnot(file.exists(infile))
load(infile) # provides 'model'

# get counts
if (is.null(opt$countfile) && ( opt$longwin > longwin(model) || opt$shortlen > shortlen(model) ) ) {
    stop("If window lengths are longer than fitted model, then need to supply counts.")
} else if (!is.null(opt$countfile)) {
    counts <- read.counts(opt$countfile,opt$leftwin)
    if ( ( opt$longwin > longwin(counts) || opt$shortwin > shortwin(counts) ) ) {
        stop("Supplied counts use a window that is too short.")
    }
    counts <- projectcounts( counts, opt$leftwin, opt$shortwin, opt$longwin-opt$leftwin-opt$shortwin )
} else {
    counts <- model@data
}
stopifnot( all( rownames(counts) == rownames(genmatrix) ) )

fitted <- fitted( model, longwin=opt$longwin, shortwin=opt$shortwin, leftwin=opt$leftwin, initcounts=rowSums(counts), genmatrix=genmatrix )

resids <- residuals( model, counts=counts, genmatrix=genmatrix )
rr <- counts-fitted

