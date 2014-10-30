#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "\
Write out the table of residuals for a fitted model, in the same format as the count data.
"

option_list <- list(
    # input/output
        make_option( c("-i","--infile"), type="character", help="Input file, .RData with fitted model object."),
        make_option( c("-o","--outfile"), type="character", help="File to save results to.  [default: like infile, but with -resids.(longwin).(shortwin).l(leftwin).tsv appended]"),
        make_option( c("-w","--longwin"), type="integer", help="Size of long window. [default: as in model]" ),
        make_option( c("-s","--shortwin"), type="integer", help="Size of short window. [default: as in model]" ),
        make_option( c("-l","--leftwin"), type="integer", help="Size of offset of short window from the left. [default: as in model]"),
        make_option( c("-p","--pretty"), type="logical", action="store_true", default=FALSE, help="Compute z-score and organize residuals by these?"),
        make_option( c("-m","--gmfile"), type="character", help="File with precomputed generator matrix."),
        make_option( c("-c","--countfile"), type="character", help="Input file with tuple counts, tab-separated, with header 'reference', 'derived', 'count'. (only necessary if longwin is longer than used to fit model)")
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$infile)) { stop("No input file.  Run\n  compute-resids.R -h\n for help.\n") }
if (is.null(opt$basedir)) { opt$basedir <- dirname(opt$infile) }
print(opt) # this will go in the pbs log

source("../context-inference-fns.R")
options(error = print.and.dump)

# and fitted model
stopifnot(file.exists(opt$infile))
load(opt$infile) # provides 'model'

if (is.null(opt$longwin)) { opt$longwin <- longwin(model) }
if (is.null(opt$shortwin)) { opt$shortwin <- shortwin(model) }
if (is.null(opt$leftwin)) { opt$leftwin <- leftwin(model) }

# set outfile name
if (is.null(opt$outfile)) {
    opt$outfile <- paste( opt$basedir, "/", gsub("\\.[^.]*","",basename(opt$infile) ), "-resids", ".", opt$longwin, ".", opt$shortwin, ".l", opt$leftwin, ".tsv", sep='' )
    print(opt$outfile)
}

if (is.null(opt$countfile)) {
    counts = NULL
} else {
    counts = read.counts(opt$countfile, opt$leftwin)
}

residframe = computeresids (model,
    pretty            = opt$pretty,
    in_longwin        = opt$longwin,
    in_shortwin       = opt$shortwin,
    in_leftwin        = opt$leftwin,
    counts            = counts,
    genmatrixfile     = opt$genmatrixfile)

options(scipen=10)
write.table(file=opt$outfile, x=format(residframe,digits=3), sep='\t', quote=FALSE, row.names=TRUE )
