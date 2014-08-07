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
        make_option( c("-p","--pretty"), type="logical", action="store_true", default=FALSE, help="Compute z-score and organize residuals by these?"),
        make_option( c("-m","--gmfile"), type="character", help="File with precomputed generator matrix."),
        make_option( c("-c","--countfile"), type="character", help="Input file with tuple counts, tab-separated, with header 'reference', 'derived', 'count'. (only necessary if longwin is longer than used to fit model)")
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$infile)) { stop("No input file.  Run\n  compute-resids.R -h\n for help.\n") }
if (is.null(opt$basedir)) { opt$basedir <- dirname(opt$infile) }
print(opt) # this will go in the pbs log
options(error = quote({dump.frames(to.file = TRUE); q()}))

source("../context-inference-fns.R")

# and fitted model
stopifnot(file.exists(opt$infile))
load(opt$infile) # provides 'model'

if (is.null(opt$longwin)){ opt$longwin <- longwin(model) }
if (is.null(opt$shortwin)){ opt$shortwin <- shortwin(model) }
if (is.null(opt$leftwin)){ opt$leftwin <- leftwin(model) }

# set outfile name
if (is.null(opt$outfile)) { 
    opt$outfile <- paste( opt$basedir, "/", gsub("\\.[^.]*","",basename(opt$infile) ), "-resids", ".", opt$longwin, ".", opt$shortwin, ".l", opt$leftwin, ".tsv", sep='' ) 
    print(opt$outfile)
}

# load generator matrix, if needed
if (!is.null(opt$gmfile)) {
    stopifnot(file.exists(opt$gmfile))
    load(opt$gmfile)  # provides 'genmatrix'
} else {
    genmatrix <- model@genmatrix
}

# get counts
if (is.null(opt$countfile) && ( opt$longwin > longwin(model) || opt$shortwin > shortwin(model) ) ) {
    stop("If window lengths are longer than fitted model, then need to supply counts.")
} else if (!is.null(opt$countfile)) {
    counts <- read.counts(opt$countfile,opt$leftwin)
    if ( ( opt$longwin > longwin(counts) || opt$shortwin > shortwin(counts) ) ) {
        stop("Supplied counts use a window that is too short.")
    }
} else {
    counts <- model@data
}
if ( (opt$longwin < longwin(counts)) || (opt$shortwin < shortwin(counts)) ) {
    counts <- projectcounts( counts, opt$leftwin, opt$shortwin, opt$longwin-opt$leftwin-opt$shortwin )
}
stopifnot( all( rownames(counts) == rownames(genmatrix) ) )

expected <- fitted( model, longwin=opt$longwin, shortwin=opt$shortwin, leftwin=opt$leftwin, initcounts=rowSums(counts), genmatrix=genmatrix )

if (opt$pretty) {
    # data frame with columns for long pattern, short pattern, observed, expected, residual, z-score
    residframe <- data.frame( inpat=rownames(counts)[row(counts)],
                        outpat=colnames(counts)[col(counts)],
                        observed=as.vector(counts),
                        expected=as.vector(expected),
                        resid=as.vector(counts)-as.vector(expected),
                        stringsAsFactors=FALSE
                    )
    residframe$z <- residframe$resid/sqrt(as.vector(expected))
    residframe <- residframe[order(residframe$z),]
    write.table(file=opt$outfile, x=residframe, sep=' ', quote=FALSE, row.names=FALSE )
} else {
    # concise matrix
    resids <- residuals( model, counts=counts, genmatrix=genmatrix )
    residframe <- as.vector(resids@counts)
    dim(residframe) <- dim(resids)
    dimnames(residframe) <- dimnames(resids)
    write.table(file=opt$outfile, x=residframe, sep=' ', quote=FALSE, row.names=TRUE )
}



