#!/usr/bin/env Rscript
library(optparse)

usage <- "\
Count Tmers from the Rdata file containing simulated (or other) sequence,\
of the following form: \
    WWWWWWWWW \
    llMMMMrrr \
where the  length of the W's is 'longwin', the length of the l's is 'leftwin', and the length of the 'M's is 'shortwin'.
 \
If --revcounts, count in the other direction, i.e. \
    llMMMMrrr \
    WWWWWWWWW \
\
If counts are on a tree (a named list of sequences),\
need to also specify the name of the clade that gets the 'long' pattern, by e.g.\
  --longclade=sp1 \
. \
 \
If --RData is not given, outputs a tab-separated file with column headers \
in which the first line, begun with a '#', gives what the shift ('leftwin') is. \
\
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", help=".RData file containing simulation." ),
        make_option( c("-o","--outfile"), type="character", help="File to direct output to."),
        make_option( c("-w","--longwin"), type="integer", help="Size of long window." ),
        make_option( c("-s","--shortwin"), type="integer", help="Size of short window." ),
        make_option( c("-l","--leftwin"), type="integer", help="Size of offset of short window from the left."),
        make_option( c("-p","--shift"), type="integer", default=0, help="Shift windows by this many sites, zero means completely overlapping (same as 1, but writes out a single file). [default=%default]"),
        make_option( c("-c","--longclade"), type="character", help="Which clade, gets the 'long' patterns. [default: root for a stick tree]"),
        make_option( c("-t","--shortclades"), type="character", help="Which clades, gets the 'short' patterns. [default: all tips except longclade]"),
        make_option( c("-R","--RData"), type="logical", action="store_true", default=FALSE, help="Output as tuplecounts object in RData? [default: as csv]")
        )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))

library(contextual)
library(simcontext)
library(contextutils)

load(opt$infile) # provides simseq.opt, simseq.config, and simseqs; config has bases, mutpats, mutrates, selpats, selcoef, fixfn, seqlen, tlen, initfreqs, 

longpats <- getpatterns(opt$longwin,simseq.config$bases)
shortpats <- getpatterns(opt$shortwin,simseq.config$bases)

if (length(simseqs)==1) { # stick tree
    if (is.null(opt$longclade)) { opt$longclade <- simseq.config$tree$node.label }
    if (is.null(opt$shortclades)) { opt$shortclades <- simseq.config$tree$tip.label }
} else {
    if (is.null(opt$longclade)) { 
        if ('root' %in% names(simseqs)) {
            opt$longclade <- "root"
        } else {
            stop("Must specify a long clade name.") 
        }
    }
    if (is.null(opt$shortclades)) { opt$shortclades <- setdiff(simseq.config$tree$tip.label,opt$longclade) }
}

if (is.null(opt$outfile)) { # default outfile
    opt$outfile <- paste(
                     gsub(".RData","",opt$infile),
                     "-", opt$longwin, 
                     "-", opt$longclade,
                     "-", opt$shortwin, 
                     "-", paste(opt$shortclades,collapse="-"),
                     "-l", opt$leftwin,
                     "-shift", opt$shift,
                     ".counts",
                 sep="")
}

counts <- counttrans.list( list(longpats,shortpats)[c(1,rep.int(2,length(opt$shortclades)))], 
    simseqs=simseqs[c(opt$longclade,opt$shortclades)], 
    leftwin=opt$leftwin, bases=simseq.config$bases,
   shift=opt$shift )
if (opt$RData) {
    outfile <- gsub("\\.counts$",".RData",opt$outfile)
    cat("Writing to:", outfile, "\n")
    save( counts, file=outfile )
} else if (opt$shift==0) {
    outfile <- opt$outfile
    cat("Writing to:", outfile, "\n")
    cframe <- countframe( counts )
    cat( paste('#', jsonlite::toJSON( list( leftwin=opt$leftwin, longwin=opt$longwin, shortwin=opt$shortwin, shift=opt$shift, offset=0 ), auto_unbox=TRUE ), "\n" ), file=outfile )
    write.table(cframe, file=opt$outfile, row.names=FALSE, sep='\t', quote=FALSE, append=TRUE)
} else {
    for (k in 1:opt$shift) {
        cframe <- countframe( counts[[k]] )
        outfile <- paste(opt$outfile,k,sep='.')
        cat("Writing to:", outfile, "\n")
        cat( paste('#', jsonlite::toJSON( list( leftwin=opt$leftwin, longwin=opt$longwin, shortwin=opt$shortwin, shift=opt$shift, offset=k ), auto_unbox=TRUE ), "\n" ), file=outfile )
        # cat( paste('# { "leftwin" : ', opt$leftwin, '}\n', sep=''), file=outfile )
        write.table(cframe, file=outfile, row.names=FALSE, sep='\t', quote=FALSE, append=TRUE)
    }
}
