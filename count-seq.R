#!/usr/bin/Rscript --vanilla
require(optparse)

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
  --longclade=sp1
.
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", help=".RData file containing simulation." ),
        make_option( c("-o","--outfile"), type="character", help="File to direct output to."),
        make_option( c("-w","--longwin"), type="integer", help="Size of long window." ),
        make_option( c("-s","--shortwin"), type="integer", help="Size of short window." ),
        make_option( c("-l","--leftwin"), type="integer", help="Size of offset of short window from the left."),
        make_option( c("-c","--longclade"), type="character", help="Which clade, in a tree, gets the 'long' patterns."),
        make_option( c("-r","--revcounts"), action="store_true", default="FALSE", help="Count reversed?")
        )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$outfile)) { # default outfile
    outfile <- paste(
                     gsub(".RData","",opt$infile),
                     if(opt$revcounts){"-rev"}else{""},
                     ".", opt$longwin, ".", opt$shortwin, ".l", opt$leftwin,
                     ".counts",
                 sep="")
}

source("../sim-context-fns.R")
source("../context-inference-fns.R")

load(opt$infile) # provides simseq.opt, simseq.config, and simseqs; config has bases, mutpats, mutrates, selpats, selcoef, fixfn, seqlen, tlen, initfreqs, 

longpats <- getpatterns(opt$longwin,simseq.config$bases)
shortpats <- getpatterns(opt$shortwin,simseq.config$bases)

# this returns a matrix
counts <- if (!revcounts) {
    counttrans( longpats, shortpats, simseqs[[1]]$initseq, simseqs[[1]]$finalseq, leftwin=opt$leftwin )
} else {
    counttrans( longpats, shortpats, simseqs[[1]]$finalseq, simseqs[[1]]$initseq, leftwin=opt$leftwin )
}

countframe <- data.frame( reference=rownames(counts)[row(counts)],
                         derived=colnames(counts)[col(counts)],
                         count=as.vector(counts)
                         )

write.table(countframe, file=outfile, row.names=FALSE, sep='\t', quote=FALSE)
