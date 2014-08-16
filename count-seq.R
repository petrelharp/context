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
"

option_list <- list(
        make_option( c("-i","--infile"), type="character", help=".RData file containing simulation." ),
        make_option( c("-o","--outfile"), type="character", help="File to direct output to."),
        make_option( c("-w","--longwin"), type="integer", help="Size of long window." ),
        make_option( c("-s","--shortwin"), type="integer", help="Size of short window." ),
        make_option( c("-l","--leftwin"), type="integer", help="Size of offset of short window from the left."),
        make_option( c("-r","--revcounts"), action="store_true", default="FALSE", help="Count reversed?")
        )
countseq.opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(countseq.opt$outfile)) { # default outfile
    outfile <- paste(
                     gsub(".RData","",countseq.opt$infile),
                     if(countseq.opt$revcounts){"-rev"}else{""},
                     ".", countseq.opt$longwin, ".", countseq.opt$shortwin, ".l", countseq.opt$leftwin,
                     ".counts",
                 sep="")
}
attach(countseq.opt)

source("../sim-context-fns.R")
source("../context-inference-fns.R")

load(infile) # provides simseq.opt, bases, mutpats, mutrates, selpats, selcoef, fixfn, seqlen, tlen, initfreqs, simseqs

longpats <- getpatterns(countseq.opt$longwin,bases)
shortpats <- getpatterns(countseq.opt$shortwin,bases)

# this returns a matrix
counts <- if (!revcounts) {
    counttrans( longpats, shortpats, simseqs[[1]]$initseq, simseqs[[1]]$finalseq, leftwin=countseq.opt$leftwin )
} else {
    counttrans( longpats, shortpats, simseqs[[1]]$finalseq, simseqs[[1]]$initseq, leftwin=countseq.opt$leftwin )
}

countframe <- data.frame( reference=rownames(counts)[row(counts)],
                         derived=colnames(counts)[col(counts)],
                         count=as.vector(counts)
                         )

write.table(countframe, file=outfile, row.names=FALSE, sep='\t', quote=FALSE)
