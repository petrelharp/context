#!/usr/bin/Rscript

library(contextual)
library(simcontext)
library(testthat)

set.seed(23)

bases <- c("X","O")

context("Projecting pattern counts.")

# Short sequences:
tlen <- .5
seqlen <- 1000
mutrates <- c(1,1)
selcoef <- c(-.5,.5)

simseqs <- list( 
                initseq=Biostrings::BString("OXXXXXXXOXXOXOOXOOXXXXOOXXXXXOXOXOXOXOOOXXXXOXXOXOXXXXOXXXXOXOOXOXXXOOOOXXOOOOOXXOXXXXOXOOOXOXOXXOXOOOXOOOOXXXOOXXOXOXXOXXXOOOOXXXXXOXOXOOOOXXOXXOOXXXXOXXXOXOOXXOXOOOXOXOOXXXXOOXOOXOXOXOOOOOXXXXOXXOOXXXXXOOOXXOXOXOOOOOXXXXOOXXXOOOOOOXXXXXXOXXXOXXXXXXXXXXOXXXOOXOXOOXOOOOXXXOOXXOXOOOOOXXOXXOXOXXOXOOOXOOOXOOOXOOXOXOXOOOXXOOOOOOOOOXOXOOOOOXOOXXOXOOOOOOOOOOXOXXOOXOXXXXXXOOXXOOOXXOOXXOOOOOXXXOOXXXOXOXOXOOOOXXXXXXXXXOXOXOXOOXOXXXXXXXOXXXOXXXXOOXXOOXXXXOOXXOOOXXOXOOOXXOXOXOXOXXOOXXOXOOXXXXOXXOXXOOXOOOOXXXOOXOOXOOXXOOXOOXXXXXOOXXOXXXXOXXOXXOXXXXXOOXXOXOOOOXOXOOXOOOOOXXXOXOXXOOOOOOXXXOOXXXOOOOOXOXXOOOXOOOOXOOXXXXOXOOXOOOXXOOOXXXXXOXOOXOXOXOXXOOXXOOXOXXXXOOXXOXOXOXOXXXOXXXXXOOOOXXXXXOOOOOXXXXOOXXXOOOOXOXOOXOXXXOOXOXXXOXOXXOXXXOOXXOXXOXOXOOXOOXXXOOXOOOOOXXXXOXXOOOOXXOXXXOXOOOXXOOXXOOXXXXXOOOXOOOXXOXOOOXXOXXOOOOOOOOOOOOOOXOXXOXXOOOXXXXOOXXXXOXOOXOXOXOOOXXOXXXOXOXXOOOXOXOOXXOXOXXXXOXOXOXOOXXOXOOOOXOOXOXOXOXXXOXOOXOXXOXOOXXOOXXOOXOOXOOOOOXOOOXXXOOOOXOOXXOOXOXXOOXOOXXOOOXXOXXXXOXXXOXXXXXOOXXOOOOOOOOXO"),
                finalseq=Biostrings::BString("XXXOXXXXOXXXOOXXOOXXXXOOOOXOXXXOXOXOOOOXXXXXOXXOXXXXXXXOXXXOXXOXOXOXOOXOXXOOOOOXOXXXXXXXOXOXOXOXXOXOOOXOOOOXXXOOOXOXXXXOOXXOXOOXOOXXOXOXXOOXOOOXXOXXXXXXXXOOXOOXXXOOOOXOXOOXXXXOOXOXXOOOOOOOOOOXXXOXXOOXXXXXXOOXXOXOXOOOOOXXXXOOOXOXOOOXOXXXXXXXXOOOXXXXXXXXXXXXXXOOXOOXOXOOXXOXXXXXXXOOOOOXXXXXOOXOXXOXOOOOXOXXOXXXOOXXXOXOXOXXOOOOOXOOOXXXOOOOOOOOXXXOXOOOOOOOOOXXXXXOXXXXXXXXOOXXOOXXXOOXXOOXOOOXXOOOXXOXXXOOOXOOXXXXXXXXOXOOXOXOOXXOXXXXXXOXXXOXXXXOOXXOOXXXXOOXXOOOOXXXOOOXXOXOXXOXXXOOXXXXXXXXXXXXXXXXOOOOOOOXXXOOXOOXOOXXOOOOOXXXXXOOXXXXXXXOXXOOXOXXXXOOOXXXXOOOOXOXOOXOXOXOXXXXXOXXXOOOXXXXXOXXXXOXXXOOOXXOOOXOOOOXOOXXXXOXXOXOOOXXOOXXXXXXOXOOXOXOXOXOOOOXXXXOXXXXOOXXOXOOOOXXXXOXXOXXXOOOXXXXXOOXXXXXXXOXXXOOOXOXOOXXXOXXXOXXOXOXOXOXXOXXXOXXXOXXOXOXOOOOOXXXOOXOOXXOXXXXOXOOOOOOXOXXXXXOXOXOOOXXXOXXXXOOOOXOOOXXOXOOXXOOXXOOOOXXOOOOOXOOOOXXOXXXOOXXXXXOXXXXOOOOXOXOXOXXXXOXXXOOXXOXOXXXXXOXOXXOXXXXXXXOOXOOXXOXOOOOOXXXXXOXOOXXXXXXXOXXOXOOXXOOXXOXXOOXOOOOXXOOOXXXXOOXXOOOXOOOOOXOOOOOOXOOXXXOXXXXOXXXOXXXXXXOXXOXOOOOOOXO"),
           maxrate = 0.817574476193644, ntrans = NULL, mutpats = list( list(c("O", "X")), list(c("X", "O"))), 
           selpats = list( c("OX", "XO"), "X"), mutrates = c(1, 1), selcoef = c(-0.5, 0.5), 
           tlen = 0.5, seqlen = 1000, bases = c("X", "O") )

long <- 6
short <- 2

longpats <- getpatterns(long,bases=bases)
shortpats <- getpatterns(short,bases=bases)
# 6-6-0 counts
long.counts <- simcontext::counttrans( ipatterns=longpats, fpatterns=longpats, simseqs=simseqs, leftwin=0 )
# 2-2-0 counts
short.counts <- simcontext::counttrans( ipatterns=shortpats, fpatterns=shortpats, simseqs=simseqs, leftwin=0 )
# 6-6-0 to 2-2-0
proj.counts <- projectcounts( counts=long.counts, new.leftwin=0, new.longwin=short, new.shortwin=short, overlapping=TRUE )
# 6-2-2 counts
long.short.counts <- simcontext::counttrans( ipatterns=longpats, fpatterns=shortpats, simseqs=simseqs, leftwin=2 )
# 6-2-2 to 2-2-0
proj.counts.2 <- projectcounts( counts=long.short.counts, new.leftwin=0, new.longwin=short, new.shortwin=short, overlapping=TRUE )

all.counts <- data.frame( inpat=rownames(proj.counts@counts)[row(proj.counts@counts)], 
    outpat=colnames(proj.counts@counts)[col(proj.counts@counts)], 
    proj=as.vector(proj.counts@counts),
    proj2=as.vector(proj.counts.2@counts),
    stringsAsFactors=FALSE
    )
all.counts$short <- short.counts@counts[ cbind( match(all.counts$inpat,rownames(short.counts@counts)), match(all.counts$outpat,colnames(short.counts@counts)) ) ]

inseq <- as.character(simseqs$initseq)
inseq.vec <- substring( inseq, first=1:(seqlen-short+1), last=short:seqlen )
outseq <- as.character(simseqs$finalseq)
outseq.vec <- substring( outseq, first=1:(seqlen-short+1), last=short:seqlen )

all.counts$truth <- sapply( 1:nrow(all.counts), function (k) {
        sum( grepl( all.counts$inpat[k], inseq.vec ) & grepl( all.counts$outpat[k], outseq.vec ) )
    } )

test_that("Counting patterns two different ways", {
    expect_true( all( all.counts$short == all.counts$truth ) )
    expect_true( sum(all.counts$proj)+long-short == sum(all.counts$short) )
    expect_true( sum(all.counts$proj2)+long-short == sum(all.counts$short) )
    expect_true( all( ( all.counts$short - all.counts$proj ) <= 1 ) )
    expect_true( all( ( all.counts$short - all.counts$proj2 ) <= 1 ) )
})


###############
context("Projecting tree pattern counts.")

# set.seed(23)
# config <- list(
#     "tree" = c( "(sp1 : 1.8, (sp2 : 1.2, sp3 : 1.0) an1 : 0.5 ) root;" ), 
#     "bases" = c( "X", "O" ),
#     "initfreqs" = c( 0.5, 0.5),
#     "sp1" = list(
#         "mutpats" = list(
#             list( c("X", "O" ) ),
#             list( c("O", "X" ) ),
#             list( c("XX", "XO" ), c( "OO", "XO" ) )
#         ),
#         "mutrates" = c( 0.1, 0.1, 0.15 )
#     ),
#     "sp2" = "sp1",
#     "sp3" = "sp1",
#     "an1" = "sp1"
# )
# config <- parse.models(treeify.config(config))
# simseqs <- simseq.tree( seqlen=100, config=config)

simseqs <- list(
                sp1=list( initseq=Biostrings::BString("XOOXXOXXXXXXOOXOXXXXXOOXOXXXOXXOOOOXOXXOXOOXXXOOOXOXOXXOOXXXOOXOOOOXOXXXXXOOXOXXXOXXXOOXOOOOXXOOXXOO"),
                           finalseq=Biostrings::BString("XXXXXOXOOOXXOOOXOXXOXOOXOXXOOOXOOOXOOOXOXOOXXOXOOXXXOOXXOXXXXOOOXXOOOOXXXXOOXOOXXOXXXXOXXOXOXXOOXOOO")),
                sp2=list( initseq=Biostrings::BString("XXOXXOXXOXXXOOXOXXXOXOOXOXXXOXXOOOOXOXXOXOOXXXXOOXOXOXXOOXXXOOOOOOOXOXXXXOOOXOXXOOXXXOOXOOOOXXOOXXOO"),
                           finalseq=Biostrings::BString("XXOXXOXXOXOXXOXOXXXOXOOXOXXXOXXOOXOXOXOOXXOOXXXOOOOXOXOXXXOXXOOOOOOXOXXXOOOOXOXXOOXXXOOXOXOXXXOOXXXO")),
                sp3=list( initseq=Biostrings::BString("XXOXXOXXOXXXOOXOXXXOXOOXOXXXOXXOOOOXOXXOXOOXXXXOOXOXOXXOOXXXOOOOOOOXOXXXXOOOXOXXOOXXXOOXOOOOXXOOXXOO"),
                           finalseq=Biostrings::BString("XOOOXOXXOXOXXXXOXOXOXXOXOOXXOXOOXOOXOXOOXOOXXXXOXXOOOXXOOXXOXOOOOXOOOXXXXOOOXOXXOOXXXOOXOOOXXXOOXXOO")),
                an1=list( initseq=Biostrings::BString("XOOXXOXXXXXXOOXOXXXXXOOXOXXXOXXOOOOXOXXOXOOXXXOOOXOXOXXOOXXXOOXOOOOXOXXXXXOOXOXXXOXXXOOXOOOOXXOOXXOO"),
                           finalseq=Biostrings::BString("XXOXXOXXOXXXOOXOXXXOXOOXOXXXOXXOOOOXOXXOXOOXXXXOOXOXOXXOOXXXOOOOOOOXOXXXXOOOXOXXOOXXXOOXOOOOXXOOXXOO")),
              root=list( finalseq=Biostrings::BString("XOOXXOXXXXXXOOXOXXXXXOOXOXXXOXXOOOOXOXXOXOOXXXOOOXOXOXXOOXXXOOXOOOOXOXXXXXOOXOXXXOXXXOOXOOOOXXOOXXOO")) )


long <- 3
short <- 2

bases <- c("X","O")
seqlen <- 100

theseqs <- list( sp1=as.character(simseqs$sp1$finalseq),
                sp2=as.character(simseqs$sp2$finalseq),
                sp3=as.character(simseqs$sp3$finalseq))
theseqvecs <- lapply(1:long, function (k) lapply( theseqs, substring, first=1:(seqlen-k+1), last=k:seqlen ) )

the.truth <- function (cf) {
    sp1sv <- theseqvecs[[nchar(cf$sp1[1])]]$sp1
    sp2sv <- theseqvecs[[nchar(cf$sp2[1])]]$sp2
    sp3sv <- theseqvecs[[nchar(cf$sp3[1])]]$sp3
    sapply( 1:nrow(cf), function (k) {
               matchvec <- grepl( cf$sp1[k], sp1sv )
               matchvec <- matchvec & grepl( cf$sp2[k], sp2sv )
               matchvec <- matchvec & grepl( cf$sp3[k], sp3sv )
               sum(matchvec)
            } )
}

longpats <- getpatterns(long,bases=bases)
long.colpatterns <- expand.grid( getpatterns(long,bases), getpatterns(long,bases), stringsAsFactors=FALSE )
colnames(long.colpatterns) <- c("sp2", "sp3")

long.counts <- simcontext::counttrans.tree( rowpatterns=longpats, rowtaxon="sp1", colpatterns=long.colpatterns, simseqs=simseqs, leftwin=0, bases=bases )
long.countframe <- countframe(long.counts)
long.countframe$truth <- the.truth(long.countframe)

shortpats <- getpatterns(short,bases=bases)
short.colpatterns <- expand.grid( getpatterns(short,bases), getpatterns(short,bases), stringsAsFactors=FALSE )
colnames(short.colpatterns) <- c("sp2", "sp3")

short.counts <- simcontext::counttrans.tree( rowpatterns=shortpats, rowtaxon="sp1", colpatterns=short.colpatterns, simseqs=simseqs, leftwin=0, bases=bases )
short.countframe <- countframe(short.counts)
short.countframe$truth <- the.truth(short.countframe)

test_that("Counting tree patterns:", {
    expect_equal( long.countframe$truth, long.countframe$count)
    expect_equal( short.countframe$truth, short.countframe$count)
} )

# Project counts:

proj.counts <- projectcounts.tree( counts=long.counts, new.leftwin=0, new.colpatterns=short.colpatterns, new.longwin=short, overlapping=TRUE )

test_that("Projecting tree counts", {
    expect_equal( dimnames(proj.counts@counts), dimnames(short.counts@counts) )
    expect_true( all( abs(as.numeric(proj.counts@counts) - as.numeric(short.counts@counts)) <= 1 ) )
    expect_true( sum(proj.counts@counts) == sum(short.counts@counts) - (long-short) )
} )
