#!/usr/bin/Rscript
source("../context-inference-fns.R")

resultsfile <- commandArgs(TRUE)[1]
# resultsfile <- "02-C-M_in_frame/win-2-2-2-2-results.RData"
if (length(commandArgs(TRUE))<2) {
    residsfile <- gsub("-results.RData","-resids.tsv",resultsfile)
} else {
    residsfile <- commandArgs(TRUE)[2]
}

load(resultsfile)

resids <- listresids(counts,expected,file=residsfile)
