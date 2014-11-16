#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "Gather and summarize results of analysis into a single Rmd file."

option_list <- list(
        make_option( c("-s","--sim"), type="character", help=".RData file containing simulation." ),
        make_option( c("-f","--fit"), type="character", help=".RData file containing fit." ),
        make_option( c("-o","--outfile"), type="character", help="File to save things in.")
        )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))

source("../context-inference-fns.R")

load(opt$sim)
load(opt$fit)

time <- as.numeric(simseq.opt$tlen);
mr_compare <- data.frame(
    fit=c( model@mutrates / time, model@selcoef, model@params),
    stringsAsFactors=FALSE );
simvalues <- data.frame(
        names=c(
            mutnames(simseq.config$tip$mutpats),
            selnames(simseq.config$tip$selpats),
            names(simseq.config$tip$fixfn.params) ),
        simulated=unlist(c(
            simseq.config$tip$mutrates,
            simseq.config$tip$selcoef,
            simseq.config$tip$fixfn.params ) )
);
mr_compare$simulated <-
    simvalues$simulated[match(rownames(mr_compare),simvalues$names)];

save(simseq.config, simseq.opt, model, mr_compare,
    file=opt$outfile);

