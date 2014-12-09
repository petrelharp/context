#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "Gather and summarize results of analysis into a single .RData file."

option_list <- list(
        make_option( c("-s","--sim"), type="character", help=".RData file containing simulation." ),
        make_option( c("-f","--fit"), type="character", help=".RData file containing fit." ),
        make_option( c("-o","--outfile"), type="character", help="File to save things in."),
        make_option( c("-R","--RData"), type="logical", default=TRUE, help="Save output in binary? (else, json) [default=%default]")
        )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))

source("../context-inference-fns.R")

load(opt$sim)  # provides "simseq.opt"    "simseq.config" "simseqs" 
load(opt$fit)  # provides "model"

time <- as.numeric(simseq.opt$tlen);
timevec <- c( rep(as.numeric(simseq.opt$tlen),nmuts(model)), rep(1,length(coef(model))-nmuts(model)) )
sim.params <- c( simseq.config$tip$mutrates*time, simseq.config$tip$selcoef, simseq.config$tip$fixfn.params )
names(sim.params) <- names(coef(model))
mr_compare <- data.frame( 
    fit=coef(model)/timevec,
    simulated=sim.params/timevec,
    stringsAsFactors=FALSE )

if ( opt$RData ) {
    save(simseq.config, simseq.opt, model, mr_compare, file=opt$outfile);
} else {
    require(jsonlite)
    json <- toJSON( list(
                    fit.coef = as.list(coef(model)),
                    sim.coef = as.list( sim.params ),
                    fit.params = list( leftwin=leftwin(model),
                            longwin=longwin(model),
                            shortwin=shortwin(model),
                            rowtaxon=rowtaxon(model@counts),
                            loglik=model@results$value,
                            convergence=model@results$convergence
                            )
                    ), auto_unbox=TRUE )
    cat( paste(json,"\n"), file=opt$outfile )
}

