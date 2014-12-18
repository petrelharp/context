#!/usr/bin/Rscript --vanilla
require(optparse)

usage <- "Gather and summarize results of analysis into a single .RData file."

option_list <- list(
        make_option( c("-s","--sim"), type="character", help=".RData file containing simulation." ),
        make_option( c("-f","--fit"), type="character", help=".RData file containing fit." ),
        make_option( c("-o","--outfile"), type="character", default='-', help="File to save things in. [default=stdout]"),
        make_option( c("-j","--json"), type="logical", action="store_true", default=FALSE, help="Save output in json? (else, RData) [default=%default]")
        )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))

source("../context-inference-fns.R")

load(opt$sim)  # provides "simseq.opt"    "simseq.config" "simseqs"
load(opt$fit)  # provides "model"

if (class(model)=="context") {
    time <- as.numeric(simseq.opt$tlen);
    timevec <- c( rep(as.numeric(simseq.opt$tlen),nmuts(model)), rep(1,length(coef(model))-nmuts(model)) )
    sim.params <- c( simseq.config$tip$mutrates*time, simseq.config$tip$selcoef, unlist(simseq.config$tip$fixfn.params) )
    names(sim.params) <- names(coef(model))
    mr_compare <- data.frame(
        fit=coef(model)/timevec,
        simulated=sim.params/timevec,
        stringsAsFactors=FALSE )
} else if (class(model)=="contextTree") {

} else { stop("unrecognized object:", class(model)) }

if ( ! opt$json ) {
    save(simseq.config, simseq.opt, model, mr_compare, file=opt$outfile)
} else {
    require(jsonlite)
    outfile <- openwrite(opt$outfile)
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
                    ), auto_unbox=TRUE, pretty=TRUE )
    cat( paste(json,"\n"), file=outfile )
    flush(outfile)
}


