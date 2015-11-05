#!/usr/bin/Rscript --vanilla
library(optparse)

usage <- "Gather and summarize results of analysis into a single .RData file."

option_list <- list(
        make_option( c("-s","--sim"), type="character", help=".RData file containing simulation." ),
        make_option( c("-f","--fit"), type="character", help=".RData file containing fit." ),
        make_option( c("-o","--outfile"), type="character", default='-', help="File to save things in. [default=stdout]"),
        make_option( c("-j","--json"), type="logical", action="store_true", default=FALSE, help="Save output in json? (else, RData) [default=%default]")
        )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))

scriptdir <-  dirname(sub("--file=","",commandArgs()[grep("--file",commandArgs())]))
source(file.path(scriptdir,"context-inference-fns.R"),chdir=TRUE)

load(opt$sim)  # provides "simseq.opt"    "simseq.config" "simseqs"
load(opt$fit)  # provides "model"

if (class(model)=="context") {
    time <- as.numeric(simseq.opt$tlen)
    timevec <- c( rep(as.numeric(simseq.opt$tlen),nmuts(model)), rep(1,length(coef(model))-nmuts(model)) )
    sim.params <- c( simseq.config$tip$mutrates*time, simseq.config$tip$selcoef, unlist(simseq.config$tip$fixfn.params) )
    names(sim.params) <- names(coef(model))
    mr_compare <- data.frame(
        fit=coef(model)/timevec,
        simulated=sim.params/timevec
    )
} else if (class(model)=="contextTree") {
    sim.tlens <- simseq.config$tree$edge.length
    names(sim.tlens) <- edge.labels(simseq.config$tree)
    sim.params <- c( sim.tlens, do.call( c, lapply( names(model@models), function (mname) {
                    c( simseq.config[[mname]]$mutrates, simseq.config[[mname]]$selcoef, unlist(simseq.config[[mname]]$fixfn.params) )
            } ) ) )
    names(sim.params) <- names(coef(model))
    mr_compare <- data.frame(
            fit=coef(model),
            simulated=sim.params
        )
} else { 
    stop(paste("unrecognized object:", class(model))) 
}

if ( ! opt$json ) {
    save(simseq.config, simseq.opt, model, mr_compare, file=opt$outfile)
} else {
    library(jsonlite)
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


