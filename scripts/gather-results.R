#!/usr/bin/env Rscript
library(optparse)

usage <- "Gather and summarize results of analysis into a single .RData file.
  gather-results.R -h
for help."

option_list <- list(
        make_option( c("-s","--sim"), type="character", help=".RData file containing simulation." ),
        make_option( c("-f","--fit"), type="character", help=".RData file containing fit." ),
        make_option( c("-o","--outfile"), type="character", default='-', help="File to save things in. [default=stdout]"),
        make_option( c("-j","--json"), type="logical", action="store_true", default=TRUE, help="Save output in json? (else, RData) [default=%default]")
        )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))

if (is.null(opt$fit)) { stop(usage) }

library(contextual)
library(contextutils)

if (!is.null(opt$sim)) load(opt$sim)  # provides "simseq.opt"    "simseq.config" "simseqs"
load(opt$fit)  # provides "model"

# mr_compare <- data.frame(fit=coef(model))

if (class(model)=="context" || class(model)=="contextMCMC") {
    if (!is.null(opt$sim))  {
        time <- as.numeric(simseq.opt$tlen)
        timevec <- c( rep(as.numeric(simseq.opt$tlen),nmuts(model)), rep(1,length(coef(model))-nmuts(model)) )
        sim.params <- c( simseq.config$tip$mutrates*time, 
                        simseq.config$tip$selcoef, 
                        unlist(simseq.config$tip$fixfn.params) )
        names(sim.params) <- c( mutnames(simseq.config$tip$mutpats), 
                               selnames(simseq.config$tip$selpats), 
                               names(simseq.config$tip$fixfn.params) )
        # mr_compare$simulated <- sim.params
    }
    # add on posterior quantiles
    if (class(model)=="contextMCMC") {
        quantiles <- lapply( c("q2.5%"=0.025, "q25%"=0.25, "q50%"=0.5, "q75%"=0.75, "q97.5%"=0.975), 
            function (q) {
                out <- rep(NA, length(coef(model)))
                for (k in 1:length(coef(model))) {
                    j <- which(model@results$use.par)[k]
                    out[j] <- quantile(model@results$batch[,k],q)/timevec[j]
                }
                return(out)
            } )
    } else {
        quantiles <- NULL
    }
} else if (class(model)=="contextTree") {
    if (!is.null(opt$sim)) {
        sim.tlens <- simseq.config$tree$edge.length
        names(sim.tlens) <- edge.labels(simseq.config$tree)
        sim.params <- c( sim.tlens, do.call( c, lapply( names(model@models), function (mname) {
                        out <- c( simseq.config[[mname]]$mutrates, simseq.config[[mname]]$selcoef, unlist(simseq.config[[mname]]$fixfn.params) )
                        names(out) <- c( mutnames(simseq.config[[mname]]$mutpats), 
                                           selnames(simseq.config[[mname]]$selpats), 
                                           names(simseq.config[[mname]]$fixfn.params) )
                        out
                } ) ) )
        # mr_compare$simulated <- sim.params
    }
    # add on posterior quantiles
    if (!is.null(model@results$batch)) {
        quantiles <- lapply( c("q2.5%"=0.025, "q25%"=0.25, "q50%"=0.5, "q75%"=0.75, "q97.5%"=0.975), 
            function (q) {
                out <- rep(NA, length(coef(model)))
                for (k in 1:sum(model@results$use.par)) {
                    j <- which(model@results$use.par)[k]
                    out[j] <- quantile(model@results$batch[,k],q)
                }
                return(out)
            } )
    } else {
        quantiles <- NULL
    }
} else { 
    stop(paste("unrecognized object:", class(model))) 
}

# mr_compare <- as.data.frame(c(mr_compare, as.data.frame(quantiles, check.names=FALSE) ))

# if ( ! opt$json ) {
#     save(simseq.config, simseq.opt, model, mr_compare, file=opt$outfile)
# } else {
    library(jsonlite)
    outfile <- openwrite(opt$outfile)
    json <- jsonlite::toJSON( list(
                    fit.coef = as.list(coef(model)),
                    sim.coef = if (is.null(opt$sim)) { NULL } else { as.list( sim.params ) },
                    fit.params = list( leftwin=leftwin(model),
                            longwin=longwin(model),
                            shortwin=shortwin(model),
                            rowtaxon=rowtaxon(model@counts),
                            loglik=model@results$value,
                            convergence=model@results$convergence
                            ),
                    posterior.quantiles = quantiles
                    ), auto_unbox=TRUE, pretty=TRUE, digits=8 )
    cat( paste(json,"\n"), file=outfile )
    flush(outfile)
# }


