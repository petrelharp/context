#!/usr/bin/Rscript

json.files <- commandArgs(TRUE)

if (length(json.files)==0) {
    cat("Usage:\
        Rscript collect-params-results.R (list of json files) > outfile\
")
    q()
}

library(jsonlite)

jsondata <- lapply( lapply( json.files, fromJSON ), function (x) {
            simc <- paste('sim',names(x[['sim.coef']]),sep=':')
            fitc <- paste('fit',names(x[['fit.coef']]),sep=':')
            those <- setdiff(names(x),c('sim.coef','fit.coef'))
            y <- do.call(cbind,lapply(x[c('sim.coef','fit.coef',those)],as.data.frame))
            names(y) <- c( simc, fitc, unlist( lapply( x[those], names ) ) )
            return(y)
        } )

all.data <- do.call( rbind, jsondata )

write.table( all.data, row.names=FALSE, sep='\t', file=stdout() )
