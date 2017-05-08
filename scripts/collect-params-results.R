#!/usr/bin/Rscript

json.files <- commandArgs(TRUE)

if (length(json.files)==0) {
    cat("Usage:\
        Rscript collect-params-results.R (list of json files) > outfile\
")
    q()
}

library(jsonlite)

# fill in empty entries with NA
adf <- function (z, znames=names(z)) { for (k in znames) { if (length(z[[k]])==0) { z[[k]] <- NA } }; as.data.frame(z, check.names=FALSE) }

jsondata <- lapply( json.files, function (jf) {
                tryCatch( {
                    x <- fromJSON(jf)            
                }, error=function (cond) { 
                    message(sprintf("Reading %s.\n",jf))
                    stop(cond)
                } )
                simc <- paste('sim',names(x[['sim.coef']]),sep=':')
                fitc <- paste('fit',names(x[['fit.coef']]),sep=':')
                those <- setdiff(names(x),c('sim.coef','fit.coef','posterior.quantiles'))
                y <- do.call(cbind,lapply(x[c('sim.coef','fit.coef',those)],adf))
                names(y) <- c( simc, fitc, unlist( lapply( x[those], names ) ) )
                if (!is.null(x$posterior.quantiles)) {
                    for (k in seq_along(x$posterior.quantiles)) {
                        names(x$posterior.quantiles[[k]]) <- names(x[['fit.coef']])
                    }
                    pq <- unlist( x$posterior.quantiles )
                    y <- as.data.frame(c(y, pq), check.names=FALSE)
                }
                return(c( list(file=jf, ctime=format(file.info(jf)$ctime, "%Y-%m-%d_%H-%M-%S")), y))
            } )

thenames <- unique(unlist(lapply(jsondata,names)))

all.data <- do.call( rbind, lapply(jsondata, function (x) adf(x, thenames) ))

write.table( all.data, row.names=FALSE, sep='\t', file=stdout() )
