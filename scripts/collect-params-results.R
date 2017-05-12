#!/usr/bin/env Rscript

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
                x <- tryCatch( {
                    fromJSON(jf)            
                }, error=function (cond) { 
                    message(sprintf("Can't read from %s.\n",jf))
                    return(NULL)
                } )
                if (is.null(x)) { return(NULL) }
                simc <- if (length(x[['sim.coef']])>0) { paste('sim',names(x[['sim.coef']]),sep=':') } else { character(0) }
                fitc <- if (length(x[['fit.coef']])>0) { paste('fit',names(x[['fit.coef']]),sep=':') } else { character(0) }
                those <- setdiff(names(x),c('sim.coef','fit.coef','posterior.quantiles'))
                the.names <- character(0)
                if (length(x[['sim.coef']])>0) { the.names <- c('sim.coef') }
                the.names <- c(the.names, c('fit.coef', those))
                y <- do.call(cbind,lapply(x[the.names],adf))
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
jsondata <- jsondata[sapply(jsondata,length)>0]

thenames <- unique(unlist(lapply(jsondata,names)))

all.data <- do.call( rbind, lapply(jsondata, function (x) adf(x, thenames) ))

write.table( all.data, row.names=FALSE, sep='\t', file=stdout() )
