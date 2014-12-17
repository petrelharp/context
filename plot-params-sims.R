#!/usr/bin/Rscript

tsvs <- commandArgs(TRUE)

if (length(tsvs)!=1) {
    cat("Usage:\
        Rscript plot-many-sims.R many-sims-results.tsv \
")
    q()
}

tsv <- tsvs[1]

x <- read.table(tsv,header=TRUE,check.names=FALSE)
pdf(paste0(tsv,".pdf"))
n_params = sum(grepl("sim:", colnames(x)))
stopifnot(n_params == sum(grepl("fit:", colnames(x))))
layout(matrix(1:(4*ceiling(n_params/4)),nrow=4))
# First come the "sim:" parameters, then the "fit:" parameters.
for (k in 1:n_params) {
    hist(x[,n_params+k], xlim=range(x[,n_params+k],x[,k]), main=names(x)[n_params+k]);
    abline(v=unique(x[,k]),col='red');
    if(any(x$longwin==6)) {
        hist(x[x$longwin==6,n_params+k], col=adjustcolor("blue",.5), add=TRUE)
    }
}
dev.off()
