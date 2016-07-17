#!/usr/bin/Rscript

source("../context-inference-fns.R",chdir=TRUE)

show(load(commandArgs[1]))

cat("acceptance: ", mrun$accept, "\n")

matplot(mrun$batch,type='l')
legend("topleft",legend=names(mrun$final),lty=1:5,col=1:6)
