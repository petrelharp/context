#!/usr/bin/Rscript


library(contextual)
library(contextutils)

show(load(commandArgs[1]))

cat("acceptance: ", mrun$accept, "\n")

matplot(mrun$batch,type='l')
legend("topleft",legend=names(mrun$final),lty=1:5,col=1:6)
