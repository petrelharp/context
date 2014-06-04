#!/bin/bash

Rscript -e "fn=commandArgs(TRUE); options(width=1000); z=do.call(rbind, lapply(fn,function(x){ load(x); c(estimates['mle',],win=win,lrwin=lwin) })); rownames(z)=fn; z" $(find . -name "*-results.RData") | sed -e "s/tmut: //g" | sed -e 's/ | /|/g'

echo <<ALLDONE
# paste this into R
x <- read.table(pipe("./quick-results.sh"))
x$sp <- gsub(".*(mm9[a-zA-Z0-9]*).*","\\1",rownames(x))
x <- x[ c( which("mm9rn5"==x$sp), 
    which("mm9oryCun2"==x$sp), 
    which("mm9hg19"==x$sp), 
    which("mm9ornAna1"==x$sp), 
    which("mm9galGal3"==x$sp) ), ]
matplot( x[x$win==3,1:18],type='l', xaxt='n', col=rainbow(24), lty=1 )
axis(1,labels=x$sp[x$win==3],at=1:sum(x$win==3))
matlines( x[x$win==1,1:18], col=rainbow(24), lty=2 )
legend("topright", legend=colnames(x)[1:18], col=rainbow(24)[1:18], lty=1 )
legend("topleft", legend=c("window 3", "window 1"), lty=1:2)
ALLDONE
