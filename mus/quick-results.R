#!/bin/Rscript


fn <- scan(pipe("find mm* -name \"*-results.RData\""),what='char')
z <- lapply(fn,function(x){ load(x); c(as.matrix(estimates)['mle',],win=win,lrwin=lwin) })
zn <- factor( gsub(" ","", sapply( lapply(z,names), paste, collapse='' ) ) )
names(z) <- fn
zz <- tapply( seq_along(z), zn, function (k) { y <- do.call( rbind, z[k] ); rownames(y) <- fn[k]; data.frame(y) } )
for (k in seq_along(zz)) {
    zz[[k]]$sp <- gsub(".*(mm9[a-zA-Z0-9]*).*","\\1",rownames(zz[[k]]))
    zz[[k]] <- zz[[k]][ c( which("mm9rn5"==zz[[k]]$sp), 
        which("mm9oryCun2"==zz[[k]]$sp), 
        which("mm9hg19"==zz[[k]]$sp), 
        which("mm9ornAna1"==zz[[k]]$sp), 
        which("mm9galGal3"==zz[[k]]$sp) ), ]
    names(zz[[k]]) <- gsub(".ab.1",".ba",gsub("..","->",gsub("tmut..","",names(zz[[k]]),fixed=TRUE),fixed=TRUE),fixed=TRUE)
    rownames(zz) <- NULL
}

x <- zz[[1]]
layout(t(1:2),width=c(3,1))
matplot(x[,1:13],type='l',xaxt='n', lty=1,col=rainbow(24))
matlines(x[,13+(1:13)],type='l',lty=2,col=rainbow(24))
axis(1,labels=x$sp[x$win==1],at=1:sum(x$win==1))
legend("topright", legend=colnames(x)[1:13], col=rainbow(24), lty=1 )
matplot(x[,26+(1:4)],type='l',xaxt='n', lty=1,col=rainbow(6))
axis(1,labels=x$sp[x$win==1],at=1:sum(x$win==1))
legend("topright", legend=colnames(x)[26+(1:4)], col=rainbow(6), lty=1 )

x <- zz[[2]]
matplot( x[x$win==1,1:18],type='l', xaxt='n', col=rainbow(24), lty=1 )
axis(1,labels=x$sp[x$win==1],at=1:sum(x$win==1))
if (any(x$win==3)) { matlines( x[x$win==3,1:18], col=rainbow(24), lty=2 ) }
legend("topright", legend=colnames(x)[1:18], col=rainbow(24)[1:18], lty=1 )
legend("topleft", legend=c("window 1", "window 3"), lty=1:2)
