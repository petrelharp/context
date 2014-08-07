#!/usr/bin/Rscript

dirs <- list.files(".",pattern="^0.*")
mle.filenames <- unlist( lapply( dirs, function (x) { list.files(x,"-results.RData$",full.names=TRUE) } ) )
mle.params <- as.numeric( t( sapply( strsplit( unlist(lapply(mle.filenames,sapply,basename)), "-" ), "[", 2:5 ) ) )
dim(mle.params) <- c(length(mle.params)/4,4)
colnames(mle.params) <- c("leftwin","shortwin","rightwin","patlen")
mle.files <- cbind( 
    data.frame( 
        dir=unlist(lapply(mle.filenames,sapply,dirname))
        ),
    mle.params,
    stringsAsFactors=FALSE
    )
mle.files$sample <- factor(gsub("_.*$","",mle.files$dir))
mle.files$frame <- factor(gsub("^[^_]*_","",mle.files$dir))

environs <- lapply( mle.filenames, function (fn) {
        env <- new.env()
        load(fn,envir=env)
        return(env)
    } )

estimates <- lapply( unique(mle.files$patlen), function (k) {
        est <- sapply( environs[mle.files$patlen==k], function (env) {
                with(env, point.estimate)
            } )
        colnames(est) <- mle.files$dir[mle.files$patlen==k]
        return(est)
    } )
allnames <- sort(colnames(estimates[[1]]))

pdf(file="point-estimates.pdf",width=12,height=8,pointsize=10)
layout(t(1:2))
# par(mar=c(max(strwidth(unlist(lapply(estimates,rownames)),"user")),4,2,1)+.1)
plot( estimates[[1]][nrow(estimates[[1]]),], estimates[[2]][nrow(estimates[[2]]),], xlab='gamma scale', ylab='gamma scale' )
abline(0,1)
text( estimates[[1]][nrow(estimates[[1]]),], estimates[[2]][nrow(estimates[[2]]),], labels=allnames )
for (est in estimates) {
    emeans <- rowMeans(est)
    #matplot( est[order(emeans)[-nrow(est)],], type='l' )
    ltys <- as.numeric(mle.files$frame)[match(colnames(est),mle.files$dir)]
    cols <- as.numeric(mle.files$sample)[match(colnames(est),mle.files$dir)]
    matplot( sweep(est[order(emeans)[-nrow(est)],],2,est[nrow(est),],"/"), type='l', xaxt='n', xlab='', ylab='mutation rate', col=cols, lty=ltys )
    axis(1,at=1:(nrow(est)-1),labels=rownames(est)[-nrow(est)],tick=FALSE,line=0,las=2)
    layout(1)
}
legend("topleft",lty=ltys,col=cols,legend=colnames(est))
dev.off()
