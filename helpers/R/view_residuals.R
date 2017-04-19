
#' List residuals in order
#'
#' Make a readable output ordered by z-score
#'  optionally writing results out to 'file'
#'  and trimming to only patterns with z-score above 'trim'
#'
#' @export
listresids <- function (counts, expected, file, trim=20, leftwin=(nchar(rownames(counts)[1])-nchar(colnames(counts)[1]))/2) {
    resids <- data.frame( inpat=rownames(counts)[row(counts)],
                         outpat=colnames(counts)[col(counts)],
                         observed=as.vector(counts),
                         expected=as.vector(expected),
                    resid=as.vector(counts-expected) )
    longwin <- nchar(rownames(counts)[1])
    resids$inpat <- with(resids, paste( tolower(substr(inpat,1,leftwin)), substr(inpat,leftwin+1,longwin-leftwin), tolower(substr(inpat,longwin-leftwin+1,longwin)), sep='' ) )
    resids$outpat <- paste(resids$outpat)
    resids$z <- resids$resid/sqrt(as.vector(expected))
    if (is.numeric(trim)) { resids <- subset( resids, is.numeric(resids$z) & (abs(resids$z) > trim) ) }
    resids <- resids[order(resids$z),]
    if (!missing(file)) {
        write.table(file=pipe(paste("column -t >", file)), x=resids, sep=' ',quote=FALSE, row.names=FALSE )
        return(invisible(resids))
    } else {
        return(resids)
    }
}

#' @export
clusterresids <- function (resids,npats=300,nclusts=12) {
    # look for motifs in resids (as above)
    invisible( lapply( c(+1,-1), function (sign) {
            resids$z <- resids$z * sign
            trim <- quantile((resids$z[is.finite(resids$z)]),1-npats/sum(is.finite(resids$z)))
            resids <- subset(resids, is.finite(z) & z>trim )
            longwin <- nchar(resids$inpat[1])
            shortwin <- nchar(resids$outpat[1])
            pats <- paste(resids$inpat,resids$outpat,sep="")
            # 10 secs for 3000x3000
            sdists <- stringdistmatrix(pats,pats,method="hamming")
            rownames(sdists) <- colnames(sdists) <- paste(resids$inpat,resids$outpat,sep="|")
            sclust <- cutree( hclust(as.dist(sdists)), k=nclusts )
            motifs <- lapply( 1:nclusts, function (k) print.motif(names(sclust)[sclust==k],print=FALSE,weights=sqrt(resids$z[sclust==k])) )
            cat( c("", paste(do.call( paste, c( motifs, list(sep="   ") ) ),"\n") ) )
            return(motifs)
        } ) )
}


#' Like table, but with weights
#'
#' @export
table.weighted <- function(x,weights=rep(1,length(x))) {
    x <- factor(x)
    tx <- table(x)
    return( tx * tapply(weights,x,sum,na.rm=TRUE) )
}

#' @export
print.motif <- function (pats,weights=1,n=24,print=TRUE,long=FALSE) {
    pats <- do.call( rbind, strsplit(pats,"") )
    freqs <- lapply( 1:ncol(pats), function(k) {
                    x <- table.weighted(pats[,k],weights); x[order(names(x))]
                } )
    if (long) {
        samps <- do.call( cbind, lapply( freqs, function (x) { names(x)[cut((1:n)/(n+1),breaks=c(0,cumsum(x))/sum(x))] } ) )
        samps <- apply(samps,1,paste,collapse='')
    } else {
        samps <- paste( sapply(freqs, function (x) {
                               themax <- which.max(x)
                               return(
                                   if (max(x)>.75*sum(x)) {
                                       (toupper(names(x)[themax]))
                                   } else if (max(x)>.5*sum(x)) {
                                       (toupper(names(x)[themax]))
                                   } else { "." }
                               )
        } ), collapse='' )
    }
    if (print) { for (x in samps) { cat(x,"\n") } }
    return(invisible(samps))
}
