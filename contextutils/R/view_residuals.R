rev.string <- function (x) { sapply(lapply(strsplit(x,''), rev), paste, collapse='') }

#' Collapse reverse complements in residual table
#'
#' @param resids A residual table, with inpat, outpat, observed, expected, resid, and z columns.
#'
#' @export
reverse.complement.resids <- function (resids) {
    rc <- resids
    rc$inpat <- chartr(rev.string(rc$inpat), old="ACGT", new="TGCA")
    rc$outpat <- chartr(rev.string(rc$outpat), old="ACGT", new="TGCA")
    x <- paste0(resids$inpat, resids$outpat)
    y <- paste0(rc$inpat, rc$outpat)
    rc.k <- match(y,x)
    rc <- rc[rc.k,]
    nonself <- (x!=y)
    resids$observed <- resids$observed + rc$observed
    resids$observed[!nonself] <- resids$observed[!nonself]/2
    resids$expected <- resids$expected + rc$expected
    resids$expected[!nonself] <- resids$expected[!nonself]/2
    resids$resid <- resids$observed - resids$expected
    resids$z <- resids$resid/sqrt(nchar(resids$inpat[1])*resids$expected)
    return(resids[rc.k >= seq_along(rc.k),])
}

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
    resids$z <- resids$resid/sqrt(nchar(resids$inpat[1])*as.vector(expected))
    if (is.numeric(trim)) { resids <- subset( resids, is.numeric(resids$z) & (abs(resids$z) > trim) ) }
    resids <- resids[order(resids$z),]
    if (!missing(file)) {
        write.table(file=pipe(paste("column -t >", file)), x=resids, sep=' ',quote=FALSE, row.names=FALSE )
        return(invisible(resids))
    } else {
        return(resids)
    }
}

#' Look for motifs in residuals
#'
#' For the top and bottom sets of `npats` residuals, does the following:
#'  - pastes 'long' and 'short' sequences together
#'  - finds the matrix of hamming distances between these
#'  - uses hclust to cluster the results
#'  - uses hclust to cluster the results
#'
#' @param resids A table of residuals containing 'inpat', 'outpat' and a 'z' score.
#' @param npats The number of top patterns to look for motifs in.
#' @param nclusts The number of clusters to report (using hclust).
#' @param top Whether to do those with large positive z scores? Otherwise, negative.
#' @param return.list Whether to output the full list of patterns in each motif.
#' @param ... Additional arguments passed to print.motif.
#'
#' @export
clusterresids <- function (resids, npats=300, nclusts=12, top=TRUE, return.list=FALSE, ...) {
    if (!top) { resids$z <- (-1)*resids$z }
    trim <- quantile((resids$z[is.finite(resids$z)]),1-npats/sum(is.finite(resids$z)))
    resids <- subset(resids, is.finite(z) & z>trim )
    sclust <- do_clusterresids(resids$inpat, resids$outpat, nclusts=nclusts)
    motifs <- lapply( 1:nclusts, function (k) 
                     print.motif(names(sclust)[sclust==k],print=FALSE,weights=resids$z[sclust==k]^2, ...) )
    # cat( c("", paste(do.call( paste, c( motifs, list(sep="   ") ) ),"\n") ) )
    combined_z <- tapply(1:nrow(resids), sclust, function (k) {
                           observed <- sum(resids[k,"observed"])
                           expected <- sum(resids[k,"expected"])
                           return( (observed-expected)/sqrt(nchar(resids$inpat[1])*expected) )
                     } )
    out <- data.frame(motif=unlist(motifs), 
                      z=combined_z,
                      npats=as.vector(table(sclust)))
    out <- out[order(out$z, decreasing=TRUE),]
    if (return.list) {
        out <- list( summary=out, patterns=lapply(1:nclusts, function(k) names(sclust)[sclust==k]) )
    }
    return(out)
}

do_clusterresids <- function (inpat, outpat, nclusts) {
    longwin <- nchar(inpat[1])
    shortwin <- nchar(outpat[1])
    pats <- paste(inpat,outpat,sep="")
    # 10 secs for 3000x3000
    sdists <- stringdist::stringdistmatrix(pats,pats,method="hamming")
    rownames(sdists) <- colnames(sdists) <- paste(inpat,outpat,sep="|")
    sclust <- cutree( hclust(as.dist(sdists)), k=nclusts )
    return(sclust)
}


#' Like table, but with weights
#'
#' @export
table.weighted <- function(x,weights=rep(1,length(x))) {
    x <- factor(x)
    tx <- table(x)
    return( tx * tapply(weights,x,sum,na.rm=TRUE) )
}

#' Print motifs
#'
#' If not 'long', then positions with a base with total weight above 75% of the total are upper-case,
#' between 50% and 75% are lower-case, otherwise are '.'.
#'
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
                                       (tolower(names(x)[themax]))
                                   } else { "." }
                               )
        } ), collapse='' )
    }
    if (print) { for (x in samps) { cat(x,"\n") } }
    return(invisible(samps))
}
