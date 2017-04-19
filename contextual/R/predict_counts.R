
predictcounts.context <- function (model, longwin=NULL, shortwin=NULL, leftwin=NULL, initcounts=rowSums(model), mutrates=model@mutrates, selcoef=model@selcoef, genmatrix=model@genmatrix, projmatrix=model@projmatrix, params=model@params, tlen=1 ) {
    # default values not cooperating with S4 methods:
    if (is.null(longwin)) { longwin <- longwin(model) }
    if (is.null(shortwin)) { shortwin <- shortwin(model) }
    if (is.null(leftwin)) { leftwin <- leftwin(model) }
    if (!missing(genmatrix) && missing(projmatrix)) {
        projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, fpatterns=getpatterns(shortwin,genmatrix@bases) )
    }
    predictcounts(longwin,shortwin,leftwin,initcounts,mutrates,selcoef,genmatrix,projmatrix,params,tlen)
}


#' Predict Counts from a Model
#'
#' Compute expected counts of paired patterns:
#'  where the actual computation happens
#'
predictcounts <- function (longwin, shortwin, leftwin, initcounts, mutrates, selcoef, genmatrix, projmatrix, params=NULL, tlen=1 ) {
    rightwin <- longwin-shortwin-leftwin
    if (!missing(mutrates)||!missing(selcoef)||!is.null(params)) { genmatrix@x <- do.call( update, c( list(genmatrix,mutrates=mutrates,selcoef=selcoef), params ) ) }
    if (missing(projmatrix)) { projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin, bases=genmatrix@bases ) }
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE, tlen=tlen)
    fullcounts <- initcounts * subtransmatrix
    return( new("tuplecounts", 
            leftwin=leftwin, 
            counts=Matrix(fullcounts), 
            bases=genmatrix@bases,
            colpatterns=data.frame(short=colnames(fullcounts))
        ) )
}

#' Compute Counts for Shorter Tmers.
#'
#'   Valid ranges for parameters are
#'    (l-lc)^+ <= k < (l+w)-(lc+wc)+(r-rc)^-
#'   where 
#'       l = leftwin, lc = new.leftwin
#'       w = shortwin, wc = new.shortwin
#'       r = rightwin, rc = new.rightwin
#' If the original counts were from overlapping windows, then this will overcount the resulting patterns:
#'    if you slide a window of length L in steps of size 1 then a subwindow of size W
#'      will be seen in ( (L-W) ) big windows;
#'    so we need to divide the counts by the factor 'overcount' below
#'      ... but patterns at the boundary of the sequence will not be overcounted. 
#'     Take the ceiling of the resulting counts to fix these.
#'
#' @export
projectcounts <- function( counts, new.leftwin, new.shortwin, new.longwin, overlapping=FALSE ) {
    leftwin <- leftwin(counts)
    longwin <- longwin(counts)
    shortwin <- shortwin(counts)
    rightwin <- longwin-shortwin-leftwin
    new.rightwin <- new.longwin-new.shortwin-new.leftwin
    if ( max(0L,leftwin-new.leftwin) > (leftwin+shortwin)-(new.leftwin+new.shortwin)+min(0L,rightwin-new.rightwin) ) {
        stop("unreconcilable windows specified.")
    }
    pcounts <- matrix(0,nrow=npatterns(new.longwin,counts@bases),ncol=npatterns(new.shortwin,counts@bases))
    for (k in max(0L,leftwin-new.leftwin):((leftwin+shortwin)-(new.leftwin+new.shortwin)+min(0L,rightwin-new.rightwin))) {
        lpmat <- collapsepatmatrix( ipatterns=rownames(counts), leftwin=k, rightwin=longwin-(k+new.longwin), bases=counts@bases )
        rpmat <- collapsepatmatrix( ipatterns=colnames(counts), leftwin=k+new.leftwin-leftwin, rightwin=shortwin-(k+new.leftwin-leftwin+new.shortwin), bases=counts@bases )
        pcounts <- pcounts + t(lpmat) %*% counts %*% (rpmat)
    }
    dimnames(pcounts) <- list( colnames(lpmat), colnames(rpmat) )
    if (overlapping) {
        overcount <- sum( 
                ( (0:(longwin-1))+new.longwin <= longwin ) &
                ( (0:(longwin-1))+new.leftwin >= leftwin ) &
                ( (0:(longwin-1))+new.leftwin+new.shortwin <= leftwin+shortwin ) )
        pcounts <- ( pcounts/overcount )
    }
    return( new("tuplecounts", 
            leftwin=new.leftwin, 
            counts=Matrix(pcounts), 
            bases=counts@bases, 
            colpatterns=data.frame(short=colnames(pcounts))
        ) )
}

#' Compute expected counts of paired patterns
#'
predicttreecounts <- function (shortwin, leftwin=0, rightwin=0, initcounts, mutrates, selcoef, tlens, genmatrix, projmatrix, initfreqs, patcomp, ... ) {
    longwin <- leftwin+shortwin+rightwin
    if (missing(genmatrix)) { genmatrix <- makegenmatrix(patlen=leftwin+shortwin+rightwin,...) }
    if (missing(projmatrix)) { projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin ) }
    if (missing(patcomp) & !is.null(names(initfreqs))) {
        patcomp <- apply( do.call(rbind, strsplit(rownames(genmatrix),'') ), 2, match, names(initfreqs) )  # which base is at each position in each pattern
    }
    patfreqs <- initfreqs[patcomp]
    dim(patfreqs) <- dim(patcomp)
    patfreqs <- apply( patfreqs, 1, prod )
    updownbranch <- getupdowntrans( genmatrix, projmatrix, mutrates=mutrates, selcoef=selcoef, initfreqs=patfreqs, tlens=tlens )
    if (missing(initcounts)) { initcounts <- 1 }
    fullcounts <- initcounts * updownbranch
    dimnames(fullcounts) <- dimnames( projmatrix )
    return( fullcounts )
}


