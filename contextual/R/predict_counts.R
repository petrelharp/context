#' @include context_classes.R

#' Predict counts for a 'context' model
predictcounts.context <- function (model, longwin=NULL, shortwin=NULL, leftwin=NULL, 
                                   initcounts=rowSums(model), mutrates=model@mutrates, 
                                   selcoef=model@selcoef, genmatrix=model@genmatrix, 
                                   projmatrix=model@projmatrix, params=model@params, tlen=1 ) {
    # default values not cooperating with S4 methods:
    if (is.null(longwin)) { longwin <- longwin(model) }
    if (is.null(shortwin)) { shortwin <- shortwin(model) }
    if (is.null(leftwin)) { leftwin <- leftwin(model) }
    if (longwin > longwin(model) && ( missing(initcounts) || missing(genmatrix) )) {
        stop("If predicting longer counts than the model was fit under, must provide genmatrix, and initcounts.")
    }
    if (missing(projmatrix) && ( longwin != longwin(model) || shortwin != shortwin(model) || leftwin != leftwin(model))) {
        projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, fpatterns=getpatterns(shortwin,genmatrix@bases) )
    }
    predictcounts(longwin=longwin,
                  shortwin=shortwin,
                  leftwin=leftwin,
                  initcounts=initcounts,
                  mutrates=mutrates,
                  selcoef=selcoef,
                  genmatrix=genmatrix,
                  projmatrix=projmatrix,
                  params=params,
                  tlen=tlen)
}

#' Predict counts for a 'contextTree' model
predictcounts.contextTree <- function (
                                   model, 
                                   rowtaxon=model@counts@rowtaxon, 
                                   coltaxa=colnames(model@counts@colpatterns),
                                   longwin=NULL, shortwin=NULL, leftwin=NULL, 
                                   initcounts=rowSums(model), 
                                   genmatrices=lapply(model@models, function (x) x@genmatrix),
                                   projmatrix=model@models[[1]]@projmatrix, 
                                   initfreqs=model@initfreqs,
                                   params=model@params, tlen=1 ) {
    # default values not cooperating with S4 methods:
    if (is.null(longwin)) { longwin <- longwin(model) }
    if (is.null(shortwin)) { shortwin <- shortwin(model) }
    if (is.null(leftwin)) { leftwin <- leftwin(model) }
    if (longwin > longwin(model) && ( missing(initcounts) || missing(genmatrices) )) {
        stop("If predicting longer counts than the model was fit under, must provide genmatrices, and initcounts.")
    }
    if (missing(projmatrix) && ( longwin != longwin(model) || shortwin != shortwin(model) || leftwin != leftwin(model))) {
        projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrices[[1]]), leftwin=leftwin, fpatterns=getpatterns(shortwin,genmatrices[[1]]@bases) )
    }
    config <- lapply(model@models, function (x) { 
                         list(mutrates=x@mutrates, selcoef=x@selcoef, params=x@params) 
                       } )
    predictcounts.tree(tree=model@tree, config=config, modelnames=model@modelnames, 
                       longwin=longwin, shortwin=shortwin, leftwin=leftwin, 
                       rowtaxon=rowtaxon, coltaxa=coltaxa,
                       initcounts=initcounts, genmatrices=genmatrices,
                       projmatrix=projmatrix, initfreqs=initfreqs)

}


#' Predict Counts from a Model
#'
#' Compute expected counts of paired patterns, either given the state at the root or not.
#'
#' @param longwin Length of the long end of the Tmer.
#' @param shortwin Length of the short end of the Tmer.
#' @param leftwin Length of the left overhang of the Tmer.
#' @param initcounts Number of occurrences of each intial (`long`) pattern. 
#'        Set to NULL to average over all states (see details).
#' @param mutrates Mutation rates optionally used to update `genmatrix`.
#' @param selcoef Selection coefficients optionally used to update `genmatrix`.
#' @param genmatrix The generator matrix used to predict counts (must correspond to a pattern length of `longwin`).
#' @param projmatrix The projection matrix used to project from long to short patterns (constructed if not supplied).
#' @param params Other parameters optionally used to update `genmatrix`.
#' @param tlen Length of time for evolution.
#' 
#' @details In a tree model, if `initcounts` is present the expected counts reported are 
#' conditioned on the given number of counts at the `rowtaxon`.  If
#' `initcounts` is of length 1, then expectations are given for that number of
#' patterns, averaging over everything (including the state at the root).
#' 
#' @return A tuplecounts object.
#' 
#' @export
predictcounts <- function (longwin, 
                           shortwin, 
                           leftwin, 
                           initcounts, 
                           mutrates, 
                           selcoef, 
                           genmatrix, 
                           projmatrix, 
                           params=NULL, 
                           tlen=1 ) {
    if (longwin != longwin(genmatrix)) { stop("Generator matrix does not have the requested longwin.") }
    if (length(initcounts) != nrow(genmatrix)) { stop("initcounts not of the correct length.") }
    rightwin <- longwin-shortwin-leftwin
    if (!missing(mutrates)||!missing(selcoef)||!is.null(params)) { genmatrix@x <- do.call( update_x, c( list(genmatrix,mutrates=mutrates,selcoef=selcoef), params ) ) }
    if (missing(projmatrix)) { projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin, bases=genmatrix@bases ) }
    subtransmatrix <- computetransmatrix( genmatrix, projmatrix, names=TRUE, tlen=tlen)
    fullcounts <- initcounts * subtransmatrix
    return( new("tuplecounts", 
            leftwin=leftwin, 
            counts=Matrix(fullcounts), 
            bases=genmatrix@bases,
            colpatterns=data.frame(short=colnames(fullcounts), stringsAsFactors=FALSE)
        ) )
}

#' Compute Counts for Shorter Tmers from longer ones.
#'
#' Valid ranges for parameters are
#'    (l-lc)^+ <= k < (l+w)-(lc+wc)+(r-rc)^-
#' where 
#'       l = leftwin,  lc = new.leftwin
#'       w = shortwin, wc = new.shortwin
#'       r = rightwin, rc = new.rightwin
#'
#' If the original counts were from overlapping windows, then this will overcount the resulting patterns:
#'    if you slide a window of length L in steps of size 1 then a subwindow of size W
#'      will be seen in ( (L-W) ) big windows;
#'    so we divide the counts by the corresponding factor 
#'      ... but patterns at the boundary of the sequence will not be overcounted;
#'    after division these will become noninteger.  There is not sufficient
#'    information in a table of counts to fix this problem.
#'
#' projectcounts.tree differs only in that the projection is specified using a new `coltaxa`,
#' rather than just left/short/long window lengths.
#'
#' @param counts A tuplecounts object.
#' @param new.leftwin The left overhang of the new T-mer.
#' @param new.shortwin The short end of the new T-mer.
#' @param new.longwin The long end of the new T-mer.
#' @param overlapping Are the original counts of overlapping T-mers?
#'
#' @return A new tuplecounts object.
#'
#' @examples
#' counts <- simcontext::counttrans(ipatterns=getpatterns(3,bases=c("X","O")), fpatterns=getpatterns(2,bases=c("X","O")), 
#'                         initseq="XOOXXXOXOXOXOOOXXO", finalseq="XOXOXOXOXXOXOOOXXO", leftwin=1)
#' projectcounts(counts, new.longwin=2, new.leftwin=0) 
#' projectcounts(counts, new.shortwin=1)
#'
#' @export
projectcounts <- function(counts, 
                          new.leftwin=leftwin(counts), 
                          new.shortwin=shortwin(counts), 
                          new.longwin=longwin(counts), 
                          overlapping=TRUE ) {
    leftwin <- leftwin(counts)
    longwin <- longwin(counts)
    shortwin <- shortwin(counts)
    rightwin <- longwin-shortwin-leftwin
    new.rightwin <- new.longwin-new.shortwin-new.leftwin
    if ( new.shortwin <= 0 || max(0L,leftwin-new.leftwin) > (leftwin+shortwin)-(new.leftwin+new.shortwin)+min(0L,rightwin-new.rightwin) ) {
        stop(sprintf("unreconcilable windows specified: %d,%d,%d to %d,%d,%d.", longwin, shortwin, leftwin, new.longwin, new.shortwin, new.leftwin))
    }
    new.colpatterns <- do.call( 
            expand.grid, 
            list(getpatterns(new.shortwin, bases=counts@bases),
                    stringsAsFactors=FALSE)[c(rep(1,ncol(counts@colpatterns)),2)])
    colnames(new.colpatterns) <- colnames(counts@colpatterns)
    return( projectcounts.tree(
            counts=counts,
            new.colpatterns=new.colpatterns,
            new.longwin=new.longwin,
            new.leftwin=new.leftwin,
            overlapping=overlapping) )
}

#' @param new.colpatterns Patterns counted in columns of the new tuplecounts object.
#' 
#' @describeIn projectcounts Project tree counts
#' @export
projectcounts.tree <- function(
                          counts, 
                          new.colpatterns,
                          new.longwin=longwin(counts),
                          new.leftwin=leftwin(counts), 
                          overlapping=TRUE ) {
    leftwin <- leftwin(counts)
    longwin <- longwin(counts)
    shortwin <- max(nchar(as.character(unlist(counts@colpatterns[1:2,]))))
    rightwin <- longwin-shortwin-leftwin
    colpatterns <- counts@colpatterns
    stopifnot(all( colnames(colpatterns)==colnames(new.colpatterns) ))
    new.shortwin <- max(nchar(as.character(unlist(new.colpatterns[1:2,]))))
    new.rightwin <- new.longwin-new.shortwin-new.leftwin
    new.cps <- lapply(1:ncol(new.colpatterns), function (k) unique( new.colpatterns[,k] ))
    if ( max(0L,leftwin-new.leftwin) > (leftwin+shortwin)-(new.leftwin+new.shortwin)+min(0L,rightwin-new.rightwin) ) {
        stop("unreconcilable windows specified.")
    }
    pcounts <- matrix(0,nrow=npatterns(new.longwin,counts@bases),ncol=nrow(new.colpatterns))
    for (k in max(0L,leftwin-new.leftwin):((leftwin+shortwin)-(new.leftwin+new.shortwin)+min(0L,rightwin-new.rightwin))) {
        lpmat <- collapsepatmatrix( ipatterns=rownames(counts), leftwin=k, rightwin=longwin-(k+new.longwin), bases=counts@bases )
        rpmat <- collapsepatmatrix( ipatterns=colpatterns[,1],
                                    fpatterns=new.cps[[1]], 
                                    leftwin=k+new.leftwin-leftwin)[,match(new.colpatterns[,1],new.cps[[1]])]
        for (j in (1:ncol(new.colpatterns))[-1]) {
            # times for AND
            rpmat <- rpmat * collapsepatmatrix( ipatterns=colpatterns[,j],
                                    fpatterns=new.cps[[j]], 
                                    leftwin=k+new.leftwin-leftwin)[,match(new.colpatterns[,j],new.cps[[j]])]
        }
        pcounts <- pcounts + Matrix::t(lpmat) %*% counts %*% (rpmat)
    }
    dimnames(pcounts) <- list( colnames(lpmat), apply(new.colpatterns,1,paste,collapse=".") )
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
            rowtaxon=counts@rowtaxon,
            colpatterns=new.colpatterns
        ) )
}

#' Compute expected counts of paired patterns on a two-taxon tree.
#'
#' @export
predicttreecounts <- function (shortwin, 
                               leftwin=0, 
                               rightwin=0, 
                               initcounts, 
                               mutrates, 
                               selcoef, 
                               tlens, 
                               genmatrix, 
                               projmatrix, 
                               initfreqs, 
                               patcomp, 
                               ... ) {
    longwin <- leftwin+shortwin+rightwin
    if (missing(genmatrix)) { genmatrix <- makegenmatrix(patlen=leftwin+shortwin+rightwin,...) }
    if (missing(projmatrix)) { projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin, rightwin=rightwin ) }
    if (missing(patcomp) & !is.null(names(initfreqs))) {
        patcomp <- apply( do.call(rbind, strsplit(rownames(genmatrix),'') ), 2, match, names(initfreqs) )  # which base is at each position in each pattern
    }
    patfreqs <- initfreqs[patcomp]
    dim(patfreqs) <- dim(patcomp)
    patfreqs <- apply( patfreqs, 1, prod )
    transmat <- getupdowntrans( genmatrix, projmatrix, mutrates=mutrates, selcoef=selcoef, initfreqs=patfreqs, tlens=tlens )
    if (!missing(initcounts)) { 
        transmat <- initcounts * transmat
    }
    dimnames(transmat) <- dimnames( projmatrix )
    return( transmat )
}

#' @param tree The tree.
#' @param config A named list of model configurations.
#' @param modelnames A list of names of models named by nodes and tips on the tree.
#' @param rowtaxon Which taxon to count long patterns in; indexed by rows of the result.
#' @param coltaxa The taxa to count in the columns.
#' @param genmatrices A named list of genmatrix objects, with names equal to the model names in `config`.
#' @param initfreqs The initial frequences of the (long) patterns at the root.
#'
#' @describeIn predictcounts Computed expected counts on a tree
#' @export
predictcounts.tree <- function (tree,
                                config,
                                modelnames,
                                longwin,
                                shortwin,
                                leftwin, 
                                rowtaxon,
                                coltaxa,
                                initcounts, 
                                genmatrices, 
                                projmatrix, 
                                initfreqs
                                ) {
    if (any(sapply(genmatrices,longwin)<longwin)) { 
        stop("Generator matrices do not have sufficiently long longwin.") 
    }
    bases <- genmatrices[[1]]@bases
    longpats <- rownames(genmatrices[[1]])
    # put the parameters in to the genmatrices
    for (mod in names(genmatrices)) {
        genmatrices[[mod]]@x <- do.call( update_x, 
                                        c( list( G=genmatrices[[mod]],
                                             mutrates=config[[modelnames[[mod]]]]$mutrates,
                                             selcoef=config[[modelnames[[mod]]]]$selcoef), 
                                           as.list(config[[modelnames[[mod]]]]$params) ) )
    }
    # calculate the distribution at the root
    initfreq.index  <- product.index( longpats=longpats, bases=bases ) # which base is at each position in each pattern
    root.distrn <- get.root.distrn( initfreqs, initfreq.index )
    if (missing(projmatrix)) { 
        projmatrix <- collapsepatmatrix(ipatterns=longpats, leftwin=leftwin, 
                                        rightwin=longwin-shortwin-leftwin, bases=bases) 
    }
    # calculate the transition matrix
    transmatrix <- peel.transmat( tree=tree, rowtaxon=rowtaxon, coltaxa=coltaxa, modelnames=modelnames, 
                                  genmatrices=genmatrices, projmatrix=projmatrix, 
                                  root.distrn=root.distrn, debug=TRUE, return.list=FALSE, 
                                  normalize=FALSE )
    if (length(initcounts)>1) {
        transmatrix <- (initcounts/rowSums(transmatrix)) * transmatrix
    } else {
        transmatrix <- initcounts * transmatrix
    }
    colpatterns <- data.frame(do.call(rbind, strsplit(colnames(transmatrix), ",")), 
                              stringsAsFactors=FALSE, check.names=FALSE)
    colnames(colpatterns) <- coltaxa
    out <- new("tuplecounts", 
            leftwin=leftwin, 
            counts=Matrix(transmatrix), 
            bases=bases,
            colpatterns=colpatterns
          )
    return(projectcounts(out, new.longwin=longwin, new.leftwin=leftwin, new.shortwin=shortwin))
}
