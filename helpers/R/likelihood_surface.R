
#' @export
gradest <- function (likfun, params, eps=mean(params)/1000) {
    # estimate gradient, crudely
    gradup <- sapply( seq_along(initpar), function (k) { likfun(initpar+ifelse(seq_along(initpar)==k,eps,0)) } )
    graddown <- sapply( seq_along(initpar), function (k) { likfun(initpar-ifelse(seq_along(initpar)==k,eps,0)) } )
    return((gradup-graddown)/(2*eps))
}



#' @export
likelihood.surface <- function (model, tlen=1, 
        plot.it=TRUE, ask=FALSE, progress=TRUE,
        mutrates=model@mutrates, selcoef=model@selcoef, 
        do.parallel=("parallel" %in% .packages(all.available=TRUE)),
        numcores=if (do.parallel) { parallel::detectCores() } else { 1 },
        ngrid=10) {
    # Compute the likelihood on a grid of points about the parameters in 'model'
    # and plot these (note: plots more than one figure)
    this.lapply <- if ( do.parallel && numcores>1 ) { function (...) { parallel::mclapply( ..., mc.cores=numcores ) } } else { lapply }

    timevec <- c( rep(as.numeric(tlen),nmuts(model)), rep(1,length(coef(model))-nmuts(model)) )
    if (!missing(mutrates)) { model@mutrates <- mutrates*tlen }
    if (!missing(selcoef)) { model@selcoef <- selcoef }
    pvec <- coef(model)/timevec
    f <- function (x) { model@likfun((x*timevec)[model@results$use.par]) }

    # contours
    jk <- combn(which(model@results$use.par),2)
    grids <- this.lapply( 1:ncol(jk), function (ii) {
            j <- jk[1,ii]
            k <- jk[2,ii]
            inds <- 1:length(pvec)
            inds[j] <- 1+length(pvec)
            inds[k] <- 2+length(pvec)
            pvals <- seq( .9*min(pvec[j]), 1.1*max(pvec[j]), length.out=ngrid )
            qvals <- seq( .9*min(pvec[k]), 1.1*max(pvec[k]), length.out=ngrid )
            lvals <- matrix( NA, nrow=length(pvals), ncol=length(qvals) )
            for (jj in seq_along(pvals)) {
                for (kk in seq_along(qvals)) {
                    lvals[jj,kk] <- f( c(pvec,pvals[jj],qvals[kk])[inds] )
                }
            }
            if (progress) { cat(".") }
            glist <- list( pvals, qvals, lvals )
            names( glist ) <- c( names(coef(model)[model@results$use.par])[c(j,k)], "lvals" )
            return(glist)
        } )
    if (plot.it) {
        for (ii in seq_along(grids)) {
            j <- jk[1,ii]
            k <- jk[2,ii]
            cols <- colorspace::diverge_hcl(64)
            plot( grids[[ii]][[1]][row(grids[[ii]]$lvals)], grids[[ii]][[2]][col(grids[[ii]]$lvals)], cex=3, pch=20, 
                    col=cols[cut(grids[[ii]]$lvals,breaks=length(cols))], ask=ask,
                    xlab=names(coef(model)[model@results$use.par])[j],
                    ylab=names(coef(model)[model@results$use.par])[k]
                )
            contour( grids[[ii]][[1]], grids[[ii]][[2]], grids[[ii]]$lvals, add=TRUE )
        }
    }
    return( grids )
}

