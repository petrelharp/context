#' Make a Likelihood Function
#'
#' Returns a log likelihood function appropriate for optimization,
#' which is a function of a single vector, concatenated, of:
#'   mutrates, selcoef, fixparams
#'
#' The likelihood is the ratio of likelihoods with (long,short) 
#' and (long,short-1) window lengths.
#'
#' The length of time is fixed to be 1.
#'
#' @param longwin The length of the longest window used in Tmers.
#' @param shortwin The length of the short window used in Tmers.
#' @param leftwin The offset of short window used.
#' @param config A list of arguments to `makegenmatrix()` except patlen.
#' @param seq The sequence to compute the likelihood of.
#'
#' @details Note this returns log-likelihood, *not* negative log-likelihood.
#'
#' @return A function.
#'
#' @examples
#'
#' config <- list(mutpats=list(list(c("XX","OO")),list(c("X","O"),c("O","X"))), selpats=list(), bases=c("X","O"), fixfn=null.fixfn, selfactors=list())
#' seq <- list( 
#'   initseq=
#'      "OOXXXXOOXXOXOXOXOOXOOOXOXOXXOOXXOOOOOOXOXOXXXXXXXOXOXOOXXXXOXOXXOOOXXXXOXOOXOXOOOXOXOXOOXXXOOOXOXOXXXOOOOXXXXOOOXOOOOXOOOXXOXXOXOOXXOOXXXOXOXOXXOOOXOOOOXOOXXXXOXXXXOOXOOOOXOOXXOXOXXOOXXXOXOXXXOOOOOOOX",
#'   finalseq=
#'      "OOXXOXOXXOXOXXOXOOOXOOXOOXXOOOXXXOOOOOOOXXXXXXXXOOXXOXOXXXOOXXXOXOOXOXXXXOOXOOXOOOXXOOXOXOXOXOXOXOXXXOOOOXXXOOXOOOOOXOXOOXXOXXOXOOXOXOXXOXXOOXXXOOOXOOOOOXOXXXXOXXXXOOOOXOOXOOXXOXOXXOOXXXOOXXXOXOOOOOOX")
#' lfn <- likelihood.function(5,3,1,config,seq)
#' # arguments to lfn are: the XX -> OO mutation rate, and the X <-> O mutation rate
#' with(environment(lfn), mutnames(genmatrix@mutpats))
#' lfn(c(1.5,2.0))
#'
#' @export likelihood.function
likelihood.function <- function (
                                 longwin,
                                 shortwin,
                                 leftwin,
                                 config,
                                 seq
                 ) {
    stopifnot(shortwin>1 && shortwin<=longwin)
    longpats <- getpatterns(longwin, config$bases)
    shortpats <- getpatterns(shortwin, config$bases)
    shorterpats <- getpatterns(shortwin-1, config$bases)
    genmatrix <- do.call( makegenmatrix, c( list(patlen=longwin), config ) )
    # for the numerator
    counts_1 <- counttrans.list( list(longpats, shortpats), 
        seqlist=seq[c("initseq", "finalseq")],
        leftwin=leftwin, bases=config$bases )
    projmatrix_1 <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin(counts), fpatterns=colnames(counts) )
    # for the denominator
    counts_2 <- counttrans.list( list(longpats, shorterpats),
        seqlist=seq[c("initseq", "finalseq")],
        leftwin=leftwin, bases=config$bases )
    projmatrix_2 <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=leftwin(counts), fpatterns=colnames(counts) )
    # the parameters
    params <- rep(1.0, nmuts(genmatrix) + nsel(genmatrix) + length(fixparams(genmatrix)))
    # and the function
    llfun <- function (params) {
        fparams <- params[seq( 1+nmuts(genmatrix)+nsel(genmatrix), length.out=length(fixparams(genmatrix)) )]
        names(fparams) <- fixparams(genmatrix)
        genmatrix@x <- do.call( update_x, c( list( G=genmatrix,mutrates=params[1:nmuts(genmatrix)],selcoef=params[seq(1+nmuts(genmatrix),length.out=nsel(genmatrix))]), as.list(fparams) ) )
        subtransmatrix_1 <- computetransmatrix( genmatrix, projmatrix_1, tlen=1.0, time="fixed")
        subtransmatrix_2 <- computetransmatrix( genmatrix, projmatrix_2, tlen=1.0, time="fixed")
        num <- sum( counts_1@counts * log(subtransmatrix_1) )
        denom <- sum( counts_2@counts * log(subtransmatrix_2) )
        return(num - denom)
    }

    return(llfun)
}
