
#' Shape Fixation Function
#'
#' Total influx of fixation given shape difference ds = (s[to] - s[from]),
#'   which assuming change is BAD, means a selection differential of (-1)*abs(ds)
#'
#' @export
shape.fixfn <- function (ds,Ne,...) {
    if (length(ds)==0) { 1 } else { ifelse( ds==0, 1, Ne*expm1(2*abs(ds))/expm1(2*Ne*abs(ds)) ) }
}

#' Population Genetics Fixation Function
#'
#' Total influx of fixation given selection coefficient (s[to] - s[from]) difference ds
#'
#' @export
popgen.fixfn <- function (ds,Ne,...) {
    if (length(ds)==0) { 1 } else { ifelse( ds==0, 1, Ne*expm1(-2*ds)/expm1(-2*Ne*ds) ) }
}

#' Ising Model Fixation Function
#'
#' The right thing for Gibbs sampling.
#'
#' @export
ising.fixfn <- function (ds,...) { 1/(1+exp(-ds)) }

#' Null Fixation Function
#'
#' Does nothing (always returns 1).
#'
#' @export
null.fixfn <- function (...) { 1 }


