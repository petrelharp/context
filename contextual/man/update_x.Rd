% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genmatrix.R
\name{update_x}
\alias{update_x}
\title{Update the Entries of a Generator Matrix With New Rates.}
\usage{
update_x(G, mutrates, selcoef, ...)
}
\arguments{
\item{G}{A genmatrix object.}

\item{mutrates}{A new vector of mutation rates.}

\item{selcoef}{A new vector of selection coefficients (length zero if not needed).}

\item{...}{Additional arguments passed to the fixation function, `G@fixfn(s,...)`}
}
\value{
A numeric vector of the same length as G@x.
}
\description{
Calculate the new entries of `x` from muttrans et al.
use like: genmatrix@x <- update_x(genmatrix,...)
}
\details{
The returned value is computed as:
  fixfn(seltrans %*% selcoef) * (muttrans %*% mutrates)
which works as `seltrans` and `muttrans` have been precomputed appropriately.
}
