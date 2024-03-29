% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transition_matrix.R
\name{downbranch}
\alias{downbranch}
\title{Collapse down branches, rather than up -- like upbranch, but transposed --}
\usage{
downbranch(genmatrix, rootmatrix, mutrates, selcoef, tlen, ...)
}
\value{
e^( tlen * t(genmatrix(mut,sel)) ) %*% rootmatrix
}
\description{
Collapse down branches, rather than up -- like upbranch, but transposed --
}
