% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patterns.R
\name{collapsepatmatrix}
\alias{collapsepatmatrix}
\title{Construct a Pattern-Substring Matrix}
\usage{
collapsepatmatrix(ipatterns, leftwin,
  shortwin = nchar(as.character(fpatterns[1])),
  rightwin = nchar(as.character(ipatterns[1])) - shortwin - leftwin,
  fpatterns = getpatterns(nchar(as.character(ipatterns[1])) - leftwin -
  rightwin, bases), bases)
}
\arguments{
\item{ipatterns}{Long patterns index the rows of the result.}

\item{leftwin}{Left offset of short patterns from long patterns.}

\item{shortwin}{Length of short patterns.}

\item{rightwin}{Right offset of short patterns from long patterns.}

\item{fpatterns}{Short patterns, which index the columns of the result.}

\item{bases}{Alphabet of bases.}
}
\value{
A (nbases)^k by (nbases)^{k-leftwin-rightwin} sparse matrix
projection matrix that maps patterns onto the shorter patterns obtained by
deleting leftwin characters at the start and rightwin characters at the end.
}
\description{
Construct the matrix U described in the tex.
ipatterns are the "input" patterns, while fpatterns are the "final" projected patterns
This function assumes that all input patterns are the same length.
}
\examples{
bases <- c("X","O")
# columns are the right-hand character
collapsepatmatrix(ipatterns=getpatterns(2,bases), leftwin=1, fpatterns=bases)
# columns are the middle two characters
collapsepatmatrix(ipatterns=getpatterns(4,bases), leftwin=1, rightwin=1, bases=bases)

}
