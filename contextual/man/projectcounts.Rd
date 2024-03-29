% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_counts.R
\name{projectcounts}
\alias{projectcounts}
\alias{projectcounts.tree}
\title{Compute Counts for Shorter Tmers from longer ones.}
\usage{
projectcounts(counts, new.leftwin = leftwin(counts),
  new.shortwin = shortwin(counts), new.longwin = longwin(counts),
  overlapping = TRUE)

projectcounts.tree(counts, new.colpatterns, new.longwin = longwin(counts),
  new.leftwin = leftwin(counts), overlapping = TRUE)
}
\arguments{
\item{counts}{A tuplecounts object.}

\item{new.leftwin}{The left overhang of the new T-mer.}

\item{new.shortwin}{The short end of the new T-mer.}

\item{new.longwin}{The long end of the new T-mer.}

\item{overlapping}{Are the original counts of overlapping T-mers?}

\item{new.colpatterns}{Patterns counted in columns of the new tuplecounts object.}
}
\value{
A new tuplecounts object.
}
\description{
Valid ranges for parameters are
   (l-lc)^+ <= k < (l+w)-(lc+wc)+(r-rc)^-
where 
      l = leftwin,  lc = new.leftwin
      w = shortwin, wc = new.shortwin
      r = rightwin, rc = new.rightwin
}
\details{
If the original counts were from overlapping windows, then this will overcount the resulting patterns:
   if you slide a window of length L in steps of size 1 then a subwindow of size W
     will be seen in ( (L-W) ) big windows;
   so we divide the counts by the corresponding factor 
     ... but patterns at the boundary of the sequence will not be overcounted;
   after division these will become noninteger.  There is not sufficient
   information in a table of counts to fix this problem.

projectcounts.tree differs only in that the projection is specified using a new `coltaxa`,
rather than just left/short/long window lengths.
}
\section{Functions}{
\itemize{
\item \code{projectcounts.tree}: Project tree counts
}}

\examples{
counts <- simcontext::counttrans(ipatterns=getpatterns(3,bases=c("X","O")), fpatterns=getpatterns(2,bases=c("X","O")), 
                        initseq="XOOXXXOXOXOXOOOXXO", finalseq="XOXOXOXOXXOXOOOXXO", leftwin=1)
projectcounts(counts, new.longwin=2, new.leftwin=0) 
projectcounts(counts, new.shortwin=1)

}
