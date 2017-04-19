# Just some commands to play with to understand things.
# This is after we've run the workflow.sh script to simulate sequences, etc.

infile <- "sim-cpg-123456.counts"
gmfile <- "genmatrices/genmatrix-4-singlebase.RData"

logfile <- "-"
logcon <- if (logfile=="-") { stdout() } else { file(logfile,open="wt") }
sink(file=logcon, type="message")
sink(file=logcon, type="output")   # send both to log file

library(contextual)
library(contextutils)
library(simcontext)

lwin <- 1

# load generator matrix
stopifnot(file.exists(gmfile))
load(gmfile)  # provides 'genmatrix'

# read in counts
counts <- read.counts(infile,lwin)
stopifnot( all( rownames(counts) == rownames(genmatrix) ) )
print(counts@lwin)
head(counts@counts)

projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=counts@lwin, fpatterns=colnames(counts) )

head(projmatrix)

# get ad hoc initial guesses at parameters
adhoc <- countmuts(counts=counts,mutpats=mutpats,lwin=counts@lwin)
adhoc.mutrates <- adhoc[1,]/adhoc[2,]
adhoc.mutrates <- ifelse( is.finite(adhoc.mutrates) & adhoc.mutrates > 0, adhoc.mutrates, 1e-4 )

head(adhoc)

genmatrix@x

adhoc.mutrates

genmatrix@x <- update(genmatrix,mutrates=adhoc.mutrates,selcoef=c())
# get an initial transition matrix
subtransmatrix <- computetransmatrix( genmatrix, projmatrix, tlen=1, time="fixed")

head(subtransmatrix)

