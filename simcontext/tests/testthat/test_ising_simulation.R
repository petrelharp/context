library(testthat)

library(contextual)
library(simcontext)

set.seed(23)

bases <- c("X","O")

patlen <- 2
mutpats <- list( 
    list( c("O","X") ),
    list( c("X","O") )
    ) 
mutrates <- c(3,5)
selpats <- list(
        c("OX","XO"),
        c("X")
    )
selcoef <- c(-2,1)


fixfn <- function (ds,...) { 1/(1+exp(-ds)) }



######
# Simulates correctly?


#  Long sequences:
tlen <- .1
seqlen <- 1000/tlen
mutrates <- c(1,1)
selcoef <- c(-.5,.5)
simseqs <- simseq( seqlen, tlen, patlen=patlen, mutpats=mutpats, selpats=selpats, mutrates=mutrates, selcoef=selcoef, bases=bases, fixfn=fixfn, quiet=TRUE )

# counts
left.win <- 2; right.win <- 2; short.win <- 2
long.win <- left.win+short.win+right.win
genmatrix <- makegenmatrix( patlen=long.win, mutpats=mutpats, selpats=selpats, mutrates=mutrates*tlen, selcoef=selcoef, bases=bases, fixfn=fixfn )
projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=left.win, rightwin=right.win, bases=genmatrix@bases )
counts <- counttrans( rownames(projmatrix), colnames(projmatrix), simseqs=simseqs, leftwin=left.win, bases=genmatrix@bases )
cwin <- 2
subcounts <- projectcounts( counts, new.shortwin=cwin, new.leftwin=0, new.longwin=cwin )

all.genmats <- list (
            "simple"=makegenmatrix( patlen=long.win, mutpats=mutpats, selpats=selpats, mutrates=mutrates*tlen, selcoef=selcoef, boundary="none", bases=bases, fixfn=fixfn ),
            "wrap"=makegenmatrix( mutpats=mutpats, selpats=selpats, patlen=long.win, mutrates=mutrates*tlen, selcoef=selcoef, boundary="wrap", bases=bases, fixfn=fixfn ),
            "mean"=meangenmatrix( leftwin=1, rightwin=1, patlen=long.win, mutpats=mutpats, selpats=selpats, mutrates=mutrates*tlen, selcoef=selcoef, boundary="wrap", bases=bases, fixfn=fixfn )
        )

all.expected <- lapply( all.genmats, function (genmatrix) {
        expected <- predictcounts( longwin=long.win, shortwin=short.win, leftwin=left.win, initcounts=rowSums(counts), mutrates=tlen*mutrates, selcoef=selcoef, genmatrix=genmatrix, projmatrix=projmatrix )
        subexpected <- projectcounts( counts=expected, new.leftwin=0, new.shortwin=cwin, new.longwin=cwin )
        return(list( expected=expected, subexpected=subexpected ) )
    } )

all.resids <- lapply( all.expected, function (x) {
            list( full = counts@counts-x[["expected"]]@counts,
                sub = subcounts@counts-x[["subexpected"]]@counts
                ) } )

vars <- sapply( all.resids, function (x) as.vector(x[["sub"]]) )
vars <-  sweep(vars,1,rowMeans(vars),"-")^2 

# different generator matrices give same answers?
stopifnot( all( rowMeans(vars) < 1e-3 ) )
expect_equal( rowMeans(vars), rep(0.0,nrow(vars)) )

# and these agree with truth pretty well?
expect_true( all( sapply( all.expected, function (x) { ( all( (0.5 - abs( ppois( as.vector(subcounts), lambda=as.vector(x[["subexpected"]]) ) - 0.5 )) > (1/2)*(1-(1-.001)^(1/length(as.vector(subcounts)))) ) ) } ) ) )


