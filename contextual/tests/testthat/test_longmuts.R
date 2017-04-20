library(testthat)

library(contextual)

bases <- c("O","X")

mutpats <- list( 
    list( c("O","X") ),
    list( c("OX","OO"), c("XO","OO") )
    ) 
mutrates <- c(3,5)
selpats <- list(
        c( "X" )
    )
selcoef <- c(1)

fixfn <- function (ds,...) { 1/(1+exp(-ds)) }


###### 
# Can make correct generator matrices?
patlen <- 3

pats <- c("OOO", "XOO", "OXO", "XXO", "OOX", "XOX", "OXX", "XXX")
nA <- rbind(   # number of O->X:
        #   OOO XOO OXO XXO OOX XOX OXX XXX
        c(   0,  1,  1,  0,  1,  0,  0,  0 ), # OOO
        c(  -1,  0,  0,  1,  0,  1,  0,  0 ), # XOO
        c(  -1,  0,  0,  1,  0,  0,  1,  0 ), # OXO
        c(   0, -1, -1,  0,  0,  0,  0,  1 ), # XXO
        c(  -1,  0,  0,  0,  0,  1,  1,  0 ), # OOX
        c(   0, -1,  0,  0, -1,  0,  0,  1 ), # XOX
        c(   0,  0, -1,  0, -1,  0,  0,  1 ), # OXX
        c(   0,  0,  0, -1,  0, -1, -1,  0 )  # XXX
    )
nB <- rbind(   # number of OX->OO:
        #   OOO XOO OXO XXO OOX XOX OXX XXX
        c(   0,  0,  0,  0,  0,  0,  0,  0 ), # OOO
        c(   0,  0,  0,  0,  0,  0,  0,  0 ), # XOO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # OXO
        c(   0,  0,  0,  0,  0,  0,  0,  0 ), # XXO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # OOX
        c(   0,  1,  0,  0,  0,  0,  0,  0 ), # XOX
        c(   0,  0,  0,  0,  1,  0,  0,  0 ), # OXX
        c(   0,  0,  0,  0,  0,  0,  0,  0 )  # XXX
    )
nC <- rbind(   # number of XO->OO:
        #   OOO XOO OXO XXO OOX XOX OXX XXX
        c(   0,  0,  0,  0,  0,  0,  0,  0 ), # OOO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # XOO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # OXO
        c(   0,  1,  0,  0,  0,  0,  0,  0 ), # XXO
        c(   0,  0,  0,  0,  0,  0,  0,  0 ), # OOX
        c(   0,  0,  0,  0,  1,  0,  0,  0 ), # XOX
        c(   0,  0,  0,  0,  0,  0,  0,  0 ), # OXX
        c(   0,  0,  0,  0,  0,  0,  0,  0 )  # XXX
    )
dimnames(nA) <- dimnames(nB) <- dimnames(nC) <- list(pats,pats)

checkit <- function (x,y) { 
    expect_equal( unlist(dimnames(x)), unlist(dimnames(y)))
    expect_equal(dim(x),dim(y))
    expect_equal(as.vector(x),as.vector(y))
}

checkit(
        makegenmatrix( mutpats=mutpats, selpats=selpats, selcoef=selcoef, mutrates=mutrates, patlen=patlen, boundary="none", bases=bases, fixfn=null.fixfn ),
        3 * (nA>0) + 5 * ( nB + nC )
    )


## wrapped

nA.wrap <- nA
nB.wrap <- rbind(   # number of OX->OO:
        #   OOO XOO OXO XXO OOX XOX OXX XXX
        c(   0,  0,  0,  0,  0,  0,  0,  0 ), # OOO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # XOO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # OXO
        c(   0,  0,  1,  0,  0,  0,  0,  0 ), # XXO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # OOX
        c(   0,  1,  0,  0,  0,  0,  0,  0 ), # XOX
        c(   0,  0,  0,  0,  1,  0,  0,  0 ), # OXX
        c(   0,  0,  0,  0,  0,  0,  0,  0 )  # XXX
    )
nC.wrap <- rbind(   # number of XO->OO:
        #   OOO XOO OXO XXO OOX XOX OXX XXX
        c(   0,  0,  0,  0,  0,  0,  0,  0 ), # OOO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # XOO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # OXO
        c(   0,  1,  0,  0,  0,  0,  0,  0 ), # XXO
        c(   1,  0,  0,  0,  0,  0,  0,  0 ), # OOX
        c(   0,  0,  0,  0,  1,  0,  0,  0 ), # XOX
        c(   0,  0,  1,  0,  0,  0,  0,  0 ), # OXX
        c(   0,  0,  0,  0,  0,  0,  0,  0 )  # XXX
    )
dimnames(nA.wrap) <- dimnames(nB.wrap) <- dimnames(nC.wrap) <- list(pats,pats)

checkit(
        makegenmatrix( mutpats=mutpats, selpats=selpats, selcoef=selcoef, mutrates=mutrates, patlen=patlen, boundary="wrap", bases=bases, fixfn=ising.fixfn ),
        ( 3 * (nA.wrap>0) + 5 * ( nB.wrap + nC.wrap ) ) * ising.fixfn( 1 * nA )
    )
