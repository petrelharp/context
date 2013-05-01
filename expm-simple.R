#!/usr/bin/R

expm.simple.poisson <- function (A,n=20) 
{
    # -> only works for stochastic matrices
    dA <- (-1) * diag(A)
    t <- max( dA )
    P <- A/t
    diag(P) <- 1-dA/t
    Pk <- exp(-t) * t * P
    expA <- Pk
    diag(expA) <- exp(-t) + diag(expA)  # identity
    for (k in 2:n) {
        Pk <- (t/k) * P%*%Pk
        expA <- expA + Pk
    }
    return( expA )
}

if (FALSE) {
    # determine number of scaling-squaring steps for the below
    #  since:
    #    ( ... ( expm.simple.poisson( A / 2^k, n ) )^2 ... )^2 takes (n-1) + k matrix multiplications
    #    and want ppois( n, lambda=lambda/2^k,lower.tail=TRUE ) < eps
    eps <- 1e-16
    lambdas <- seq(.1,100,length.out=100)
    optvals <- sapply( lambdas, function (lambda) { ( sapply( 1:20, function (k) { k + sum(ppois(q=1:200,lambda=lambda/2^k,log.p=TRUE,lower.tail=FALSE)>log(eps)) - 1 } ) ) } )
    layout(matrix(1:4,nrow=2))
    matplot( optvals, type='l', ylim=c(0,30), xlab="number of scaling-sqarings", ylab="number of matrix mults" )
    points( apply(optvals,2,which.min), apply(optvals,2,min), pch=20 )
    plot( lambdas, lambdas/2^(apply(optvals,2,which.min)), xlab="lambda", ylab="scale-to-this" )
    plot( lambdas, apply(optvals,2,min), ylab="number of matrix mults" )
    plot( lambdas, apply(optvals,2,which.min), ylab="number of rescalings" )
    lines( lambdas, ceiling( log2(lambdas) - log2(.05) ) )
    # Seems to rescale lambda to be less than .2 (usually, .4) ... (with eps=1e-8); which is enough so that we can take n=4.
    #   Then lambda/2^k < .4   <=>   k = ceiling( log2(lambda/.4) )
    #   Let's go with this. 
    #  Might think that:
    #    need dpois( x=2, lambda=t/2^k, lower.tail=FALSE ) < eps  (see above)
    #      ... which is 1 - exp(-t/2^k) ( 1 + t/2^k + t^2/2^(2k+1) ) 
    #             \approx 1 - (1-t/2^k+t^2/2^(2k+1)-t^3/(3*2^(3k+1)))(1 + t/2^k + t^2/2^(2k+1))
    #             \approx t^3/(3*2^(3k+1))
    #    so get k = (1/3) * log2( t^3/(6*eps) )
    # but something's not quite right...
    lines( lambdas, (1/3) * ( 3 * log2( lambdas ) - log2(6*eps) ) )
}

expm.poisson <- function (A,n=5)
{
    # -> only works for stochastic matrices
    # as above, but determine n automatrically
    #  and do scaling-and-squaring
    # 
    # --> the constants are for eps=1e-8, maybe??  <--
    dA <- (-1) * diag(A)
    dmax <- max( dA )
    P <- A/dmax
    diag(P) <- 1-dA/dmax
    nscales <- max(0, ceiling(log2(dmax) - log2(.2)))  # see above
    t <- dmax/2^(nscales)
    Pk <- exp(-t) * t * P
    expA <- Pk
    diag(expA) <- exp(-t) + diag(expA)  # identity
    for (k in 2:n) {
        Pk <- (t/k) * P%*%Pk
        expA <- expA + Pk
    }
    while (nscales>0) {
        expA <- expA %*% expA
        nscales <- nscales-1
    }
    return( expA )
}
