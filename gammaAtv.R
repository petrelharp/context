# Let X be a Markov chain with generator matrix A
#  and T ~ Gamma(shape=alpha,scale=beta)
#  so that
#    P(X_T=y|X_0=x) = \int_0^\infty ( \beta^\alpha / \Gamma(\alpha) ) t^{\alpha-1} \exp( - \beta t + A t ) dt
#  which if (\lambda,v) is the eigendecomposition of A is
#                   = sum_i v_i v_i^T \beta^\alpha / (\beta - \lambda_i)^\alpha
# 
# Note that if N(t) ~ Pois(lambda t)
#  and T ~ Gamma(shape=alpha,scale=beta)
#  then N(T) ~ NegativeBinomial(alpha,lambda/(lambda+beta))
#  i.e. P(N(T)=n) = Gamma(n+alpha)/(Gamma(alpha)*factorial(n)) * (beta/(beta+lambda))^alpha * (lambda/(beta+lambda))^n
#  and that
#      P(N(T)=n) = (lambda * (n+alpha-1)) / ( n * (beta+lambda) ) * P(N(T)=n-1) .
#
require(Matrix)

gammaAtv <- function (A,scale,shape,v,tol=1e-6) {
    # evaluate transition matrix after gamma-distributed time multiplied by the matrix v
    # here shape=alpha, scale=beta, scale.t=lambda
    totalrates <- rowSums(A)-diag(A)
    scale.t <- max(totalrates)
    stopifnot(scale.t>0)
    P <- (1/scale.t) * A
    diag(P) <- (1-totalrates/scale.t)
    pp <- scale.t/(scale.t+scale)
    ans <- kterm <- (1-pp)^shape * v
    m <- qnbinom(p=1-tol,size=shape,prob=1-pp)
    for (n in 1:m) {
        kterm <- pp * (n+shape-1)/n * ( P %*% kterm )
        ans <- ans + kterm
    }
    return(ans)
}

gammam <- function (A,scale,shape,...) {
    # return whole transition matrix
    gammaAtv(A=A,scale=scale,shape=shape,v=Diagonal(n=nrow(A)))
}

expAtv.poisson <- function (A,t,v,tol=1e-6) {
    # use a similar strategy to evaluate exp(Atv) for generator matrices
    #   for comparison:
    # number of jumps after time t is Pois( t * scale.t )
    totalrates <- rowSums(A)-diag(A)
    scale.t <- max(totalrates)
    stopifnot(scale.t>0)
    P <- (1/scale.t) * A
    diag(P) <- (1-totalrates/scale.t)
    ans <- kterm <- exp(-scale.t*t) * v
    m <- qpois(p=1-tol,lambda=scale.t*t)
    for (n in 1:m) {
        kterm <- ( t * scale.t / n ) * ( P %*% kterm )
        ans <- ans + kterm
    }
    return(ans)
}

## test
if (FALSE) {

    M <- matrix(rexp(100),nrow=10,ncol=10)
    M <- (M + t(M))/2
    diag(M) <- (-1)*( rowSums(M) - diag(M) )
    eM <- eigen(M)
    stopifnot(all(eM$values <= 0))

    scale <- 2
    shape <- 3
    scale.t <- max((-1)*diag(M))
    gM <- gammam(M,scale=scale,shape=shape)

    true.gM <- matrix(0,nrow=nrow(M),ncol=ncol(M))
    for (k in 1:nrow(M)) { true.gM <- true.gM + outer(eM$vectors[,k],eM$vectors[,k]) * scale^shape / (scale-eM$values[k])^shape }
    stopifnot( all( abs(rowSums(true.gM)) < 1e-8 ) )

    stopifnot( all( abs(true.gM - gM) < 1e-6 ) )

    lcoefs <- shape*log(scale) + (0:1000)*log(scale.t) + lgamma((0:1000)+shape) - lfactorial(0:1000) - lgamma(shape) - ((0:1000)+shape)*log(scale+scale.t)
    stopifnot( any( (exp(lcoefs) - dnbinom(x=0:1000,size=shape,prob=1-scale.t/(scale.t+scale))) < 1e-8 ) )
    stopifnot( abs(sum(exp(lcoefs)) - 1) < 1e-8 )

    # compare expAtv to expAtv.poisson
    source("codon-inference-fns.R")
    load("bcells/genmatrices/genmatrix-6-none-0-2.RData")
    projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), lwin=2, rwin=2 )
    for (tlen in c(.1,1,10)) {
        show(system.time( {
                totalrates <- rowSums(genmatrix)
                scale.t <- mean(totalrates)
                A <- (1/scale.t) * ( genmatrix - Diagonal(nrow(genmatrix),totalrates) )
                a <- sapply( 1:ncol(projmatrix), function (k) { expAtv( A=A, t=tlen*scale.t, v=projmatrix[,k] )$eAtv} )
        } ))
        show(system.time( b <- expAtv.poisson( A=genmatrix, t=tlen, v=projmatrix ) ))
        stopifnot(all(abs(a-b)<1e-7))
    }
    # slower!! and scales much worse.  better to figure out how to do gamma with krylov...

}
