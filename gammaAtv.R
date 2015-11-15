# Let X be a Markov chain with generator matrix A
#  and T ~ Gamma(shape=alpha,rate=beta)
#  so that
#    P(X_T=y|X_0=x) = \int_0^\infty ( \beta^\alpha / \Gamma(\alpha) ) t^{\alpha-1} \exp( - \beta t + A t ) dt
#  which if (\lambda,v) is the eigendecomposition of A is
#                   = sum_i v_i v_i^T \beta^\alpha / (\beta - \lambda_i)^\alpha
# 
# Note that if N(t) ~ Pois(lambda t)
#  and T ~ Gamma(shape=alpha,rate=beta)
#  then N(T) ~ NegativeBinomial(alpha,lambda/(lambda+beta))
#  i.e. P(N(T)=n) = Gamma(n+alpha)/(Gamma(alpha)*factorial(n)) * (beta/(beta+lambda))^alpha * (lambda/(beta+lambda))^n
#  and that
#      P(N(T)=n) = (lambda * (n+alpha-1)) / ( n * (beta+lambda) ) * P(N(T)=n-1) .
#
library(Matrix)

gammaAtv <- function (A,scale,shape,v,tol=1e-6,verbose=FALSE) {
    # evaluate transition matrix after gamma-distributed time,
    #    multiplied by the matrix v
    # here shape=alpha, scale=1/beta, scale.t=lambda
    totalrates <- rowSums(A)-diag(A)
    scale.t <- max(totalrates)
    stopifnot(scale.t>0)
    P <- (1/scale.t) * A
    diag(P) <- (1-totalrates/scale.t)
    pp <- scale.t/(scale.t+1/scale)
    nb.prob <- total.prob <- (1-pp)^shape
    ans <- kterm <- nb.prob * v   # N(T)=0
    m <- qnbinom(p=1-tol,size=shape,prob=1-pp)
    if (verbose) { cat("gammaAtv: summing ", m, "terms\n") }
    for (n in 1:m) {   # N(T)=n
        nb.prob <- nb.prob * pp * (n+shape-1)/n  # neg binom prob
        total.prob <- total.prob + nb.prob  # total prob so far
        kterm <- pp * (n+shape-1)/n * ( P %*% kterm )
        ans <- ans + kterm
    }
    return(ans/total.prob)  # normalize by remaining probability so it sums to 1
}

gammam <- function (A,scale,shape,...) {
    # return whole transition matrix
    gammaAtv(A=A,scale=scale,shape=shape,v=Diagonal(n=nrow(A)),...)
}

extend.gammam <- function (P,a,da,eps=da/(a+da),tol=1e-8) {
    # if P is e^{TG} where T is Exp(a)
    #   return e^{SG} where S is Exp(a+da), with da>0,
    #   a+da = a/(1-eps), so eps=(da/a)/(1+da/a)
    ans <- Pm <- (1-eps)*P
    m <- ceiling(log(tol)/log(eps))
    # cat("m=",m,"\n")
    for (k in seq(2,length.out=(m-1))) {
        Pm <- eps * P %*% Pm
        ans <- ans + Pm
    }
    return(ans)
}

expAtv.poisson <- function (A,t,v,tol=1e-8) {
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

    lcoefs <- shape*log(scale) + (0:1000)*log(scale.t) + lgamma((0:1000)+shape) - lfactorial(0:1000) - lgamma(shape) - ((0:1000)+shape)*log(scale+scale.t)
    stopifnot( any( (exp(lcoefs) - dnbinom(x=0:1000,size=shape,prob=1-scale.t/(scale.t+scale))) < 1e-8 ) )
    stopifnot( abs(sum(exp(lcoefs)) - 1) < 1e-8 )

    # compare expAtv to expAtv.poisson
    source("context-inference-fns.R",chdir=TRUE)
    load("bcells/genmatrices/genmatrix-6-none-0-2.RData")
    projmatrix <- collapsepatmatrix( ipatterns=rownames(genmatrix), leftwin=2, rightwin=2 )
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
