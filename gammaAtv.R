# compute gamma(shape,tA) and gamma(shape,tA)v for a generator matrix A and a vector v
#
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
#                 = 1/Beta(n,alpha) * (beta^alpha * lambda^n)/(beta+lambda)^{n+alpha} .
#  and that
#      P(N(T)=n) = (lambda * (n+alpha-1)) / ( n * (beta+lambda) ) * P(N(T)=n-1) .
#


gammam <- function (A,t,shape,m=10) {
    totalrates <- rowSums(A)-diag(A)
    scale.t <- max(totalrates)
    P <- (1/scale.t) * A
    diag(P) <- (1-totalrates/scale.t)
    ans <- kterm <- Diagonal(nrow(A))
    for (n in 1:m) {
        kterm <- scale.t * (n+shape-1) / ( n * (t+scale.t) ) * ( P %*% kterm )
        ans <- ans + kterm
    }
    return(ans)
}

gammaAtv <- function (A,t,v,shape,m=10) {
    totalrates <- rowSums(A)
    scale.t <- max(totalrates)
    P <- (1/scale.t) * A
    diag(P) <- (1-totalrates/scale.t)
    ans <- kterm <- v
    for (n in 1:m) {
        kterm <- scale.t * (n+shape-1) / ( n * (t+scale.t) ) * ( P %*% kterm )
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

    t <- 2
    shape <- 3
    scale.t <- max((-1)*diag(M))
    gM <- gammam(M,t=t,shape=shape)

    true.gM <- matrix(0,nrow=nrow(M),ncol=ncol(M))
    for (k in 1:nrow(M)) { true.gM <- true.gM + outer(eM$vectors[,k],eM$vectors[,k]) * t^shape / (t-eM$values[k])^shape }
    range(rowSums(true.gM))

    dnbinom(x=0:10,size=.5,prob=4/(4+3))
    4^.5 * 3^(0:10) / ( beta((0:10),4) * (4+3)^(.5+(0:10)) )

    lcoefs <- shape*log(t) + (0:1000)*log(scale.t) - lbeta(0:1000,shape) + lgamma(0:1000+shape) - (0:1000+shape)*log(t+scale.t)
    range(exp(lcoefs) - dnbinom(x=0:1000,size=shape,prob=scale.t/(scale.t+t)))
    sum(exp(lcoefs))

}
