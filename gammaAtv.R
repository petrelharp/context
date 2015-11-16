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
# Furthermore, 
#  if  T ~ Exponential(rate=beta),
#  we've already sampled N(T),
#  and want to get a sample of N(S), 
#  where S ~ Exponential(rate=beta+u)
#  then we can let
#    S = T + sum_{k=0}^M T_k
#  where M is Geometric(beta/(beta+u)) and T_k are independent copies of T.
# This implies that we can write
#    N(S) = N(T) + N' 
#  where N' ~ NegBinom(M,lambda/(lambda+beta))
#  which is the same as
#        N' ~ NegBinom(M,lambda/(lambda+beta+u))

require(Matrix)

gammaAtv <- function (A,scale,shape,v,tol=1e-6,verbose=FALSE) {
    # evaluate transition matrix after gamma-distributed time,
    #    multiplied by the matrix v
    # here shape=alpha, scale=1/beta, scale.t=lambda
    #
    # Note that diag(A) is not used.
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


extend.expAtv <- function (A,scale,dscale,gAtv,tol=1e-8,verbose=FALSE,shape) {
    # Suppose that gAtv is the output of gammaAtv(A,shape=1,scale=scale,v)
    # and we want to return gammaAtv(A,shape=1,scale=scale+dscale,v)
    #
    # The number of extra events that occur is zero with probability scale/(scale+dscale)
    # and otherwise is a Geometric with parameter (scale+dscale)*lambda/(1+lambda*(scale+dscale)),
    # (which can still be zero);
    # in total the probability of being zero is
    #   P(N'=0) = scale/(scale+dscale) + dscale/(scale+dscale) * 1/(1+lambda*(scale+dscale))
    # and 
    #   P(N'=1) = dscale/(scale+dscale) * 1/(1+lambda*(scale+dscale)) * lambda*(scale+dscale)/(1+lambda*(scale+dscale))
    #   P(N'=n+1) = P(N'=n) * lambda*(scale+dscale)/(1+lambda*(scale+dscale))  for n>1
    #
    # Below, scale.t is lambda.
    # 
    # Note that diag(A) is not used.
    if (!missing(shape)) { stop("not implemented for shape != 1: need to extend to a binomial number of these.") }
    totalrates <- rowSums(A)-diag(A)
    scale.t <- max(totalrates)
    stopifnot(scale.t>0)
    P <- (1/scale.t) * A
    diag(P) <- (1-totalrates/scale.t)
    pp <- scale.t*(scale+dscale)/(1+scale.t*(scale+dscale))
    bg.prob <- dscale/(scale+dscale) * (1-pp)
    kterm <- bg.prob * gAtv
    total.prob <- (scale/(scale+dscale) + bg.prob)
    ans <- (scale/(scale+dscale) + bg.prob) * gAtv
    m <- qgeom(p=1-tol,prob=1-pp)
    if (verbose) { cat("gammaAtv: summing ", m, "terms\n") }
    for (n in 1:m) {   # N(T)=n
        bg.prob <- bg.prob * pp
        total.prob <- total.prob + bg.prob  # total prob so far
        kterm <- pp * ( P %*% kterm )
        ans <- ans + kterm
    }
    stopifnot((1-total.prob)<tol)
    return(ans/total.prob)  # normalize by remaining probability so it sums to 1
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

    # Claim: if
    #   N ~ NegBinom(M,q)  
    #    (in R's notation, so P(N=x) = Gamma(x+M)/(Gamma(M) x!) q^M (1-q)^x )
    #   M ~ Geo(p)
    #    (also in R's notation, so P(M=x) = p (1-p)^x )
    # then P(N=0) = 1-(1-p)(1-q)/(1-(1-p)q)
    # and  P(N=n) = p((1-q)/(1-(1-p)q))^n
    rneggeo <- function (n,p,q) {
        M <- rgeom(n,p)
        ans <- numeric(n)
        ans[M>0] <- rnbinom(sum(M>0),size=M[M>0],prob=rep_len(q,n)[M>0])
        ans
    }
    dneggeo <- function (x,p,q) {
        u <- rep_len((1-p)*q,length(x))
        ans <- p*u/(1-u)*((1-q)/(1-u))^x
        ans[x==0] <- ans[x==0]/u[x==0]
        ans
    }
    p <- 0.2
    q <- 0.3
    xx <- rneggeo(1e6,p,q)
    comp <- data.frame( n=seq(0,max(xx)), obs=tabulate(1+xx)/length(xx), exp=dneggeo(0:max(xx),p,q) )
    with( comp, plot( exp ~ obs, log='xy' ) )
    abline(0,1)
}
