
#' Evaluate Transition Matrix After a Gamma-Distributed Time
#'
#' Evaluate transition matrix after gamma-distributed time, multiplied by the
#' matrix `v`.
#'
#' The math behind this: let X be a Markov chain with generator matrix A
#'  and T ~ Gamma(shape=alpha,rate=beta)
#'  so that
#'    P(X_T=y|X_0=x) 
#'        = \int_0^\infty ( \beta^\alpha / \Gamma(\alpha) ) t^{\alpha-1} \exp( - \beta t + A t ) dt
#'  which if (\lambda,v) is the eigendecomposition of A is
#'                   = sum_i v_i v_i^T \beta^\alpha / (\beta - \lambda_i)^\alpha
#' 
#' Note that if N(t) ~ Pois(lambda t)
#'  and T ~ Gamma(shape=alpha,rate=beta)
#'  then N(T) ~ NegativeBinomial(alpha,lambda/(lambda+beta))
#'  i.e. P(N(T)=n) = Gamma(n+alpha)/(Gamma(alpha)*factorial(n)) 
#'                   * (beta/(beta+lambda))^alpha * (lambda/(beta+lambda))^n
#'  and that
#'      P(N(T)=n) = (lambda * (n+alpha-1)) / ( n * (beta+lambda) ) * P(N(T)=n-1) .
#'
#' Furthermore, 
#'  if  T ~ Exponential(rate=beta),
#'  we've already sampled N(T),
#'  and want to get a sample of N(S), 
#'  where S ~ Exponential(rate=beta+u)
#'  then we can let
#'    S = T + sum_{k=0}^M T_k
#'  where M is Geometric(beta/(beta+u)) and T_k are independent copies of T.
#' This implies that we can write
#'    N(S) = N(T) + N' 
#'  where N' ~ NegBinom(M,lambda/(lambda+beta))
#'  which is the same as
#'        N' ~ NegBinom(M,lambda/(lambda+beta+u))
#'
#' @param shape Shape parameter for Gamma (=alpha).
#' @param scale Scale parameter (=1/beta). 
#' @param scale.t Overall rate of the transitions (a time scaling).
#'
#' Note that diag(A) is not used.
#'
#' @export
gammaAtv <- function (A,scale,shape,v,tol=1e-6,verbose=FALSE) {
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

#' @describeIn gammaAtv Return the whole transition matrix
#' @export
gammam <- function (A,scale,shape,...) {
    gammaAtv(A=A,scale=scale,shape=shape,v=Diagonal(n=nrow(A)),...)
}

#' @describeIn gammaAtv Extend a transition matrix from gammam() to a longer time.
#' @export
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


#' @describeIn gammaAtv Extend a computed gammaAtv to a longer time.
#' @export
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
    # Note this is a generalization of gammaAtv().
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

