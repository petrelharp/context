source("../gammaAtv.R")

# randomly generated matrix
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
for (k in 1:nrow(M)) { true.gM <- true.gM + outer(eM$vectors[,k],eM$vectors[,k]) * (1/scale)^shape / ((1/scale)-eM$values[k])^shape }
stopifnot( all( abs(rowSums(true.gM)-1) < 1e-8 ) )

stopifnot( all( abs(true.gM - gM) < 1e-6 ) )
stopifnot( all( abs(true.gM - gammam(M,scale=scale,shape=shape,tol=1e-12)) < 1e-12 ) )

# construct an answer:
#    gM[1,1] should be the probability that an Exponential(M[1,2]) still hasn't happened by an Exponential(1/scale) time,
#    which is 1/(1+scale*M[1,2])
M <- matrix( c(-1,0,1,0), nrow=2 )
scale <- 2
gM <- gammam(M, scale=scale, shape=1, tol=1e-8)
stopifnot( all( abs( gM - matrix( c(1/(1+scale*M[1,2]), 0, 1-1/(1+scale*M[1,2]), 1.0), nrow=2) ) < 1e-8 ) )
gM <- gammam(M, scale=scale, shape=1, tol=1e-14)
stopifnot( all( abs( gM - matrix( c(1/(1+scale*M[1,2]), 0, 1-1/(1+scale*M[1,2]), 1.0), nrow=2) ) < 1e-14 ) )

# check extend.gammam:
# prob an exponential(theta) is greater than an exponential(1) is 1/(1+theta), so:
stopifnot( abs( 1/(1+0.25/.8) - extend.gammam( 1/(1+0.25), eps=0.2 ) ) <= 1e-8 )
stopifnot( abs( 1/(1+0.25/.8) - extend.gammam( 1/(1+0.25), eps=0.2, tol=1e-15 ) ) <= 1e-15 )

# and, this should work for matrices, too
extend.fac <- 0.8
egM <- extend.gammam( gM, eps=1-extend.fac, tol=1e-8 )
stopifnot( all( abs( egM - matrix( c(1/(1+scale/extend.fac*M[1,2]), 0, 1-1/(1+scale/extend.fac*M[1,2]), 1.0), nrow=2) ) < 1e-8 ) )
stopifnot( all( abs( egM - gammam( M, scale=scale/extend.fac, shape=1, tol=1e-8 ) ) < 2e-8 ) )

# try longer extensions
extend.fac <- 0.5
egM <- extend.gammam( gM, eps=1-extend.fac, tol=1e-8 )
stopifnot( all( abs( egM - matrix( c(1/(1+scale/extend.fac*M[1,2]), 0, 1-1/(1+scale/extend.fac*M[1,2]), 1.0), nrow=2) ) < 1e-8 ) )
stopifnot( all( abs( egM - gammam( M, scale=scale/extend.fac, shape=1, tol=1e-8 ) ) < 2e-8 ) )

# bigger matrix
M <- matrix( c(
    0,1,1,1,1,0,0,0,1,0,0,0,1,0,0,0,#AA
    1,0,1,1,0,1,0,0,0,1,0,0,0,1,0,0,#CA
    1,1,0,1,0,0,1,0,0,0,1,0,0,0,1,0,#GA
    1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,#TA
    1,0,0,0,0,1,1,1,1,0,0,0,1,0,0,0,#AC
    0,1,0,0,1,0,1,1,0,1,0,0,0,1,0,0,#CC
    0,0,1,0,1,1,0,1,0,0,1,0,0,0,1,0,#GC
    0,0,0,1,1,1,1,0,0,0,0,1,0,0,0,1,#TC
    1,0,0,0,1,0,0,0,0,1,1,1,1,0,0,0,#AG
    0,1,0,0,0,1,0,0,1,0,1,1,0,1,0,0,#CG
    0,0,1,0,0,0,1,0,1,1,0,1,0,0,1,0,#GG
    0,0,0,1,0,0,0,1,1,1,1,0,0,0,0,1,#TG
    1,0,0,0,1,0,0,0,1,0,0,0,0,1,1,1,#AT
    0,1,0,0,0,1,0,0,0,1,0,0,1,0,1,1,#CT
    0,0,1,0,0,0,1,0,0,0,1,0,1,1,0,1,#GT
    0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,0 #TT
    ), nrow=16, byrow=TRUE )
diag(M) <- (-1)*rowSums(M)

scale <- 0.3
gM <- gammam(M, scale=scale, shape=1, tol=1e-8)
extend.fac <- 0.8
egM <- extend.gammam( gM, eps=1-extend.fac, tol=1e-8 )
stopifnot( all( abs( egM - gammam( M, scale=scale/extend.fac, shape=1, tol=1e-8 ) ) < 2e-8 ) )

###
# test extend.expAtv
scale <- 0.3
dscale <- 0.1
gAtv <- gammaAtv(M, scale=scale, shape=1, v=Diagonal(nrow(M)), tol=1e-8)
extend.fac <- 0.8
egAtv <- extend.expAtv( A=M, scale=scale, dscale=dscale, gAtv=gAtv, tol=1e-8 )
stopifnot( all( abs( egAtv - gammaAtv( M, scale=scale+dscale, shape=1, v=Diagonal(nrow(M)), tol=1e-8 ) ) < 2e-8 ) )

scale <- 2.1
dscale <- 0.05
gAtv <- gammaAtv(M, scale=scale, shape=1, v=Diagonal(nrow(M)), tol=1e-8)
extend.fac <- 0.8
egAtv <- extend.expAtv( A=M, scale=scale, dscale=dscale, gAtv=gAtv, tol=1e-8 )
stopifnot( all( abs( egAtv - gammaAtv( M, scale=scale+dscale, shape=1, v=Diagonal(nrow(M)), tol=1e-8 ) ) < 2e-8 ) )


###
# check extend.expAtv agrees with gammaAtv when scale=0
scale <- 2.1
dscale <- 0.5
gAtv <- gammaAtv(M, scale=scale, shape=1, v=Diagonal(nrow(M)), tol=1e-8)
egAtv <- extend.expAtv(A=M, scale=0, dscale=scale, gAtv=Diagonal(nrow(M)), tol=1e-8)
stopifnot(all(abs(gAtv-egAtv)<.Machine$double.eps))  # should be less than machine tolerance as they do the same thing
