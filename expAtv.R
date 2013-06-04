# fixes issues in expm:::expAtv

expAtv <- function (A, v, t = 1, method = "Sidje98", tol = 1e-07, btol = 1e-07, 
    m.max = 30, mxrej = 10, verbose = getOption("verbose")) 
{
    stopifnot(length(v) == (n <- nrow(A)), m.max >= 2)
    if (n <= 1) {
        if (n == 1) 
            return(list(eAtv = exp(A * t) * v, error = 0, nstep = 0L, 
                n.reject = 0L))
        stop("nrow(A) must be >= 1")
    }
    method <- match.arg(method)
    m <- min(n, m.max)
    gamma <- 0.9
    delta <- 1.2
    mb <- m
    t_1 <- abs(t)
    sgn <- sign(t)
    t_now <- 0
    s_error <- 0
    nA <- norm(A, "I")
    rndoff <- nA * .Machine$double.eps
    k1 <- 2
    xm <- 1/m
    beta <- sqrt(sum(v * v))
    fact <- (((m + 1)/exp(1))^(m + 1)) * sqrt(2 * pi * (m + 1))
    myRound <- function(tt) {
        s <- 10^(floor(log10(tt)) - 1)
        ceiling(tt/s) * s
    }
    t_new <- myRound((1/nA) * (fact * tol/(4 * beta * nA))^xm)
    V <- matrix(0, n, m + 1)
    H <- matrix(0, m + 2, m + 2)
    nstep <- n.rej <- 0L
    w <- v
    while (t_now < t_1) {
        nstep <- nstep + 1L
        t_step <- min(t_1 - t_now, t_new)
        if (verbose) 
            cat(sprintf("while(t_now = %g < ..): nstep=%d, t_step=%g\n", 
                t_now, nstep, t_step))
        V[, 1] <- (1/beta) * w
        for (j in 1:m) {
            p <- as.vector(A %*% V[, j])
            for (i in 1:j) {
                H[i, j] <- s <- sum(V[, i] * p)
                p <- p - s * V[, i]
            }
            s <- sqrt(sum(p * p))
            if (s < btol) {
                k1 <- 0
                mb <- j
                t_step <- t_1 - t_now
                break
            }
            H[j + 1, j] <- s
            V[, j + 1] <- p/s
        }
        if (k1 != 0) {
            H[m + 2, m + 1] <- 1
            av <- A %*% V[, m + 1]
            avnorm <- sqrt(sum(av * av))
        }
        i.rej <- 0L
        while (i.rej <= mxrej) {
            mx <- mb + k1
            imx <- seq_len(mx)
            if (verbose) 
                cat(sprintf("\tinner while: k1=%d -> mx=%d\n", k1, mx))
            F <- expm(sgn * t_step * H[imx, imx, drop=FALSE])
            if (k1 == 0) {
                err_loc <- btol
                break
            }
            else {
                phi1 <- abs(beta * F[m + 1, 1])
                phi2 <- abs(beta * F[m + 2, 1] * avnorm)
                if (phi1 > 10 * phi2) {
                  err_loc <- phi2
                  xm <- 1/m
                }
                else if (phi1 > phi2) {
                  err_loc <- (phi1 * phi2)/(phi1 - phi2)
                  xm <- 1/m
                }
                else {
                  err_loc <- phi1
                  xm <- 1/(m - 1)
                }
            }
            if (err_loc <= delta * t_step * tol) 
                break
            else {
                if (i.rej == mxrej) 
                  stop(gettextf("The requested tolerance (tol=%g) is too small for mxrej=%d.", 
                    tol, mxrej))
                t_step <- gamma * t_step * (t_step * tol/err_loc)^xm
                s <- 10^(floor(log10(t_step)) - 1)
                t_step <- s * ceiling(t_step/s)
                i.rej <- i.rej + 1L
            }
        }
        n.rej <- n.rej + i.rej
        mx <- mb + max(0, k1 - 1)
        imx <- seq_len(mx)
        w <- as.vector(V[, imx, drop=FALSE] %*% (beta * F[imx, 1]))
        beta <- sqrt(sum(w * w))
        t_now <- t_now + t_step
        t_new <- myRound(gamma * t_step * (t_step * tol/err_loc)^xm)
        err_loc <- max(err_loc, rndoff)
        s_error <- s_error + err_loc
    }
    list(eAtv = w, error = s_error, nstep = nstep, n.reject = n.rej)
}
