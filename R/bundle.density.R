bundle.density <-
function (bundle) 
{
    if (!inherits(bundle, "bundle") && !inherits(bundle, "restricted")) 
        stop("Function needs output from 'expectile.bundle()' or 'expectile.restricted()'.")
    basis = bundle$design
    np = length(bundle$intercepts)
    pp <- bundle$expectiles
    rst <- (bundle$response - (basis %*% bundle$trend.coef))/(basis %*% 
        bundle$residual.coef)
    u <- seq(1.2 * min(rst), 1.2 * max(rst), length = 100)
    m <- length(u)
    bundle$asymmetry = sort(bundle$asymmetry)
    A <- matrix(0, np + 1, m)
    A[np + 1, ] <- 1
    for (k in 1:np) {
        a1 <- (1 - pp[k]) * (u - bundle$asymmetry[k]) * (u <= 
            bundle$asymmetry[k])
        a2 <- pp[k] * (u - bundle$asymmetry[k]) * (u > bundle$asymmetry[k])
        A[k, ] <- a1 + a2
    }
    D <- diff(diag(m), diff = 2)
    lambda <- 100
    P <- lambda * t(D) %*% D
    v1 <- solve(t(A) %*% A + P, t(A) %*% c(rep(0, np), 1))
    q <- c(rep(0, np), 1)
    lambda2 <- 1
    D2 <- diff(diag(m), diff = 3)
    P2 <- lambda2 * t(D2) %*% D2
    z <- log(v1 - min(v1) + 0.02 * max(v1))
    for (it in 1:20) {
        g <- exp(z)
        r <- q - A %*% g
        B <- A * outer(rep(1, np + 1), as.vector(g))
        Q <- t(B) %*% B
        znew <- solve(Q + P2, t(B) %*% r + Q %*% z)
        dz <- max(abs(z - znew))
        z <- znew
        cat("iteration: ", it, ", convergence: ", dz, "\n")
        if (dz < 1e-06) 
            break
    }
    dev.new()
    hist(rst, breaks = seq(min(rst), max(rst), length = 30), 
        freq = F, xlim = range(u), main = "density")
    lines(u, g/(u[2] - u[1]), col = "red")
    result = list(random = rst, density = g/(u[2] - u[1]), x = u)
    class(result) = "bundledensity"
    result
}
