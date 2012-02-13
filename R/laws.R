laws <-
function (B, DD, yy, pp, lambda, smooth, nb, center, constmat) 
{
    nterms = length(nb)
    m = length(yy)
    np = length(pp)
    if (length(lambda) < nterms) 
        lala = rep(lambda[1], nterms)
    else lala = lambda
    dummy.reg <- function(pp, lala, smooth, yy, B, DD, nb, nterms, 
        center) {
        print(paste("Expectile:", pp, sep = " "))
        if (smooth == "schall") {
            dc = 1
            dw = 1
            w <- rep(1, times = m)
            it = 1
            while ((dc >= 0.01 || dw != 0) && it < 100) {
                aa <- asyregpen.lsfit(yy, B, pp, lala, DD, nb, 
                  constmat)
                vector.a.ma.schall <- aa$a
                w0 <- w
                l0 <- lala
                for (i in 1:nterms) {
                  partbasis = (sum(nb[0:(i - 1)]) + 1):(sum(nb[0:i]))
                  if (center) {
                    partB = B[, -1][, partbasis, drop = FALSE]
                    partDD = DD[, -1][-1, ][, partbasis, drop = FALSE]
                    partaa = aa$a[-1][partbasis]
                  }
                  else {
                    partB = B[, partbasis, drop = FALSE]
                    partDD = DD[, partbasis, drop = FALSE]
                    partaa = aa$a[partbasis]
                  }
                  if (any(partDD != 0)) {
                    v <- partDD %*% partaa
                    z <- aa$fitted
                    w <- aa$weight
                    H = solve(t(partB) %*% (w * partB) + lala[i] * 
                      t(partDD) %*% partDD)
                    H = apply(sqrt(w) * partB, 1, function(x) {
                      t(x) %*% H %*% x
                    })
                    sig2 <- sum(w * (yy - z)^2, na.rm = TRUE)/(m - 
                      sum(aa$diag.hat.ma))
                    tau2 <- sum(v^2)/sum(diag(H)) + 1e-06
                    lala[i] = max(sig2/tau2, 1e-10)
                  }
                }
                dc <- max(abs(log10(l0) - log10(lala)))
                dw <- sum(w != w0, na.rm = TRUE)
                it = it + 1
            }
            if (it == 100) 
                warning("Schall algorithm did not converge. Stopping after 100 iterations.")
        }
        else if (smooth == "acv") {
            acv.min = nlm(acv, p = lala, yy = yy, B = B, quantile = pp, 
                DD = DD, nb = nb, constmat = constmat, iterlim = 100)
            aa <- asyregpen.lsfit(yy, B, pp, abs(acv.min$estimate), 
                DD, nb, constmat)
            vector.a.ma.schall <- aa$a
            lala <- abs(acv.min$estimate)
        }
        else {
            aa <- asyregpen.lsfit(yy, B, pp, lala, DD, nb, constmat)
            vector.a.ma.schall <- aa$a
        }
        list(vector.a.ma.schall, lala)
    }
    if (.Platform$OS.type == "unix") 
        coef.vector = mclapply(pp, function(pp) dummy.reg(pp, 
            lala, smooth, yy, B, DD, nb, nterms, center), mc.cores = min(getOption("cores"), 
            4))
    else if (.Platform$OS.type == "windows") 
        coef.vector = mclapply(pp, function(pp) dummy.reg(pp, 
            lala, smooth, yy, B, DD, nb, nterms, center), mc.cores = 1)
    lala <- matrix(lambda, nrow = nterms, ncol = np)
    vector.a.ma.schall <- matrix(NA, nrow = sum(nb) + (1 * center), 
        ncol = np)
    for (i in 1:np) {
        vector.a.ma.schall[, i] = coef.vector[[i]][[1]]
        lala[, i] = coef.vector[[i]][[2]]
    }
    return(list(vector.a.ma.schall, lala))
}
