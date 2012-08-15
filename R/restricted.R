restricted <-
function (B, DD, yy, pp, lambda, smooth, nb, center, constmat) 
{
    nterms = length(nb)
    m = length(yy)
    np = length(pp)
    lala <- matrix(lambda, nrow = nterms, ncol = 2, dimnames = list(1:nterms, 
        c("mean", "residual")))
    vector.a.ma.schall <- matrix(NA, nrow = sum(nb) + (1 * center), 
        ncol = np)
    if (smooth == "schall") {
        dc = 1
        dw = 1
        w <- matrix(0.5, nrow = m, ncol = nterms)
        it = 1
        while (dc >= 0.01 && it < 100) {
            aa <- asyregpen.lsfit(yy, B, 0.5, lala[, 1], DD, 
                nb, constmat)
            mean.coefficients <- aa$a
            diag.hat = aa$diag.hat.ma
            sig2 <- vector()
            tau2 <- vector()
            l0 <- lala[, 1]
            for (i in 1:nterms) {
                partbasis = (sum(nb[0:(i - 1)]) + 1):(sum(nb[0:i]))
                if (center) {
                  partB = B[, -1, drop = FALSE][, partbasis, 
                    drop = FALSE]
                  partDD = DD[, -1, drop = FALSE][-1, , drop = FALSE][, 
                    partbasis, drop = FALSE]
                  partaa = aa$a[-1][partbasis]
                }
                else {
                  partB = B[, partbasis, drop = FALSE]
                  partDD = DD[, partbasis, drop = FALSE]
                  partaa = aa$a[partbasis]
                }
                v <- partDD %*% partaa
                z <- aa$fitted
                H = solve(t(partB) %*% (w[, i] * partB) + lala[i, 
                  1] * t(partDD) %*% partDD)
                H = apply(sqrt(w[, i]) * partB, 1, function(x) {
                  t(x) %*% H %*% x
                })
                sig2[i] <- sum(w[, i] * (yy - z)^2, na.rm = TRUE)/(m - 
                  sum(aa$diag.hat.ma, na.rm = TRUE))
                tau2[i] <- sum(v^2, na.rm = TRUE)/sum(H) + 1e-06
                lala[i, 1] <- max(sig2[i]/tau2[i], 1e-10, na.rm = TRUE)
            }
            dc <- max(abs(log10(l0 + 1e-06)) - log10(lala[, 1] + 
                1e-06))
            it = it + 1
        }
        if (it == 100) 
            warning("Schall algorithm did not converge. Stopping after 100 iterations.")
    }
    else if (smooth == "acv") {
        acv.min = nlm(acv, p = lala[, 1], yy = yy, B = B, quantile = 0.5, 
            DD = DD, nb = nb, constmat = constmat, ndigit = 8, 
            iterlim = 50, gradtol = 1e-04)
        aa <- asyregpen.lsfit(yy, B, 0.5, abs(acv.min$estimate), 
            DD, nb, constmat)
        mean.coefficients <- aa$a
        lala[, 1] <- abs(acv.min$estimate)
        diag.hat = aa$diag.hat.ma
    }
    else {
        aa <- asyregpen.lsfit(yy, B, 0.5, lala[, 1], DD, nb, 
            constmat)
        mean.coefficients <- aa$a
        diag.hat = aa$diag.hat.ma
    }
    residuals = yy - B %*% mean.coefficients
    constmat[, ] = 0
    gg = asyregpen.lsfit(abs(residuals), B, 0.5, lala[, 1], DD, 
        nb, constmat)
    cc = NULL
    for (q in 1:np) {
        ca = asyregpen.lsfit(residuals, gg$fitted, pp[q], NULL, 
            matrix(0, nrow = 1, ncol = 1), 0, matrix(0, nrow = 1, 
                ncol = 1))
        cc[q] = ca$a
        vector.a.ma.schall[, q] = mean.coefficients + ca$a * 
            gg$a
    }
    diag.hat = matrix(diag.hat, nrow = length(diag.hat), ncol = np)
    return(list(vector.a.ma.schall, lala, diag.hat, mean.coefficients, 
        gg$a, cc))
}
