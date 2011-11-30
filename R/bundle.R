bundle <-
function (B, DD, yy, pp, lambda, smooth, nb, center, constmat) 
{
    nterms = length(nb)
    m = length(yy)
    np = length(pp)
    lala <- matrix(lambda, nrow = nterms, ncol = 2, dimnames = list(1:nterms, 
        c("mean", "residual")))
    vector.a.ma.schall <- matrix(NA, nrow = sum(nb) + 1 * center, 
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
            residuals = yy - aa$fitted
            b <- rep(1, ncol(B))
            cc <- pp - 0.5
            if (any(cc != 0)) 
                for (i in 1:20) {
                  mo <- fitampllsfit(residuals, B, b, pp, cc, 
                    DD, lala[, 2], nb)
                  b <- mo$b/((sum(mo$b^2)/length(mo$b))^(1/2))
                  c0 <- cc
                  cc <- fitasy(residuals, B, b, pp, cc)
                  dc <- max(abs(cc - c0))
                  if (dc < 1e-06) 
                    break
                }
            for (q in 1:np) {
                vector.a.ma.schall[, q] = mean.coefficients + 
                  cc[q] * b
            }
            sig.med <- vector()
            tau.med <- vector()
            sig.res <- vector()
            tau.res <- vector()
            lmed <- lala[, 1]
            lres <- lala[, 2]
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
                v <- partDD %*% partaa
                z <- aa$fitted
                w[, i] <- 0.5 * (yy > z) + 0.5 * (yy <= z)
                H = solve(t(partB) %*% (w[, i] * partB) + lala[i, 
                  1] * t(partDD) %*% partDD)
                H = apply(sqrt(w[, i]) * partB, 1, function(x) {
                  t(x) %*% H %*% x
                })
                sig.med[i] <- sum(w[, i] * (yy - z)^2, na.rm = TRUE)/(m - 
                  sum(aa$diag.hat.ma))
                tau.med[i] <- sum(v^2)/sum(diag(H)) + 1e-06
                lala[i, 1] <- max(sig.med[i]/tau.med[i], 1e-10)
                partb = b[-1][partbasis]
                v = partDD %*% partb
                z = B %*% b
                w[, i] <- 0.5 * (residuals > z) + 0.5 * (residuals <= 
                  z)
                H = solve(t(partB) %*% (w[, i] * partB) + lala[i, 
                  2] * t(partDD) %*% partDD)
                H = apply(sqrt(w[, i]) * partB, 1, function(x) {
                  t(x) %*% H %*% x
                })
                sig.res[i] <- sum(w[, i] * (residuals - z)^2, 
                  na.rm = TRUE)/(m - sum(aa$diag.hat.ma))
                tau.res[i] <- sum(v^2)/sum(diag(H)) + 1e-06
                lala[i, 2] <- max(sig.res[i]/tau.res[i], 1e-10)
            }
            dc <- max(abs(log10(lmed) - log10(lala[, 1]))) + 
                max(abs(log10(lres) - log10(lala[, 2])))
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
        residuals = yy - B %*% mean.coefficients
        constmat[, ] = 0
        acv.min = nlm(acv, p = lala[, 2], yy = residuals, B = B, 
            quantile = 0.5, DD = DD, nb = nb, constmat = constmat, 
            ndigit = 8, iterlim = 50, gradtol = 1e-04)
        lala[, 2] <- abs(acv.min$estimate)
        b <- rep(1, ncol(B))
        cc <- pp - 0.5
        if (any(cc != 0)) 
            for (i in 1:20) {
                mo <- fitampllsfit(residuals, B, b, pp, cc, DD, 
                  abs(acv.min$estimate), nb)
                b <- mo$b/((sum(mo$b^2)/length(mo$b))^(1/2))
                c0 <- cc
                cc <- fitasy(residuals, B, b, pp, cc)
                dc <- max(abs(cc - c0))
                if (dc < 1e-06) 
                  break
            }
        for (q in 1:np) {
            vector.a.ma.schall[, q] = mean.coefficients + cc[q] * 
                b
        }
    }
    else {
        aa <- asyregpen.lsfit(yy, B, 0.5, lala[, 1], DD, nb, 
            constmat)
        mean.coefficients <- aa$a
        residuals = yy - B %*% mean.coefficients
        b <- rep(1, ncol(B))
        cc <- pp - 0.5
        if (any(cc != 0)) 
            for (i in 1:20) {
                mo <- fitampllsfit(residuals, B, b, pp, cc, DD, 
                  lala[, 2], nb)
                b <- mo$b/((sum(mo$b^2)/length(mo$b))^(1/2))
                c0 <- cc
                cc <- fitasy(residuals, B, b, pp, cc)
                dc <- max(abs(cc - c0))
                if (dc < 1e-06) 
                  break
            }
        for (q in 1:np) {
            vector.a.ma.schall[, q] = mean.coefficients + cc[q] * 
                b
        }
    }
    return(list(vector.a.ma.schall, lala, mean.coefficients, 
        b, cc))
}
