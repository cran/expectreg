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
        it = 1
        while (dc >= 0.01 && it < 100) {
            aa <- asyregpen.lsfit(yy, B, 0.5, lala[, 1], DD, 
                nb, constmat)
            mean.coefficients <- aa$a
            diag.hat = aa$diag.hat.ma
            sig.med <- vector()
            tau.med <- vector()
            lmed <- lala[, 1]
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
                H = solve(t(partB) %*% (partB) + lala[i, 1] * 
                  t(partDD) %*% partDD)
                H = apply(partB, 1, function(x) {
                  t(x) %*% H %*% x
                })
                sig.med[i] <- sum(0.5 * (yy - z)^2, na.rm = TRUE)/(m - 
                  sum(aa$diag.hat.ma, na.rm = TRUE))
                tau.med[i] <- sum(v^2, na.rm = TRUE)/sum(H, na.rm = TRUE) + 
                  1e-06
                lala[i, 1] <- max(sig.med[i]/tau.med[i], 1e-10, 
                  na.rm = TRUE)
            }
            dc <- max(abs(log10(lmed + 1e-06) - log10(lala[, 
                1] + 1e-06)))
            it = it + 1
        }
        if (it == 100) 
            warning("Schall algorithm did not converge. Stopping after 100 iterations.")
        residuals = yy - aa$fitted
        dc = 1
        it = 1
        while (dc >= 0.01 && it < 100) {
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
            sig.res <- vector()
            tau.res <- vector()
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
                partb = b[-1][partbasis]
                v = partDD %*% partb
                z = B %*% b
                H = solve(t(partB) %*% (partB) + lala[i, 2] * 
                  t(partDD) %*% partDD)
                H = apply(partB, 1, function(x) {
                  t(x) %*% H %*% x
                })
                sig.res[i] <- sum(0.5 * (residuals - z)^2, na.rm = TRUE)/(m - 
                  sum(aa$diag.hat.ma, na.rm = TRUE))
                tau.res[i] <- sum(v^2, na.rm = TRUE)/sum(H, na.rm = TRUE) + 
                  1e-06
                lala[i, 2] <- max(sig.res[i]/tau.res[i], 1e-10, 
                  na.rm = TRUE)
            }
            dc <- max(abs(log10(lres + 1e-06) - log10(lala[, 
                2] + 1e-06)))
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
        diag.hat = aa$diag.hat.ma
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
    diag.hat = matrix(diag.hat, nrow = length(diag.hat), ncol = np)
    return(list(vector.a.ma.schall, lala, diag.hat, mean.coefficients, 
        b, cc))
}
