expectile.bundle <-
function (formula, smooth = c("schall", "acv", "none"), lambda = 0.1, 
    density = FALSE) 
{
    smooth = match.arg(smooth)
    if (density) 
        pp <- seq(0.01, 0.99, by = 0.01)
    else pp <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 
        0.98, 0.99)
    np <- length(pp)
    yy = eval(parse(text = formula[2]))
    m = length(yy)
    design = list()
    x = list()
    types = list()
    bnd = list()
    Zspathelp = list()
    nb = vector()
    krig.phi = list()
    for (i in 1:length(labels(terms(formula)))) {
        design[[i]] = eval(parse(text = labels(terms(formula))[i]))
    }
    nterms = length(design)
    B = design[[1]][[1]]
    DD = design[[1]][[2]]
    x[[1]] = design[[1]][[3]]
    types[[1]] = design[[1]][[4]]
    bnd[[1]] = design[[1]][[5]]
    Zspathelp[[1]] = design[[1]][[6]]
    nb[1] = ncol(design[[1]][[1]])
    krig.phi[[1]] = design[[1]][[7]]
    if (length(labels(terms(formula))) > 1) 
        for (i in 2:length(labels(terms(formula)))) {
            B = cbind(B, design[[i]][[1]])
            DD = rbind(cbind(DD, matrix(0, nrow = dim(DD)[1], 
                ncol = dim(design[[i]][[2]])[1])), cbind(matrix(0, 
                nrow = dim(design[[i]][[2]])[1], ncol = dim(DD)[2]), 
                design[[i]][[2]]))
            x[[i]] = design[[i]][[3]]
            types[[i]] = design[[i]][[4]]
            bnd[[i]] = design[[i]][[5]]
            Zspathelp[[i]] = design[[i]][[6]]
            nb[i] = ncol(design[[i]][[1]])
            krig.phi[[i]] = design[[i]][[7]]
        }
    B = cbind(1, B)
    DD = rbind(0, cbind(0, DD))
    lala <- matrix(lambda, nrow = nterms, ncol = 2, dimnames = list(1:nterms, 
        c("mean", "residual")))
    vector.a.ma.schall <- matrix(NA, nrow = sum(nb) + 1, ncol = np)
    med = which(pp == 0.5)
    penalty <- lambda
    if (smooth == "schall") {
        dc = 1
        dw = 1
        w <- matrix(0.5, nrow = m, ncol = nterms)
        it = 1
        while (dc >= 0.01 && it < 100) {
            aa <- asyregpen.lsfit(yy, B, 0.5, lala[, 1], DD, 
                nb)
            vector.a.ma.schall[, med] <- aa$a
            residuals = yy - B %*% vector.a.ma.schall[, med]
            b <- rep(1, ncol(B))
            cc <- pp - 0.5
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
                vector.a.ma.schall[, q] = vector.a.ma.schall[, 
                  med] + cc[q] * b
            }
            sig.med <- vector()
            tau.med <- vector()
            sig.res <- vector()
            tau.res <- vector()
            lmed <- lala[, 1]
            lres <- lala[, 2]
            for (i in 1:nterms) {
                partbasis = (sum(nb[0:(i - 1)]) + 1):(sum(nb[0:i]))
                partB = B[, -1][, partbasis]
                partDD = DD[, -1][-1, ][partbasis, partbasis]
                partaa = aa$a[-1][partbasis]
                partw = aa$weight[-1][partbasis]
                v <- partDD %*% partaa
                z <- B %*% aa$a
                w[, i] <- 0.5 * (yy > z) + 0.5 * (yy <= z)
                H = solve(t(partB) %*% (w[, i] * partB) + lala[i, 
                  1] * t(partDD) %*% partDD)
                H = apply(sqrt(w[, i]) * partB, 1, function(x) {
                  t(x) %*% H %*% x
                })
                sig.med[i] <- sum(w[, i] * (yy - z)^2, na.rm = TRUE)/(m - 
                  sum(aa$diag.hat.ma))
                tau.med[i] <- sum(v^2)/sum(H) + 1e-06
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
                tau.res[i] <- sum(v^2)/sum(H) + 1e-06
                lala[i, 2] <- max(sig.res[i]/tau.res[i], 1e-10)
            }
            dc <- max(abs(log10(lmed) - log10(lala[, 1]))) + 
                max(abs(log10(lres) - log10(lala[, 2])))
            it = it + 1
        }
    }
    else if (smooth == "acv") {
        acv.min = nlm(acv.bundle, p = lala, yy = yy, B = B, pp = pp, 
            DD = DD, nb = nb, ndigit = 8, iterlim = 50, gradtol = 1e-04)
        min.lambda = matrix(abs(acv.min$estimate), ncol = 2)
        aa <- asyregpen.lsfit(yy, B, 0.5, min.lambda[, 1], DD, 
            nb)
        vector.a.ma.schall[, med] <- aa$a
        lala[, 1] <- min.lambda[, 1]
        lala[, 2] <- min.lambda[, 2]
        residuals = yy - B %*% vector.a.ma.schall[, med]
        b <- rep(1, ncol(B))
        cc <- pp - 0.5
        for (i in 1:20) {
            mo <- fitampllsfit(residuals, B, b, pp, cc, DD, min.lambda[, 
                2], nb)
            b <- mo$b/((sum(mo$b^2)/length(mo$b))^(1/2))
            c0 <- cc
            cc <- fitasy(residuals, B, b, pp, cc)
            dc <- max(abs(cc - c0))
            if (dc < 1e-06) 
                break
        }
        for (q in 1:np) {
            vector.a.ma.schall[, q] = vector.a.ma.schall[, med] + 
                cc[q] * b
        }
    }
    else {
        if (length(penalty) < length(nb)) 
            penalty = rep(penalty[1], length(nb))
        aa <- asyregpen.lsfit(yy, B, 0.5, penalty, DD, nb)
        vector.a.ma.schall[, med] <- aa$a
        residuals = yy - B %*% vector.a.ma.schall[, med]
        b <- rep(1, ncol(B))
        cc <- pp - 0.5
        for (i in 1:20) {
            mo <- fitampllsfit(residuals, B, b, pp, cc, DD, penalty, 
                nb)
            b <- mo$b/((sum(mo$b^2)/length(mo$b))^(1/2))
            c0 <- cc
            cc <- fitasy(residuals, B, b, pp, cc)
            dc <- max(abs(cc - c0))
            if (dc < 1e-06) 
                break
        }
        for (q in 1:np) {
            vector.a.ma.schall[, q] = vector.a.ma.schall[, med] + 
                cc[q] * b
        }
    }
    Z <- list()
    coefficients <- list()
    final.lambdas <- list()
    intercept = vector.a.ma.schall[1, ]
    B = B[, -1]
    vector.a.ma.schall = vector.a.ma.schall[-1, ]
    for (k in 1:length(design)) {
        final.lambdas[[k]] = lala[k, ]
        partbasis = (sum(nb[0:(k - 1)]) + 1):(sum(nb[0:k]))
        dev.new()
        if (types[[k]] == "pspline") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            for (i in 1:np) {
                Z[[k]][, i] <- B[, partbasis] %*% vector.a.ma.schall[partbasis, 
                  i] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i]
            }
            Z[[k]] = Z[[k]][order(x[[k]]), ]
            plot(x[[k]], yy, cex = 0.5, pch = 20, col = "grey42", 
                xlab = "x", ylab = "y", ylim = range(cbind(yy, 
                  Z[[k]])))
            matlines(sort(x[[k]])[seq(1, m, length = min(m, 100))], 
                Z[[k]][seq(1, m, length = min(m, 100)), ], col = rainbow(np + 
                  1)[1:np], lty = 1)
            legend(x = "bottomright", pch = 19, cex = 1, col = rev(rainbow(np + 
                1)[1:np]), legend = rev(pp), bg = "white", bty = "n")
        }
        else if (types[[k]] == "markov") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k] + 1, 
                ncol = np)
            z = NULL
            for (i in 1:np) {
                Z[[k]][, i] = B[, partbasis] %*% vector.a.ma.schall[partbasis, 
                  i] + intercept[i]
                coefficients[[k]][, i] = Zspathelp[[k]] %*% vector.a.ma.schall[partbasis, 
                  i]
                z = cbind(z, diag(dim(Zspathelp[[k]])[1]) %*% 
                  Zspathelp[[k]] %*% vector.a.ma.schall[partbasis, 
                  i] + intercept[i])
            }
            if (class(bnd[[k]]) != "bnd") {
                plot(seq(0, 1.1 * max(x[[k]]), length = 10), 
                  seq(0, max(z[, i]), length = 10), type = "n", 
                  xlab = "District", ylab = "coefficients")
                points(rep(as.numeric(attr(bnd[[k]], "regions")), 
                  times = np), z[, i], col = rainbow(np + 1)[1:np])
                legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np + 
                  1)[1:np]), legend = rev(pp), bg = "white", 
                  bty = "n")
            }
            else {
                par(mfrow = (c(3, 4)))
                plot.limits = range(z)
                for (i in 1:np) {
                  re = data.frame(cbind(as.numeric(attr(bnd[[k]], 
                    "regions")), z[, i]))
                  drawmap(re, bnd[[k]], regionvar = 1, plotvar = 2, 
                    mar.min = NULL, limits = plot.limits, main = pp[i])
                }
            }
        }
        else if (types[[k]] == "2dspline") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            x.min = apply(x[[k]], 2, min)
            x.max = apply(x[[k]], 2, max)
            x.gitter = cbind(rep(seq(x.min[1], x.max[1], length = 50), 
                times = 50), rep(seq(x.min[2], x.max[2], length = 50), 
                each = 50))
            B.gitter = base(x.gitter, "2dspline")[[1]]
            par(mfrow = (c(3, 4)))
            for (i in 1:np) {
                z <- B.gitter %*% vector.a.ma.schall[partbasis, 
                  i] + intercept[i]
                Z[[k]][, i] = B[, partbasis] %*% vector.a.ma.schall[partbasis, 
                  i] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i]
                z = t(matrix(z, nrow = 50, ncol = 50))
                persp(seq(x.min[1], x.max[1], length = 50), seq(x.min[2], 
                  x.max[2], length = 50), z, ticktype = "detailed", 
                  phi = 40, zlim = range(yy), col = "lightblue", 
                  xlab = "X", ylab = "Y", zlab = "Z", main = pp[i])
            }
        }
        else if (types[[k]] == "radial") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            x.min = apply(x[[k]], 2, min)
            x.max = apply(x[[k]], 2, max)
            x.gitter = cbind(rep(seq(x.min[1], x.max[1], length = 50), 
                times = 50), rep(seq(x.min[2], x.max[2], length = 50), 
                each = 50))
            x[[k]] = x[[k]][order(x[[k]][, 1]), ]
            knots = x[[k]][seq(1, dim(x[[k]])[1], length = min(50, 
                dim(x[[k]])[1])), ]
            B.gitter = matrix(NA, nrow = dim(x.gitter)[1], ncol = dim(knots)[1])
            for (i in 1:dim(x.gitter)[1]) for (j in 1:dim(knots)[1]) {
                r = sqrt(sum((x.gitter[i, ] - knots[j, ])^2))
                B.gitter[i, j] = r^2 * log(r)
            }
            par(mfrow = (c(3, 4)))
            for (i in 1:np) {
                z <- B.gitter %*% vector.a.ma.schall[partbasis, 
                  i] + intercept[i]
                Z[[k]][, i] = B[, partbasis] %*% vector.a.ma.schall[partbasis, 
                  i] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i]
                z = t(matrix(z, nrow = 50, ncol = 50))
                persp(seq(x.min[1], x.max[1], length = 50), seq(x.min[2], 
                  x.max[2], length = 50), z, ticktype = "detailed", 
                  phi = 40, zlim = range(yy), col = "lightblue", 
                  xlab = "X", ylab = "Y", zlab = "Z", main = pp[i])
            }
        }
        else if (types[[k]] == "krig") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            x.min = apply(x[[k]], 2, min)
            x.max = apply(x[[k]], 2, max)
            x.gitter = cbind(rep(seq(x.min[1], x.max[1], length = 50), 
                times = 50), rep(seq(x.min[2], x.max[2], length = 50), 
                each = 50))
            x[[k]] = x[[k]][order(x[[k]][, 1]), ]
            knots = x[[k]][seq(1, dim(x[[k]])[1], length = min(50, 
                dim(x[[k]])[1])), ]
            B.gitter = matrix(NA, nrow = dim(x.gitter)[1], ncol = dim(knots)[1])
            for (i in 1:dim(x.gitter)[1]) for (j in 1:dim(knots)[1]) {
                r = sqrt(sum((x.gitter[i, ] - knots[j, ])^2))/krig.phi[[k]]
                B.gitter[i, j] = exp(-r) * (1 + r)
            }
            par(mfrow = (c(3, 4)))
            for (i in 1:np) {
                z <- B.gitter %*% vector.a.ma.schall[partbasis, 
                  i] + intercept[i]
                Z[[k]][, i] = B[, partbasis] %*% vector.a.ma.schall[partbasis, 
                  i] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i]
                z = t(matrix(z, nrow = 50, ncol = 50))
                persp(seq(x.min[1], x.max[1], length = 50), seq(x.min[2], 
                  x.max[2], length = 50), z, ticktype = "detailed", 
                  phi = 40, zlim = range(yy), col = "lightblue", 
                  xlab = "X", ylab = "Y", zlab = "Z", main = pp[i])
            }
        }
        else if (types[[k]] == "random") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            for (i in 1:np) {
                Z[[k]][, i] <- B[, partbasis] %*% vector.a.ma.schall[partbasis, 
                  i] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i]
            }
            plot(seq(0, 1.1 * max(x[[k]]), length = 10), seq(0, 
                max(coefficients[[k]] + intercept), length = 10), 
                type = "n", xlab = "Group", ylab = "coefficients")
            points(rep(sort(unique(x[[k]])), times = np), coefficients[[k]] + 
                intercept, col = rainbow(np + 1)[1:np])
            legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np + 
                1)[1:np]), legend = rev(pp), bg = "white", bty = "n")
        }
        else if (types[[k]] == "ridge") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            for (i in 1:np) {
                Z[[k]][, i] <- B[, partbasis] %*% vector.a.ma.schall[partbasis, 
                  i] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i]
            }
            plot(seq(0, 1.1 * dim(x[[k]])[2], length = 10), seq(0, 
                max(coefficients[[k]] + intercept), length = 10), 
                type = "n", xlab = "X variables", ylab = "coefficients")
            points(rep(1:dim(x[[k]])[2], times = np), coefficients[[k]] + 
                intercept, col = rainbow(np + 1)[1:np])
            legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np + 
                1)[1:np]), legend = rev(pp), bg = "white", bty = "n")
        }
    }
    result = list(lambda = final.lambdas, intercepts = intercept, 
        coefficients = coefficients, residual = b, asymmetry = cc, 
        values = Z, response = yy, covariates = x, formula = formula)
    class(result) = c("expectreg", "bundle")
    result
}
