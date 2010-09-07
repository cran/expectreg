expectile.sheets <-
function (formula, data = NULL, smooth = c("acv", "none"), lambda = 0.1, 
    lambdap = 5, expectiles = NA, density = FALSE) 
{
    smooth = match.arg(smooth)
    if (density) {
        pp <- seq(0.01, 0.99, by = 0.01)
        pp.plot <- c(1, 2, 5, 10, 20, 50, 80, 90, 95, 98, 99)
        row.grid = 3
        col.grid = 4
    }
    else if (any(is.na(expectiles)) || !is.vector(expectiles) || 
        any(expectiles > 1) || any(expectiles < 0)) {
        pp <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 
            0.98, 0.99)
        pp.plot <- 1:length(pp)
        row.grid = 3
        col.grid = 4
    }
    else {
        pp <- expectiles
        pp.plot <- 1:length(pp)
        row.grid = floor(sqrt(length(pp)))
        col.grid = ceiling(sqrt(length(pp)))
        if (length(pp) > row.grid * col.grid) 
            row.grid = row.grid + 1
    }
    np <- length(pp)
    np.plot <- length(pp.plot)
    yy = eval(parse(text = formula[2]), envir = data, enclos = .GlobalEnv)
    m = length(yy)
    design = list()
    x = list()
    types = list()
    bnd = list()
    Zspathelp = list()
    nb = vector()
    krig.phi = list()
    B = list()
    DD = list()
    for (i in 1:length(labels(terms(formula)))) {
        types[[i]] = strsplit(labels(terms(formula))[i], "(", 
            fixed = TRUE)[[1]][1]
        if (types[[i]] == labels(terms(formula))[i]) {
            design[[i]] = base(types[[i]], "parametric")
            types[[i]] = "parametric"
        }
        else design[[i]] = eval(parse(text = labels(terms(formula))[i]), 
            envir = data, enclos = .GlobalEnv)
    }
    nterms = length(design)
    B[[1]] = design[[1]][[1]]
    DD[[1]] = design[[1]][[2]]
    x[[1]] = design[[1]][[3]]
    types[[1]] = design[[1]][[4]]
    bnd[[1]] = design[[1]][[5]]
    Zspathelp[[1]] = design[[1]][[6]]
    nb[1] = ncol(design[[1]][[1]])
    krig.phi[[1]] = design[[1]][[7]]
    if (length(labels(terms(formula))) > 1) 
        for (i in 2:length(labels(terms(formula)))) {
            B[[i]] = design[[i]][[1]]
            DD[[i]] = design[[i]][[2]]
            x[[i]] = design[[i]][[3]]
            types[[i]] = design[[i]][[4]]
            bnd[[i]] = design[[i]][[5]]
            Zspathelp[[i]] = design[[i]][[6]]
            nb[i] = ncol(design[[i]][[1]])
            krig.phi[[i]] = design[[i]][[7]]
        }
    lala <- matrix(c(rep(lambda, nterms), rep(lambdap, nterms)), 
        nrow = nterms, ncol = 2, dimnames = list(1:nterms, c("curve", 
            "sheet")))
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
        acv.min = nlm(acv.sheets, p = lala, yy = yy, B = B, pp = pp, 
            DD = DD, nb = nb, ndigit = 8, iterlim = 50, gradtol = 1e-04)
        min.lambda = matrix(abs(acv.min$estimate), ncol = 2)
        lala[, 1] <- min.lambda[, 1]
        lala[, 2] <- min.lambda[, 2]
        ynp <- rep(yy, np)
        ps <- rep(pp, each = m)
        vv <- log(ps/(1 - ps))
        w <- runif(m * np)
        p2f.new <- pspfit2d.new(B, DD, ps, ynp, w, lala[, 1], 
            lala[, 2])
        vector.a.ma.schall <- matrix(p2f.new$coef, ncol = np, 
            byrow = T)
    }
    else {
        ynp <- rep(yy, np)
        ps <- rep(pp, each = m)
        vv <- log(ps/(1 - ps))
        w <- runif(m * np)
        p2f.new <- pspfit2d.new(B, DD, ps, ynp, w, lala[, 1], 
            lala[, 2])
    }
    curves = p2f.new$curves
    Z <- list()
    coefficients <- list()
    final.lambdas <- list()
    intercept = p2f.new$intercept
    for (k in 1:length(design)) {
        final.lambdas[[k]] = lala[k, ]
        partbasis = (sum(nb[0:(k - 1)]) + 1):(sum(nb[0:k]))
        dev.new()
        if (types[[k]] == "pspline") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            plot(x[[k]], yy, cex = 0.5, pch = 20, col = "grey42", 
                xlab = "x", ylab = "y", ylim = range(cbind(yy, 
                  curves[[k]]), na.rm = TRUE))
            matlines(sort(x[[k]])[seq(1, m, length = min(m, 100))], 
                curves[[k]][order(x[[k]])[seq(1, m, length = min(m, 
                  100))], pp.plot], col = rainbow(np.plot + 1)[1:np.plot], 
                lty = 1)
            legend(x = "topright", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                1)[1:np.plot]), legend = rev(pp), bg = "white", 
                bty = "n")
        }
        else if (types[[k]] == "markov") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k] + 1, 
                ncol = np)
            if (class(bnd[[k]]) != "bnd") {
                plot(seq(0, 1.1 * max(x[[k]]), length = 10), 
                  seq(0, max(z[, pp.plot]), length = 10), type = "n", 
                  xlab = "District", ylab = "coefficients")
                matpoints(sort(x[[k]]), curves[[k]][order(x[[k]]), 
                  pp.plot], col = rainbow(np.plot + 1)[1:np.plot])
                legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                  1)[1:np.plot]), legend = rev(pp.plot/100), 
                  bg = "white", bty = "n")
            }
            else {
                par(mfrow = (c(row.grid, col.grid)))
                plot.limits = range(z)
                for (i in 1:np.plot) {
                  re = data.frame(cbind(x[[k]], curves[[k]][, 
                    pp.plot[i]]))
                  drawmap(re, bnd[[k]], regionvar = 1, plotvar = 2, 
                    mar.min = NULL, limits = plot.limits, main = pp.plot[i]/100)
                }
            }
        }
        else if (types[[k]] == "2dspline") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            par(mfrow = (c(row.grid, col.grid)))
            for (i in 1:np) {
                if (i %in% pp.plot) {
                  z = interp(x[[k]][, 1], x[[k]][, 2], curves[[k]][, 
                    i])
                  persp(z[[1]], z[[2]], z[[3]], ticktype = "detailed", 
                    phi = 40, zlim = range(yy), col = "lightblue", 
                    xlab = "X", ylab = "Y", zlab = "Z", main = pp.plot[i]/100)
                }
            }
        }
        else if (types[[k]] == "radial") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            par(mfrow = (c(row.grid, col.grid)))
            for (i in 1:np) {
                if (i %in% pp.plot) {
                  z = interp(x[[k]][, 1], x[[k]][, 2], curves[[k]][, 
                    i])
                  persp(z[[1]], z[[2]], z[[3]], ticktype = "detailed", 
                    phi = 40, zlim = range(yy), col = "lightblue", 
                    xlab = "X", ylab = "Y", zlab = "Z", main = pp.plot[i]/100)
                }
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
            par(mfrow = (c(row.grid, col.grid)))
            for (i in 1:np) {
                Z[[k]][, i] = B[[k]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
                if (i %in% pp.plot) {
                  z <- B.gitter %*% vector.a.ma.schall[partbasis, 
                    i, drop = FALSE] + intercept[i]
                  z = t(matrix(z, nrow = 50, ncol = 50))
                  persp(seq(x.min[1], x.max[1], length = 50), 
                    seq(x.min[2], x.max[2], length = 50), z, 
                    ticktype = "detailed", phi = 40, zlim = range(yy), 
                    col = "lightblue", xlab = "X", ylab = "Y", 
                    zlab = "Z", main = pp.plot[i]/100)
                }
            }
        }
        else if (types[[k]] == "random") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            for (i in 1:np) {
                Z[[k]][, i] <- B[[k]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
            }
            plot(seq(0, 1.1 * max(x[[k]]), length = 10), seq(0, 
                max(coefficients[[k]] + intercept), length = 10), 
                type = "n", xlab = "Group", ylab = "coefficients")
            points(rep(sort(unique(x[[k]])), times = np.plot), 
                (coefficients[[k]] + intercept)[, pp.plot], col = rainbow(np.plot + 
                  1)[1:np.plot])
            legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                1)[1:np.plot]), legend = rev(pp.plot/100), bg = "white", 
                bty = "n")
        }
        else if (types[[k]] == "ridge") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            for (i in 1:np) {
                Z[[k]][, i] <- B[[k]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
            }
            plot(seq(0, 1.1 * dim(x[[k]])[2], length = 10), seq(0, 
                max(coefficients[[k]] + intercept), length = 10), 
                type = "n", xlab = "X variables", ylab = "coefficients")
            points(rep(1:dim(x[[k]])[2], times = np.plot), (coefficients[[k]] + 
                intercept)[, pp.plot], col = rainbow(np.plot + 
                1)[1:np.plot])
            legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                1)[1:np.plot]), legend = rev(pp.plot/100), bg = "white", 
                bty = "n")
        }
        else if (types[[k]] == "parametric") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            for (i in 1:np) {
                Z[[k]][, i] <- B[[k]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
            }
            plot(x[[k]], yy, cex = 0.5, pch = 20, col = "grey42", 
                xlab = "x", ylab = "y", ylim = range(cbind(yy, 
                  Z[[k]])))
            matlines(sort(x[[k]]), Z[[k]][order(x[[k]]), pp.plot], 
                col = rainbow(np.plot + 1)[1:np.plot], lty = 1)
            legend(x = "bottomright", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                1)[1:np.plot]), legend = rev(pp), bg = "white", 
                bty = "n")
        }
    }
    result = list(lambda = final.lambdas, intercepts = intercept, 
        values = curves, response = yy, covariates = x, formula = formula)
    class(result) = c("expectreg", "sheets")
    result
}
