expectile.restricted <-
function (formula, data = NULL, smooth = c("schall", "acv", "fixed"), 
    lambda = 0.1, expectiles = NA, density = FALSE) 
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
    center = TRUE
    varying = list()
    for (i in 1:length(labels(terms(formula)))) {
        types[[i]] = strsplit(labels(terms(formula))[i], "(", 
            fixed = TRUE)[[1]][1]
        if (types[[i]] == labels(terms(formula))[i]) {
            design[[i]] = base(matrix(eval(parse(text = labels(terms(formula))[i]), 
                envir = data, enclos = .GlobalEnv), nrow = m), 
                "parametric")
            formula = eval(substitute(update(formula, . ~ variable2 + 
                . - variable1), list(variable1 = as.name(types[[i]]), 
                variable2 = as.name(paste("base(", types[[i]], 
                  ",'parametric')", sep = "")))))
            types[[i]] = "parametric"
        }
        else design[[i]] = eval(parse(text = labels(terms(formula))[i]), 
            envir = data, enclos = .GlobalEnv)
    }
    nterms = length(design)
    varying[[1]] = design[[1]][[9]]
    if (any(!is.na(varying[[1]]))) 
        B = design[[1]][[1]] * varying[[1]]
    else B = design[[1]][[1]]
    DD = as.matrix(design[[1]][[2]])
    x[[1]] = design[[1]][[3]]
    types[[1]] = design[[1]][[4]]
    bnd[[1]] = design[[1]][[5]]
    Zspathelp[[1]] = design[[1]][[6]]
    nb[1] = ncol(design[[1]][[1]])
    krig.phi[[1]] = design[[1]][[7]]
    center = center && design[[1]][[8]]
    if (length(labels(terms(formula))) > 1) 
        for (i in 2:length(labels(terms(formula)))) {
            varying[[i]] = design[[i]][[9]]
            if (any(!is.na(varying[[i]]))) 
                B = cbind(B, design[[i]][[1]] * varying[[i]])
            else B = cbind(B, design[[i]][[1]])
            design[[i]][[2]] = as.matrix(design[[i]][[2]])
            DD = rbind(cbind(DD, matrix(0, nrow = nrow(DD), ncol = ncol(design[[i]][[2]]))), 
                cbind(matrix(0, nrow = nrow(design[[i]][[2]]), 
                  ncol = ncol(DD)), design[[i]][[2]]))
            x[[i]] = design[[i]][[3]]
            types[[i]] = design[[i]][[4]]
            bnd[[i]] = design[[i]][[5]]
            Zspathelp[[i]] = design[[i]][[6]]
            nb[i] = ncol(design[[i]][[1]])
            krig.phi[[i]] = design[[i]][[7]]
            center = center && design[[i]][[8]]
        }
    if (center) {
        B = cbind(1, B)
        DD = rbind(0, cbind(0, DD))
    }
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
                nb)
            mean.coefficients <- aa$a
            sig2 <- vector()
            tau2 <- vector()
            l0 <- lala[, 1]
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
                z <- B %*% aa$a
                H = solve(t(partB) %*% (w[, i] * partB) + lala[i, 
                  1] * t(partDD) %*% partDD)
                H = apply(sqrt(w[, i]) * partB, 1, function(x) {
                  t(x) %*% H %*% x
                })
                sig2[i] <- sum(w[, i] * (yy - z)^2, na.rm = TRUE)/(m - 
                  sum(aa$diag.hat.ma))
                tau2[i] <- sum(v^2)/sum(H) + 1e-06
                lala[i, 1] <- max(sig2[i]/tau2[i], 1e-10)
            }
            dc <- max(abs(log10(l0) - log10(lala[, 1])))
            it = it + 1
        }
        if (it == 100) 
            warning("Schall algorithm did not converge. Stopping after 100 iterations.")
    }
    else if (smooth == "acv") {
        acv.min = nlm(acv, p = lala[, 1], yy = yy, B = B, quantile = 0.5, 
            DD = DD, nb = nb, ndigit = 8, iterlim = 50, gradtol = 1e-04)
        aa <- asyregpen.lsfit(yy, B, 0.5, abs(acv.min$estimate), 
            DD, nb)
        mean.coefficients <- aa$a
        lala[, 1] <- abs(acv.min$estimate)
    }
    else {
        aa <- asyregpen.lsfit(yy, B, 0.5, lala[, 1], DD, nb)
        mean.coefficients <- aa$a
    }
    residuals = yy - B %*% mean.coefficients
    gg = asyregpen.lsfit(abs(residuals), B, 0.5, lala[, 1], DD, 
        nb)
    cc = NULL
    for (q in 1:np) {
        ca = asyregpen.lsfit(residuals, gg$fitted, pp[q], NULL, 
            matrix(0, nrow = 1, ncol = 1), 0)
        cc[q] = ca$a
        vector.a.ma.schall[, q] = mean.coefficients + ca$a * 
            gg$a
    }
    Z <- list()
    coefficients <- list()
    final.lambdas <- list()
    helper <- list()
    fitted = B %*% vector.a.ma.schall
    if (center) {
        intercept = vector.a.ma.schall[1, ]
        B = B[, -1, drop = FALSE]
        vector.a.ma.schall = vector.a.ma.schall[-1, , drop = FALSE]
    }
    else intercept = rep(0, np)
    for (k in 1:length(design)) {
        final.lambdas[[k]] = lala[k, ]
        partbasis = (sum(nb[0:(k - 1)]) + 1):(sum(nb[0:k]))
        dev.new()
        if (types[[k]] == "pspline") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = NA
            for (i in 1:np) {
                Z[[k]][, i] <- design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
            }
            plot(x[[k]], yy, cex = 0.5, pch = 20, col = "grey42", 
                xlab = "x", ylab = "y", ylim = range(cbind(yy, 
                  Z[[k]]), na.rm = TRUE))
            matlines(sort(x[[k]])[seq(1, m, length = min(m, 100))], 
                Z[[k]][order(x[[k]])[seq(1, m, length = min(m, 
                  100))], pp.plot], col = rainbow(np.plot + 1)[1:np.plot], 
                lty = 1)
            legend(x = "topright", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                1)[1:np.plot]), legend = rev(pp[pp.plot]), bg = "white", 
                bty = "n")
        }
        else if (types[[k]] == "markov") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k] + 1, 
                ncol = np)
            z = NULL
            helper[[k]] = list(bnd[[k]], Zspathelp[[k]])
            for (i in 1:np) {
                Z[[k]][, i] = design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = Zspathelp[[k]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
                z = cbind(z, diag(dim(Zspathelp[[k]])[1]) %*% 
                  Zspathelp[[k]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i])
            }
            if (class(bnd[[k]]) != "bnd") {
                plot(seq(0, 1.1 * max(x[[k]]), length = 10), 
                  seq(0, max(z[, pp.plot]), length = 10), type = "n", 
                  xlab = "District", ylab = "coefficients")
                points(rep(as.numeric(attr(bnd[[k]], "regions")), 
                  times = np), z[, pp.plot], col = rainbow(np.plot + 
                  1)[1:np.plot])
                legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                  1)[1:np.plot]), legend = rev(pp[pp.plot]), 
                  bg = "white", bty = "n")
            }
            else {
                par(mfrow = (c(row.grid, col.grid)))
                plot.limits = range(z)
                for (i in 1:np.plot) {
                  re = data.frame(cbind(as.numeric(attr(bnd[[k]], 
                    "regions")), z[, pp.plot[i]]))
                  drawmap(re, bnd[[k]], regionvar = 1, plotvar = 2, 
                    mar.min = NULL, limits = plot.limits, main = pp[pp.plot][i])
                }
            }
        }
        else if (types[[k]] == "2dspline") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = NA
            x.min = apply(x[[k]], 2, min)
            x.max = apply(x[[k]], 2, max)
            x.gitter = cbind(rep(seq(x.min[1], x.max[1], length = 50), 
                times = 50), rep(seq(x.min[2], x.max[2], length = 50), 
                each = 50))
            B.gitter = base(x.gitter, "2dspline")[[1]]
            par(mfrow = (c(row.grid, col.grid)))
            for (i in 1:np) {
                Z[[k]][, i] = design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
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
                    zlab = "Z", main = pp[pp.plot][i])
                }
            }
        }
        else if (types[[k]] == "radial") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = NA
            x.min = apply(x[[k]], 2, min)
            x.max = apply(x[[k]], 2, max)
            x.gitter = cbind(rep(seq(x.min[1], x.max[1], length = 50), 
                times = 50), rep(seq(x.min[2], x.max[2], length = 50), 
                each = 50))
            x.ord = x[[k]][order(x[[k]][, 1]), ]
            knots = x.ord[seq(1, dim(x[[k]])[1], length = min(50, 
                dim(x[[k]])[1])), ]
            B.gitter = matrix(NA, nrow = dim(x.gitter)[1], ncol = dim(knots)[1])
            for (i in 1:dim(x.gitter)[1]) for (j in 1:dim(knots)[1]) {
                r = sqrt(sum((x.gitter[i, ] - knots[j, ])^2))
                B.gitter[i, j] = r^2 * log(r)
            }
            par(mfrow = (c(row.grid, col.grid)))
            for (i in 1:np) {
                Z[[k]][, i] = design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
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
                    zlab = "Z", main = pp[pp.plot][i])
                }
            }
        }
        else if (types[[k]] == "krig") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = krig.phi[[k]]
            x.min = apply(x[[k]], 2, min)
            x.max = apply(x[[k]], 2, max)
            x.gitter = cbind(rep(seq(x.min[1], x.max[1], length = 50), 
                times = 50), rep(seq(x.min[2], x.max[2], length = 50), 
                each = 50))
            x.ord = x[[k]][order(x[[k]][, 1]), ]
            knots = x.ord[seq(1, dim(x[[k]])[1], length = min(50, 
                dim(x[[k]])[1])), ]
            B.gitter = matrix(NA, nrow = dim(x.gitter)[1], ncol = dim(knots)[1])
            for (i in 1:dim(x.gitter)[1]) for (j in 1:dim(knots)[1]) {
                r = sqrt(sum((x.gitter[i, ] - knots[j, ])^2))/krig.phi[[k]]
                B.gitter[i, j] = exp(-r) * (1 + r)
            }
            par(mfrow = (c(row.grid, col.grid)))
            for (i in 1:np) {
                Z[[k]][, i] = design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
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
                    zlab = "Z", main = pp[pp.plot][i])
                }
            }
        }
        else if (types[[k]] == "random") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = NA
            for (i in 1:np) {
                Z[[k]][, i] <- design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
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
                1)[1:np.plot]), legend = rev(pp[pp.plot]), bg = "white", 
                bty = "n")
        }
        else if (types[[k]] == "ridge") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = NA
            for (i in 1:np) {
                Z[[k]][, i] <- design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
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
                1)[1:np.plot]), legend = rev(pp[pp.plot]), bg = "white", 
                bty = "n")
        }
        else if (types[[k]] == "parametric") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = NA
            for (i in 1:np) {
                Z[[k]][, i] <- design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
            }
            matplot(1:nb[k], coefficients[[k]][, pp.plot], col = rainbow(np.plot + 
                1)[1:np.plot], pch = 15)
            legend(x = "bottomright", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                1)[1:np.plot]), legend = rev(pp[pp.plot]), bg = "white", 
                bty = "n")
        }
        else if (types[[k]] == "special") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = NA
            for (i in 1:np) {
                Z[[k]][, i] <- design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
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
                1)[1:np.plot]), legend = rev(pp[pp.plot]), bg = "white", 
                bty = "n")
        }
    }
    result = list(lambda = final.lambdas, intercepts = intercept, 
        coefficients = coefficients, trend.coef = mean.coefficients, 
        residual.coef = gg$a, asymmetry = cc, values = Z, response = yy, 
        covariates = x, formula = formula, expectiles = pp, effects = types, 
        helper = helper, design = cbind(1, B), fitted = fitted)
    class(result) = c("expectreg", "restricted")
    result
}
