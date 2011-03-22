expectile.bundle <-
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
    yy = eval(parse(text = formula[2]), envir = data, enclos = environment(formula))
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
                envir = data, enclos = environment(formula)), 
                nrow = m), "parametric")
            formula = eval(substitute(update(formula, . ~ variable2 + 
                . - variable1), list(variable1 = as.name(types[[i]]), 
                variable2 = as.name(paste("base(", types[[i]], 
                  ",'parametric')", sep = "")))))
            types[[i]] = "parametric"
        }
        else design[[i]] = eval(parse(text = labels(terms(formula))[i]), 
            envir = data, enclos = environment(formula))
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
    vector.a.ma.schall <- matrix(NA, nrow = sum(nb) + 1 * center, 
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
            DD = DD, nb = nb, ndigit = 8, iterlim = 50, gradtol = 1e-04)
        aa <- asyregpen.lsfit(yy, B, 0.5, abs(acv.min$estimate), 
            DD, nb)
        mean.coefficients <- aa$a
        lala[, 1] <- abs(acv.min$estimate)
        residuals = yy - B %*% mean.coefficients
        acv.min = nlm(acv, p = lala[, 2], yy = residuals, B = B, 
            quantile = 0.5, DD = DD, nb = nb, ndigit = 8, iterlim = 50, 
            gradtol = 1e-04)
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
        aa <- asyregpen.lsfit(yy, B, 0.5, lala[, 1], DD, nb)
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
        }
        else if (types[[k]] == "markov") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k] + 1 * 
                center, ncol = np)
            helper[[k]] = list(bnd[[k]], Zspathelp[[k]])
            for (i in 1:np) {
                Z[[k]][, i] = design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = Zspathelp[[k]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
            }
        }
        else if (types[[k]] == "2dspline") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = NA
            for (i in 1:np) {
                Z[[k]][, i] = design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
            }
        }
        else if (types[[k]] == "radial") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = NA
            for (i in 1:np) {
                Z[[k]][, i] = design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
            }
        }
        else if (types[[k]] == "krig") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = krig.phi[[k]]
            for (i in 1:np) {
                Z[[k]][, i] = design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
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
        }
    }
    result = list(lambda = final.lambdas, intercepts = intercept, 
        coefficients = coefficients, trend.coef = mean.coefficients, 
        residual.coef = b, asymmetry = cc, values = Z, response = yy, 
        covariates = x, formula = formula, expectiles = pp, effects = types, 
        helper = helper, design = cbind(1, B), fitted = fitted)
    class(result) = c("expectreg", "bundle")
    result
}
