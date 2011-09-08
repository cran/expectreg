expectile.laws <-
function (formula, data = NULL, smooth = c("schall", "acv", "fixed"), 
    lambda = 0.1, expectiles = NA, parallel = FALSE) 
{
    smooth = match.arg(smooth)
    if (any(is.na(expectiles)) || !is.vector(expectiles) || any(expectiles > 
        1) || any(expectiles < 0)) {
        pp <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 
            0.98, 0.99)
        row.grid = 3
        col.grid = 4
    }
    else {
        pp <- expectiles
        row.grid = floor(sqrt(length(pp)))
        col.grid = ceiling(sqrt(length(pp)))
        if (length(pp) > row.grid * col.grid) 
            row.grid = row.grid + 1
    }
    np <- length(pp)
    yy = eval(as.expression(formula[[2]]), envir = data, enclos = environment(formula))
    attr(yy, "name") = deparse(formula[[2]])
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
    if (formula[[3]] == "1") {
        design[[1]] = base(matrix(1, nrow = m, ncol = 1), "parametric", 
            center = F)
        smooth = "fixed"
    }
    else if (formula[[3]] == ".") {
        design[[1]] = base(data[, names(data) != all.vars(formula[[2]])], 
            "parametric")
        smooth = "fixed"
    }
    else for (i in 1:length(labels(terms(formula)))) {
        types[[i]] = strsplit(labels(terms(formula))[i], "(", 
            fixed = TRUE)[[1]][1]
        if (types[[i]] == labels(terms(formula))[i]) {
            design[[i]] = base(matrix(eval(parse(text = labels(terms(formula))[i]), 
                envir = data, enclos = environment(formula)), 
                nrow = m), "parametric")
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
    names(x)[1] = design[[1]]$xname
    types[[1]] = design[[1]][[4]]
    bnd[[1]] = design[[1]][[5]]
    Zspathelp[[1]] = design[[1]][[6]]
    nb[1] = ncol(design[[1]][[1]])
    krig.phi[[1]] = design[[1]][[7]]
    center = center && design[[1]][[8]]
    if (length(design) > 1) 
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
            names(x)[i] = design[[i]]$xname
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
    myapply <- lapply
    if (parallel && .Platform$OS.type == "unix" && require("multicore")) {
        if (!multicore:::isChild()) {
            myapply <- mclapply
        }
    }
    dummy.reg <- function(pp, lambda, smooth, yy, B, DD, nb, 
        nterms, center) {
        print(paste("Expectile:", pp, sep = " "))
        if (length(lambda) < nterms) 
            lala = rep(lambda[1], nterms)
        else lala = lambda
        if (smooth == "schall") {
            dc = 1
            dw = 1
            w <- rep(1, times = m)
            it = 1
            while ((dc >= 0.01 || dw != 0) && it < 100) {
                aa <- asyregpen.lsfit(yy, B, pp, lala, DD, nb)
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
                dc <- max(abs(log10(l0) - log10(lala)))
                dw <- sum(w != w0, na.rm = TRUE)
                it = it + 1
            }
            if (it == 100) 
                warning("Schall algorithm did not converge. Stopping after 100 iterations.")
        }
        else if (smooth == "acv") {
            acv.min = nlm(acv, p = lala, yy = yy, B = B, quantile = pp, 
                DD = DD, nb = nb, ndigit = 8, iterlim = 100, 
                gradtol = 1e-04)
            aa <- asyregpen.lsfit(yy, B, pp, abs(acv.min$estimate), 
                DD, nb)
            vector.a.ma.schall <- aa$a
            lala <- abs(acv.min$estimate)
        }
        else {
            aa <- asyregpen.lsfit(yy, B, pp, lala, DD, nb)
            vector.a.ma.schall <- aa$a
        }
        list(vector.a.ma.schall, lala)
    }
    coef.vector = myapply(pp, function(pp) dummy.reg(pp, lambda, 
        smooth, yy, B, DD, nb, nterms, center))
    lala <- matrix(lambda, nrow = nterms, ncol = np)
    vector.a.ma.schall <- matrix(NA, nrow = sum(nb) + (1 * center), 
        ncol = np)
    for (i in 1:np) {
        vector.a.ma.schall[, i] = coef.vector[[i]][[1]]
        lala[, i] = coef.vector[[i]][[2]]
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
        names(final.lambdas)[k] = design[[k]]$xname
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
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = list(bnd[[k]], Zspathelp[[k]])
            for (i in 1:np) {
                Z[[k]][, i] = design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
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
        names(Z)[k] = design[[k]]$xname
        names(coefficients)[k] = design[[k]]$xname
    }
    result = list(lambda = final.lambdas, intercepts = intercept, 
        coefficients = coefficients, values = Z, response = yy, 
        covariates = x, formula = formula, expectiles = pp, effects = types, 
        helper = helper, design = cbind(1, B), fitted = fitted)
    result$predict <- function(newdata = NULL) {
        BB = list()
        values = list()
        bmat = NULL
        for (k in 1:length(coefficients)) {
            BB[[k]] = predict(design[[k]], newdata)
            values[[k]] <- BB[[k]] %*% coefficients[[k]]
            values[[k]] = t(apply(values[[k]], 1, function(x) {
                x + intercept
            }))
            bmat = cbind(bmat, BB[[k]])
        }
        if (center) 
            bmat = cbind(1, bmat)
        fitted = bmat %*% rbind(intercept, vector.a.ma.schall)
        names(values) = names(coefficients)
        list(fitted = fitted, values = values)
    }
    class(result) = c("expectreg", "laws")
    result
}
