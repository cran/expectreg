expectile.boost <-
function (formula, data = NULL, mstop = NA, expectiles = NA, 
    parallel = FALSE, cv = TRUE) 
{
    require(mboost)
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
    if (any(is.na(mstop))) 
        mstop = rep(4000, np)
    else if (length(mstop) == 1) 
        mstop = rep(mstop, np)
    else if (length(mstop) < np) 
        mstop = c(mstop, rep(max(mstop), times = np - length(mstop)))
    blsstr <- labels(terms(formula))
    types = list()
    bls <- list()
    x <- list()
    bnd = list()
    yy = eval(parse(text = formula[2]), envir = data, enclos = environment(formula))
    for (i in 1:length(blsstr)) {
        types[[i]] = strsplit(blsstr[i], "(", fixed = TRUE)[[1]][1]
        if (types[[i]] == blsstr[i]) {
            types[[i]] = "parametric"
            x[[i]] = eval(parse(text = blsstr[i]), envir = data, 
                enclos = environment(formula))
            bls[[i]] = NULL
        }
        else {
            bls[[i]] <- eval(parse(text = blsstr[i]), envir = data, 
                enclos = environment(formula))
            x[[i]] <- eval(parse(text = bls[[i]]$get_names()), 
                envir = data, enclos = environment(formula))
        }
        if (types[[i]] != "bmrf") 
            bnd[[i]] = NA
        else bnd[[i]] = environment(bls[[i]]$dpp)$args$bnd
    }
    m = length(yy)
    kfld <- 10
    ntest <- floor(m/kfld)
    cv10f <- matrix(c(rep(c(rep(0, ntest), rep(1, m)), kfld - 
        1), rep(0, m * kfld - (kfld - 1) * (m + ntest))), nrow = m)
    values <- list()
    z <- list()
    helper <- list()
    for (k in 1:length(blsstr)) {
        values[[k]] = matrix(NA, nrow = m, ncol = np)
        fitted = matrix(NA, nrow = m, ncol = np)
        helper[[k]] = NA
        if (types[[k]] == "bbs") {
            z[[k]] = values[[k]]
            types[[k]] = "pspline"
        }
        else if (types[[k]] == "parametric") 
            z[[k]] = values[[k]]
        else if (types[[k]] == "bmrf") {
            z[[k]] = matrix(NA, nrow = length(attr(bnd[[k]], 
                "regions")), ncol = np)
            types[[k]] = "markov"
        }
        else if (types[[k]] == "brandom") {
            z[[k]] = matrix(NA, nrow = length(unique(x[[k]])), 
                ncol = np)
            types[[k]] = "random"
        }
        else if (types[[k]] == "bspatial") {
            z[[k]] = list()
            types[[k]] = "2dspline"
        }
    }
    myapply <- lapply
    if (parallel && .Platform$OS.type == "unix" && require("multicore")) {
        if (!multicore:::isChild()) {
            myapply <- mclapply
        }
    }
    dummy.reg <- function(p, formula, data, mstop, pp, cv10f, 
        types, x, blsstr, bnd) {
        values <- list()
        z <- list()
        inb <- gamboost(formula = formula, data = data, control = boost_control(mstop = mstop[p], 
            nu = 0.1, risk = "inbag"), family = ExpectReg(pp[p]))
        if (cv) {
            cvr = cvrisk(inb, folds = cv10f)
            print(paste("Expectile", pp[p], ", Boosting iterations:", 
                mstop(cvr), sep = " "))
            inb = inb[mstop(cvr)]
        }
        else print(paste("Expectile", pp[p], sep = " "))
        fitted = fitted(inb)
        for (k in 1:length(blsstr)) {
            if (types[[k]] == "pspline") {
                independent = which(dimnames(data)[[2]] == bls[[k]]$get_names())
                d.tmp = data
                for (j in (1:dim(d.tmp)[2])[-independent]) d.tmp[, 
                  j] = NA
                values[[k]] = predict(inb, d.tmp, which = k) + 
                  inb$offset
                z[[k]] = values[[k]]
            }
            else if (types[[k]] == "2dspline") {
                x[[k]] = data[, bls[[k]]$get_names()]
                x.min = apply(x[[k]], 2, min)
                x.max = apply(x[[k]], 2, max)
                x.gitter = cbind(rep(seq(x.min[1], x.max[1], 
                  length = 20), times = 20), rep(seq(x.min[2], 
                  x.max[2], length = 20), each = 20))
                d.tmp = matrix(0, nrow = 400, ncol = dim(data)[2])
                for (j in (1:dim(d.tmp)[2])[-independent]) d.tmp[, 
                  j] = NA
                d.val = data
                independent = c(which(dimnames(data)[[2]] == 
                  bls[[k]]$get_names()[1]), which(dimnames(data)[[2]] == 
                  bls[[k]]$get_names()[2]))
                d.tmp[, independent] = x.gitter
                dimnames(d.tmp) = list(1:400, dimnames(data)[[2]])
                for (j in 1:dim(d.val)[2][-independent]) d.val[, 
                  j] = NA
                z[[k]] <- predict(inb, data.frame(d.tmp), which = k) + 
                  inb$offset
                values[[k]] = predict(inb, which = k) + inb$offset
                z[[k]] = t(matrix(z[[k]], nrow = 20, ncol = 20))
            }
            else if (types[[k]] == "random") {
                districts = unique(x[[k]])
                d.tmp = matrix(0, nrow = length(districts), ncol = dim(data)[2])
                dimnames(d.tmp) = list(1:length(districts), dimnames(data)[[2]])
                independent = which(dimnames(data)[[2]] == bls[[k]]$get_names())
                d.tmp[, independent] = districts
                for (j in (1:dim(d.tmp)[2])[-independent]) d.tmp[, 
                  j] = NA
                d.val = data
                for (j in 1:dim(d.val)[2][-independent]) d.val[, 
                  j] = NA
                values[[k]] = predict(inb, which = k) + inb$offset
                z[[k]] <- predict(inb, d.tmp, which = k) + inb$offset
            }
            else if (types[[k]] == "markov") {
                districts = as.numeric(attr(bnd[[k]], "regions"))
                d.tmp = data.frame(matrix(0, nrow = length(districts), 
                  ncol = dim(data)[2]))
                dimnames(d.tmp) = list(1:length(districts), dimnames(data)[[2]])
                independent = which(dimnames(data)[[2]] == bls[[k]]$get_names())
                d.tmp[, independent] = districts
                for (j in (1:dim(d.tmp)[2])[-independent]) d.tmp[, 
                  j] = NA
                d.val = data
                for (j in 1:dim(d.val)[2][-independent]) d.val[, 
                  j] = NA
                values[[k]] = predict(inb, d.val, which = k) + 
                  inb$offset
                z[[k]] <- predict(inb, d.tmp, which = k) + inb$offset
            }
            else if (types[[k]] == "parametric") {
                independent = which(dimnames(data)[[2]] == blsstr[k])
                d.tmp = data
                for (j in (1:dim(d.tmp)[2])[-independent]) d.tmp[, 
                  j] = NA
                values[[k]] = predict(inb, which = k) + inb$offset
                z[[k]] = values[[k]]
            }
        }
        gc()
        list(values, z, fitted, inb)
    }
    coef.vector = myapply(1:np, function(i) dummy.reg(i, formula, 
        data, mstop, pp, cv10f, types, x, blsstr, bnd))
    boost.object = list()
    for (i in 1:np) {
        boost.object[[i]] = coef.vector[[i]][[4]]
        fitted[, i] = coef.vector[[i]][[3]]
        for (k in 1:length(blsstr)) {
            values[[k]][, i] = coef.vector[[i]][[1]][[k]]
            if (types[[k]] == "2dspline") 
                z[[k]][[i]] = coef.vector[[i]][[2]][[k]]
            else z[[k]][, i] = coef.vector[[i]][[2]][[k]]
        }
    }
    result = list(values = values, response = yy, covariates = x, 
        formula = formula, expectiles = pp, effects = types, 
        helper = helper, fitted = fitted)
    result$predict <- function(newdata = NULL) {
        values = list()
        fitted = matrix(NA, nrow = nrow(newdata), ncol = np)
        for (i in 1:np) {
            fitted[, i] = predict(boost.object[[i]], newdata = newdata)
        }
        for (k in 1:length(blsstr)) {
            values[[k]] = matrix(NA, nrow = nrow(newdata), ncol = np)
            for (i in 1:np) {
                values[[k]][, i] = predict(boost.object[[i]], 
                  which = k, newdata = newdata) + boost.object[[i]]$offset
            }
        }
        list(fitted = fitted, values = values)
    }
    class(result) = c("expectreg", "boost")
    result
}
