quant.boost <-
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
    blsstr <- labels(terms(formula))
    types = list()
    bls <- list()
    x <- list()
    bnd = list()
    yy = eval(parse(text = formula[2]), envir = data, enclos = .GlobalEnv)
    if (is.null(data)) 
        data = yy
    for (i in 1:length(blsstr)) {
        types[[i]] = strsplit(blsstr[i], "(", fixed = TRUE)[[1]][1]
        if (types[[i]] == blsstr[i]) {
            types[[i]] = "parametric"
            x[[i]] = eval(parse(text = blsstr[i]), envir = data, 
                enclos = .GlobalEnv)
            bls[[i]] = NULL
        }
        else {
            bls[[i]] <- eval(parse(text = blsstr[i]), envir = data, 
                enclos = .GlobalEnv)
            x[[i]] <- as.matrix(bls[[i]]$get_data())
        }
        if (types[[i]] != "bmrf") 
            bnd[[i]] = NA
        else bnd[[i]] = environment(bls[[i]]$dpp)$args$bnd
        data = cbind(data, x[[i]])
    }
    data = data.frame(data)
    m = length(yy)
    kfld <- 10
    ntest <- floor(m/kfld)
    cv10f <- matrix(c(rep(c(rep(0, ntest), rep(1, m)), kfld - 
        1), rep(0, m * kfld - (kfld - 1) * (m + ntest))), nrow = m)
    values <- list()
    z <- list()
    for (k in 1:length(blsstr)) {
        values[[k]] = matrix(NA, nrow = m, ncol = np)
        if (types[[k]] == "bbs") 
            z[[k]] = values[[k]]
        else if (types[[k]] == "bmrf") 
            z[[k]] = matrix(NA, nrow = length(attr(bnd[[k]], 
                "regions")), ncol = np)
        else if (types[[k]] == "brandom") 
            z[[k]] = matrix(NA, nrow = length(unique(x[[k]])), 
                ncol = np)
        else if (types[[k]] == "bspatial") 
            z[[k]] = list()
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
            nu = 0.1, risk = "inbag"), family = QuantReg(pp[p]))
        if (cv) {
            cvr = cvrisk(inb, folds = cv10f)
            print(paste("Quantile", pp[p], ", Boosting iterations:", 
                mstop(cvr), sep = " "))
            inb = inb[mstop(cvr)]
        }
        for (k in 1:length(blsstr)) {
            if (types[[k]] == "bbs") {
                independent = which(dimnames(data)[[2]] == bls[[k]]$get_names())
                d.tmp = data
                for (j in (1:dim(d.tmp)[2])[-independent]) if (is.factor(d.tmp[, 
                  j])) 
                  d.tmp[, j] = as.numeric(names(which.max(table(d.tmp[, 
                    j]))))
                else d.tmp[, j] = mean(d.tmp[, j])
                values[[k]] = predict(inb, d.tmp)
                z[[k]] = values[[k]]
            }
            else if (types[[k]] == "bspatial") {
                x[[k]] = data[, bls[[k]]$get_names()]
                x.min = apply(x[[k]], 2, min)
                x.max = apply(x[[k]], 2, max)
                x.gitter = cbind(rep(seq(x.min[1], x.max[1], 
                  length = 50), times = 50), rep(seq(x.min[2], 
                  x.max[2], length = 50), each = 50))
                d.tmp = matrix(0, nrow = 2500, ncol = dim(data)[2])
                for (j in (1:dim(d.tmp)[2])[-independent]) if (is.factor(d.tmp[, 
                  j])) 
                  d.tmp[, j] = as.numeric(names(which.max(table(d.tmp[, 
                    j]))))
                else d.tmp[, j] = mean(d.tmp[, j])
                d.val = data
                independent = c(which(dimnames(data)[[2]] == 
                  bls[[k]]$get_names()[1]), which(dimnames(data)[[2]] == 
                  bls[[k]]$get_names()[2]))
                d.tmp[, independent] = x.gitter
                dimnames(d.tmp) = list(1:2500, dimnames(data)[[2]])
                for (j in 1:dim(d.val)[2][-independent]) if (is.factor(d.val[, 
                  j])) 
                  d.val[, j] = as.numeric(names(which.max(table(d.val[, 
                    j]))))
                else d.val[, j] = mean(d.val[, j])
                z[[k]] <- predict(inb, data.frame(d.tmp))
                values[[k]] = predict(inb, d.val)
                z[[k]] = t(matrix(z[[k]], nrow = 50, ncol = 50))
            }
            else if (types[[k]] == "brandom") {
                districts = unique(x[[k]])
                d.tmp = matrix(0, nrow = length(districts), ncol = dim(data)[2])
                dimnames(d.tmp) = list(1:length(districts), dimnames(data)[[2]])
                independent = which(dimnames(data)[[2]] == bls[[k]]$get_names())
                d.tmp[, independent] = districts
                for (j in (1:dim(d.tmp)[2])[-independent]) if (is.factor(d.tmp[, 
                  j])) 
                  d.tmp[, j] = as.numeric(names(which.max(table(d.tmp[, 
                    j]))))
                else d.tmp[, j] = mean(d.tmp[, j])
                d.val = data
                for (j in 1:dim(d.val)[2][-independent]) if (is.factor(d.val[, 
                  j])) 
                  d.val[, j] = as.numeric(names(which.max(table(d.val[, 
                    j]))))
                else d.val[, j] = mean(d.val[, j])
                values[[k]] = predict(inb, d.val)
                z[[k]] <- predict(inb, d.tmp)
            }
            else if (types[[k]] == "bmrf") {
                districts = as.numeric(attr(bnd[[k]], "regions"))
                d.tmp = data.frame(matrix(0, nrow = length(districts), 
                  ncol = dim(data)[2]))
                dimnames(d.tmp) = list(1:length(districts), dimnames(data)[[2]])
                independent = which(dimnames(data)[[2]] == bls[[k]]$get_names())
                d.tmp[, independent] = districts
                d.tmp[, independent] = as.factor(d.tmp[, independent])
                for (j in (1:dim(d.tmp)[2])[-independent]) if (is.factor(d.tmp[, 
                  j])) 
                  d.tmp[, j] = as.numeric(names(which.max(table(d.tmp[, 
                    j]))))
                else d.tmp[, j] = mean(d.tmp[, j])
                d.val = data
                for (j in 1:dim(d.val)[2][-independent]) if (is.factor(d.val[, 
                  j])) 
                  d.val[, j] = as.numeric(names(which.max(table(d.val[, 
                    j]))))
                else d.val[, j] = mean(d.val[, j])
                values[[k]] = predict(inb, d.val)
                z[[k]] <- predict(inb, d.tmp)
            }
            else if (types[[k]] == "parametric") {
                independent = which(dimnames(data)[[2]] == blsstr[k])
                d.tmp = data
                for (j in (1:dim(d.tmp)[2])[-independent]) if (is.factor(d.tmp[, 
                  j])) 
                  d.tmp[, j] = as.numeric(names(which.max(table(d.tmp[, 
                    j]))))
                else d.tmp[, j] = mean(d.tmp[, j])
                values[[k]] = predict(inb, d.tmp)
                z[[k]] = values[[k]]
            }
        }
        rm(inb)
        gc()
        list(values, z)
    }
    coef.vector = myapply(1:np, function(i) dummy.reg(i, formula, 
        data, mstop, pp, cv10f, types, x, blsstr, bnd))
    for (i in 1:np) for (k in 1:length(blsstr)) {
        values[[k]][, i] = coef.vector[[i]][[1]][[k]]
        if (types[[k]] == "bspatial") 
            z[[k]][[i]] = coef.vector[[i]][[2]][[k]]
        else z[[k]][, i] = coef.vector[[i]][[2]][[k]]
    }
    for (k in 1:length(blsstr)) {
        dev.new()
        if (types[[k]] == "bbs") {
            plot(x[[k]], yy, cex = 0.5, pch = 20, col = "grey42", 
                xlab = "x", ylab = "y", ylim = range(cbind(yy, 
                  values[[k]])))
            matlines(sort(x[[k]]), z[[k]][order(x[[k]]), ], col = rainbow(np + 
                1)[1:np], lty = 1)
            legend(x = "bottomright", pch = 19, cex = 1, col = rev(rainbow(np + 
                1)[1:np]), legend = rev(pp), bg = "white", bty = "n")
        }
        else if (types[[k]] == "bspatial") {
            x.min = apply(x[[k]], 2, min)
            x.max = apply(x[[k]], 2, max)
            par(mfrow = (c(3, 4)))
            for (i in 1:np) {
                persp(seq(x.min[1], x.max[1], length = 50), seq(x.min[2], 
                  x.max[2], length = 50), z[[k]][[i]], ticktype = "detailed", 
                  phi = 40, zlim = range(yy), col = "lightblue", 
                  xlab = "X", ylab = "Y", zlab = "Z", main = pp[i])
            }
        }
        else if (types[[k]] == "brandom") {
            plot(seq(0, 1.1 * max(x[[k]]), length = 10), seq(0, 
                max(z), length = 10), type = "n", xlab = "Group", 
                ylab = "coefficients")
            points(rep(sort(unique(x[[k]])), times = np), z[[k]], 
                col = rainbow(np + 1)[1:np])
            legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np + 
                1)[1:np]), legend = rev(pp), bg = "white", bty = "n")
        }
        else if (types[[k]] == "bmrf") {
            if (class(bnd[[k]]) != "bnd") {
                plot(seq(0, 1.1 * max(x[[k]]), length = 10), 
                  seq(0, max(z), length = 10), type = "n", xlab = "Districts", 
                  ylab = "coefficients")
                points(rep(as.numeric(attr(bnd[[k]], "regions")), 
                  times = np), z[[k]], col = rainbow(np + 1)[1:np])
                legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np + 
                  1)[1:np]), legend = rev(pp), bg = "white", 
                  bty = "n")
            }
            else {
                par(mfrow = (c(3, 4)))
                plot.limits = range(z[[k]])
                n = as.numeric(attr(bnd[[k]], "regions"))
                for (i in 1:np) {
                  re = data.frame(cbind(n, z[[k]][, i]))
                  drawmap(re, bnd[[k]], regionvar = 1, plotvar = 2, 
                    mar.min = NULL, limits = plot.limits, main = pp[i])
                }
            }
        }
        else if (types[[k]] == "parametric") {
            plot(x[[k]], yy, cex = 0.5, pch = 20, col = "grey42", 
                xlab = "x", ylab = "y", ylim = range(cbind(yy, 
                  values[[k]])))
            matlines(sort(x[[k]]), z[[k]][order(x[[k]]), ], col = rainbow(np + 
                1)[1:np], lty = 1)
            legend(x = "bottomright", pch = 19, cex = 1, col = rev(rainbow(np + 
                1)[1:np]), legend = rev(pp), bg = "white", bty = "n")
        }
    }
    result = list(values = values, response = yy, covariates = x, 
        formula = formula, expectiles = pp)
    class(result) = c("expectreg", "boost")
    result
}
