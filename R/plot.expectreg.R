plot.expectreg <-
function (x, ...) 
{
    if (inherits(x, "boost")) 
        stop("Plotting function for laws, restricted and bundle.")
    expectreg = x
    yy = expectreg$response
    x = expectreg$covariates
    Z = expectreg$values
    coefficients = expectreg$coefficients
    formula = expectreg$formula
    intercept = expectreg$intercepts
    m = length(yy)
    if (ncol(coefficients[[1]]) != 11) 
        pp <- seq(0.01, 0.99, by = 0.01)
    else pp <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 
        0.98, 0.99)
    np <- length(pp)
    design <- list()
    types <- list()
    Zspathelp <- list()
    bnd <- list()
    krig.phi = list()
    nb = vector()
    for (i in 1:length(labels(terms(formula)))) {
        design[[i]] = eval(parse(text = labels(terms(formula))[i]))
        nb[i] = ncol(design[[i]][[1]])
        types[[i]] = design[[i]][[4]]
        bnd[[i]] = design[[i]][[5]]
        Zspathelp[[i]] = design[[i]][[6]]
        krig.phi[[i]] = design[[i]][[7]]
    }
    for (k in 1:length(design)) {
        partbasis = (sum(nb[0:(k - 1)]) + 1):(sum(nb[0:k]))
        dev.new()
        if (types[[k]] == "pspline") {
            plot(x[[k]], yy, cex = 0.5, pch = 20, col = "grey42", 
                xlab = "x", ylab = "y", ylim = range(cbind(yy, 
                  Z[[k]])))
            matlines(sort(x[[k]])[seq(1, m, length = min(m, 100))], 
                Z[[k]][seq(1, m, length = min(m, 100)), ], col = rainbow(np + 
                  1)[1:11], lty = 1)
            legend(x = "bottomright", pch = 19, cex = 1, col = rev(rainbow(np + 
                1)[1:np]), legend = rev(pp), bg = "white", bty = "n")
        }
        else if (types[[k]] == "markov") {
            z = NULL
            for (i in 1:np) {
                z = cbind(z, diag(dim(Zspathelp[[k]])[1]) %*% 
                  coefficients[[k]] + intercept[i])
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
            x.min = apply(x[[k]], 2, min)
            x.max = apply(x[[k]], 2, max)
            x.gitter = cbind(rep(seq(x.min[1], x.max[1], length = 50), 
                times = 50), rep(seq(x.min[2], x.max[2], length = 50), 
                each = 50))
            B.gitter = base(x.gitter, "2dspline")[[1]]
            par(mfrow = (c(3, 4)))
            for (i in 1:np) {
                z <- B.gitter %*% coefficients[[k]][, i] + intercept[i]
                z = t(matrix(z, nrow = 50, ncol = 50))
                persp(seq(x.min[1], x.max[1], length = 50), seq(x.min[2], 
                  x.max[2], length = 50), z, ticktype = "detailed", 
                  phi = 40, zlim = range(yy), col = "lightblue", 
                  xlab = "X", ylab = "Y", zlab = "Z", main = pp[i])
            }
        }
        else if (types[[k]] == "radial") {
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
                z <- B.gitter %*% coefficients[[k]][, i] + intercept[i]
                z = t(matrix(z, nrow = 50, ncol = 50))
                persp(seq(x.min[1], x.max[1], length = 50), seq(x.min[2], 
                  x.max[2], length = 50), z, ticktype = "detailed", 
                  phi = 40, zlim = range(yy), col = "lightblue", 
                  xlab = "X", ylab = "Y", zlab = "Z", main = pp[i])
            }
        }
        else if (types[[k]] == "krig") {
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
                z <- B.gitter %*% coefficients[[k]][, i] + intercept[i]
                z = t(matrix(z, nrow = 50, ncol = 50))
                persp(seq(x.min[1], x.max[1], length = 50), seq(x.min[2], 
                  x.max[2], length = 50), z, ticktype = "detailed", 
                  phi = 40, zlim = range(yy), col = "lightblue", 
                  xlab = "X", ylab = "Y", zlab = "Z", main = pp[i])
            }
        }
        else if (types[[k]] == "random") {
            plot(seq(0, 1.1 * max(x[[k]]), length = 10), seq(0, 
                max(coefficients[[k]] + intercept), length = 10), 
                type = "n", xlab = "Group", ylab = "coefficients")
            points(rep(sort(unique(x[[k]])), times = np), coefficients[[k]] + 
                intercept, col = rainbow(np + 1)[1:np])
            legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np + 
                1)[1:np]), legend = rev(pp), bg = "white", bty = "n")
        }
        else if (types[[k]] == "ridge") {
            plot(seq(0, 1.1 * dim(x[[k]])[2], length = 10), seq(0, 
                max(coefficients[[k]] + intercept), length = 10), 
                type = "n", xlab = "X variables", ylab = "coefficients")
            points(rep(1:dim(x[[k]])[2], times = np), coefficients[[k]] + 
                intercept, col = rainbow(np + 1)[1:np])
            legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np + 
                1)[1:np]), legend = rev(pp), bg = "white", bty = "n")
        }
    }
}
