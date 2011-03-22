plot.expectreg <-
function (x, ...) 
{
    yy = x$response
    cov = x$covariates
    Z = x$values
    coefficients = x$coefficients
    formula = x$formula
    intercept = x$intercepts
    m = length(yy)
    types = x$effects
    helper = x$helper
    pp = x$expectiles
    np <- length(pp)
    if (identical(pp, seq(0.01, 0.99, by = 0.01))) {
        pp.plot <- c(1, 2, 5, 10, 20, 50, 80, 90, 95, 98, 99)
        row.grid = 3
        col.grid = 4
    }
    else if (identical(pp, c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 
        0.8, 0.9, 0.95, 0.98, 0.99))) {
        pp.plot <- 1:length(pp)
        row.grid = 3
        col.grid = 4
    }
    else {
        if (np > 15) 
            pp.plot = seq(1, np, length = 15)
        else pp.plot <- 1:length(pp)
        row.grid = floor(sqrt(length(pp)))
        col.grid = ceiling(sqrt(length(pp)))
        if (length(pp) > row.grid * col.grid) 
            row.grid = row.grid + 1
    }
    np.plot <- length(pp.plot)
    for (k in 1:length(types)) {
        dev.new()
        if (types[[k]] == "pspline") {
            ZZZ = Z[[k]][order(cov[[k]])[seq(1, m, length = min(m, 
                100))], pp.plot]
            plot(cov[[k]], yy, cex = 0.5, pch = 20, col = "grey42", 
                xlab = "x", ylab = "y", ylim = range(cbind(yy, 
                  Z[[k]])))
            matlines(sort(cov[[k]])[seq(1, m, length = min(m, 
                100))], ZZZ, col = rainbow(np.plot + 1)[1:np.plot], 
                lty = 1)
            legend(x = "topright", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                1)[1:np.plot]), legend = rev(pp[pp.plot]), bg = "white", 
                bty = "n")
        }
        else if (types[[k]] == "markov") {
            z = NULL
            Zspathelp = helper[[k]][[2]]
            bnd = helper[[k]][[1]]
            for (i in 1:np) {
                z = cbind(z, coefficients[[k]][, i] + intercept[i])
            }
            if (inherits(x, "sheets") || inherits(x, "boost")) {
                if (class(bnd) != "bnd") {
                  plot(seq(0, 1.1 * max(cov[[k]]), length = 10), 
                    seq(0, max(z), length = 10), type = "n", 
                    xlab = "Districts", ylab = "coefficients")
                  matpoints(cov[[k]], Z[[k]], col = rainbow(np + 
                    1)[1:np])
                  legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                    1)[1:np.plot]), legend = rev(pp[pp.plot]), 
                    bg = "white", bty = "n")
                }
                else {
                  par(mfrow = (c(row.grid, col.grid)))
                  plot.limits = range(z[[k]])
                  n = as.numeric(attr(bnd, "regions"))
                  for (i in 1:np.plot) {
                    re = data.frame(cbind(cov[[k]], Z[[k]][, 
                      i]))
                    drawmap(re, bnd, regionvar = 1, plotvar = 2, 
                      mar.min = NULL, limits = plot.limits, main = pp[pp.plot[i]], 
                      cols = "grey", swapcolors = TRUE)
                  }
                }
            }
            else {
                if (class(bnd) != "bnd") {
                  plot(seq(0, 1.1 * max(cov[[k]]), length = 10), 
                    seq(0, max(z[, pp.plot]), length = 10), type = "n", 
                    xlab = "District", ylab = "coefficients")
                  points(rep(as.numeric(attr(bnd, "regions")), 
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
                    re = data.frame(cbind(as.numeric(attr(bnd, 
                      "regions")), z[, pp.plot[i]]))
                    drawmap(re, bnd, regionvar = 1, plotvar = 2, 
                      mar.min = NULL, limits = plot.limits, main = pp[pp.plot[i]])
                  }
                }
            }
        }
        else if (types[[k]] == "2dspline") {
            if (inherits(x, "sheets") || inherits(x, "boost")) {
                par(mfrow = (c(row.grid, col.grid)))
                for (i in 1:np) {
                  if (i %in% pp.plot) {
                    z = interp(cov[[k]][, 1], cov[[k]][, 2], 
                      Z[[k]][, i])
                    persp(z[[1]], z[[2]], z[[3]], ticktype = "detailed", 
                      phi = 40, zlim = range(yy), col = "lightblue", 
                      xlab = "X", ylab = "Y", zlab = "Z", main = pp.plot[i]/100)
                  }
                }
            }
            else {
                x.min = apply(cov[[k]], 2, min)
                x.max = apply(cov[[k]], 2, max)
                x.gitter = cbind(rep(seq(x.min[1], x.max[1], 
                  length = 50), times = 50), rep(seq(x.min[2], 
                  x.max[2], length = 50), each = 50))
                B.gitter = base(x.gitter, "2dspline")[[1]]
                par(mfrow = (c(row.grid, col.grid)))
                for (i in 1:np) {
                  if (i %in% pp.plot) {
                    z <- B.gitter %*% coefficients[[k]][, i] + 
                      intercept[i]
                    z = t(matrix(z, nrow = 50, ncol = 50))
                    persp(seq(x.min[1], x.max[1], length = 50), 
                      seq(x.min[2], x.max[2], length = 50), z, 
                      ticktype = "detailed", phi = 40, zlim = range(yy), 
                      col = "lightblue", xlab = "X", ylab = "Y", 
                      zlab = "Z", main = pp.plot[i]/100)
                  }
                }
            }
        }
        else if (types[[k]] == "radial") {
            if (inherits(x, "sheets") || inherits(x, "boost")) {
                par(mfrow = (c(row.grid, col.grid)))
                for (i in 1:np) {
                  if (i %in% pp.plot) {
                    z = interp(cov[[k]][, 1], cov[[k]][, 2], 
                      Z[[k]][, i])
                    persp(z[[1]], z[[2]], z[[3]], ticktype = "detailed", 
                      phi = 40, zlim = range(yy), col = "lightblue", 
                      xlab = "X", ylab = "Y", zlab = "Z", main = pp.plot[i]/100)
                  }
                }
            }
            else {
                x.min = apply(cov[[k]], 2, min)
                x.max = apply(cov[[k]], 2, max)
                x.gitter = cbind(rep(seq(x.min[1], x.max[1], 
                  length = 50), times = 50), rep(seq(x.min[2], 
                  x.max[2], length = 50), each = 50))
                cov[[k]] = cov[[k]][order(cov[[k]][, 1]), ]
                knots = cov[[k]][seq(1, dim(cov[[k]])[1], length = min(50, 
                  dim(cov[[k]])[1])), ]
                B.gitter = matrix(NA, nrow = dim(x.gitter)[1], 
                  ncol = dim(knots)[1])
                for (j in 1:dim(knots)[1]) {
                  r = sqrt(rowSums((x.gitter - matrix(unlist(knots[j, 
                    ]), nrow = nrow(x.gitter), ncol = ncol(knots), 
                    byrow = T))^2))
                  r[r == 0] = 1
                  B.gitter[, j] = r^2 * log(r)
                }
                par(mfrow = (c(row.grid, col.grid)))
                for (i in 1:np) {
                  if (i %in% pp.plot) {
                    z <- B.gitter %*% coefficients[[k]][, i] + 
                      intercept[i]
                    z = t(matrix(z, nrow = 50, ncol = 50))
                    persp(seq(x.min[1], x.max[1], length = 50), 
                      seq(x.min[2], x.max[2], length = 50), z, 
                      ticktype = "detailed", phi = 40, zlim = range(yy), 
                      col = "lightblue", xlab = "X", ylab = "Y", 
                      zlab = "Z", main = pp.plot[i]/100)
                  }
                }
            }
        }
        else if (types[[k]] == "krig") {
            if (inherits(x, "sheets") || inherits(x, "boost")) {
                par(mfrow = (c(row.grid, col.grid)))
                for (i in 1:np) {
                  if (i %in% pp.plot) {
                    z = interp(cov[[k]][, 1], cov[[k]][, 2], 
                      Z[[k]][, i])
                    persp(z[[1]], z[[2]], z[[3]], ticktype = "detailed", 
                      phi = 40, zlim = range(yy), col = "lightblue", 
                      xlab = "X", ylab = "Y", zlab = "Z", main = pp.plot[i]/100)
                  }
                }
            }
            else {
                krig.phi = helper[[k]]
                x.min = apply(cov[[k]], 2, min)
                x.max = apply(cov[[k]], 2, max)
                x.gitter = cbind(rep(seq(x.min[1], x.max[1], 
                  length = 50), times = 50), rep(seq(x.min[2], 
                  x.max[2], length = 50), each = 50))
                cov[[k]] = cov[[k]][order(cov[[k]][, 1]), ]
                knots = cov[[k]][seq(1, dim(cov[[k]])[1], length = min(50, 
                  dim(cov[[k]])[1])), ]
                B.gitter = matrix(NA, nrow = dim(x.gitter)[1], 
                  ncol = dim(knots)[1])
                for (j in 1:dim(knots)[1]) {
                  r = sqrt(rowSums((x.gitter - matrix(unlist(knots[j, 
                    ]), nrow = nrow(x.gitter), ncol = ncol(knots), 
                    byrow = T))^2))/krig.phi
                  B.gitter[, j] = exp(-r) * (1 + r)
                }
                par(mfrow = (c(row.grid, col.grid)))
                for (i in 1:np) {
                  if (i %in% pp.plot) {
                    z <- B.gitter %*% coefficients[[k]][, i] + 
                      intercept[i]
                    z = t(matrix(z, nrow = 50, ncol = 50))
                    persp(seq(x.min[1], x.max[1], length = 50), 
                      seq(x.min[2], x.max[2], length = 50), z, 
                      ticktype = "detailed", phi = 40, zlim = range(yy), 
                      col = "lightblue", xlab = "X", ylab = "Y", 
                      zlab = "Z", main = pp.plot[i]/100)
                  }
                }
            }
        }
        else if (types[[k]] == "random") {
            if (inherits(x, "sheets") || inherits(x, "boost")) {
                plot(seq(0, 1.1 * max(cov[[k]]), length = 10), 
                  seq(0, max(z), length = 10), type = "n", xlab = "Group", 
                  ylab = "coefficients")
                points(rep(sort(unique(cov[[k]])), times = np), 
                  Z[[k]], col = rainbow(np + 1)[1:np])
                legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np + 
                  1)[1:np]), legend = rev(pp), bg = "white", 
                  bty = "n")
            }
            else {
                plot(seq(0, 1.1 * max(cov[[k]]), length = 10), 
                  seq(0, max(coefficients[[k]] + intercept), 
                    length = 10), type = "n", xlab = "Group", 
                  ylab = "coefficients")
                points(rep(sort(unique(cov[[k]])), times = np.plot), 
                  (coefficients[[k]] + intercept)[, pp.plot], 
                  col = rainbow(np.plot + 1)[1:np.plot])
                legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                  1)[1:np.plot]), legend = rev(pp[pp.plot]), 
                  bg = "white", bty = "n")
            }
        }
        else if (types[[k]] == "ridge") {
            plot(seq(0, 1.1 * dim(cov[[k]])[2], length = 10), 
                seq(0, max(coefficients[[k]] + intercept), length = 10), 
                type = "n", xlab = "X variables", ylab = "coefficients")
            points(rep(1:dim(cov[[k]])[2], times = np.plot), 
                (coefficients[[k]] + intercept)[, pp.plot], col = rainbow(np.plot + 
                  1)[1:np.plot])
            legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                1)[1:np.plot]), legend = rev(pp[pp.plot]), bg = "white", 
                bty = "n")
        }
        else if (types[[k]] == "parametric") {
            if (inherits(x, "sheets") || inherits(x, "boost")) {
                plot(cov[[k]], yy, cex = 0.5, pch = 20, col = "grey42", 
                  xlab = "x", ylab = "y", ylim = range(cbind(yy, 
                    Z[[k]])))
                matlines(sort(cov[[k]]), Z[[k]][order(cov[[k]]), 
                  pp.plot], col = rainbow(np.plot + 1)[1:np.plot], 
                  lty = 1)
                legend(x = "bottomright", pch = 19, cex = 1, 
                  col = rev(rainbow(np.plot + 1)[1:np.plot]), 
                  legend = rev(pp[pp.plot]), bg = "white", bty = "n")
            }
            else {
                matplot(1:nrow(coefficients[[k]]), coefficients[[k]][, 
                  pp.plot], col = rainbow(np.plot + 1)[1:np.plot], 
                  pch = 15)
                legend(x = "bottomright", pch = 19, cex = 1, 
                  col = rev(rainbow(np.plot + 1)[1:np.plot]), 
                  legend = rev(pp[pp.plot]), bg = "white", bty = "n")
            }
        }
        else if (types[[k]] == "special") {
            plot(cov[[k]], yy, cex = 0.5, pch = 20, col = "grey42", 
                xlab = "x", ylab = "y", ylim = range(cbind(yy, 
                  Z[[k]])))
            matlines(sort(cov[[k]]), Z[[k]][order(cov[[k]]), 
                pp.plot], col = rainbow(np.plot + 1)[1:np.plot], 
                lty = 1)
            legend(x = "bottomright", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                1)[1:np.plot]), legend = rev(pp[pp.plot]), bg = "white", 
                bty = "n")
        }
    }
}
