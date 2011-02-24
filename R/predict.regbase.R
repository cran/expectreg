predict.regbase <-
function (object, newdata = NULL, ...) 
{
    type = object$type
    bnd = object$bnd
    Zspathelp = object$Zspathelp
    phi = object$phi
    center = object$center
    x = object$x
    if (is.null(newdata)) 
        newdata = x
    if (type == "pspline") {
        require(splines)
        B.deg = 2
        B.size = 20
        diff.size = 2
        x0 <- min(x) - 0.001
        x1 <- max(x) + 0.001
        dx = (x1 - x0)/(B.size - 1)
        B = splineDesign(knots = seq(x0 - dx * B.deg, x1 + dx * 
            B.deg, by = dx), x = newdata, ord = B.deg + 1)
        P <- diag(dim(B)[2])
        P <- diff(P, diff = diff.size)
        if (center) {
            tildeU <- matrix(0, dim(B)[2], diff.size)
            for (i in 1:diff.size) tildeU[, i] <- (1:(dim(B)[2]))^(i - 
                1)
            tildeZ <- t(P) %*% solve(P %*% t(P))
            U <- B %*% tildeU
            Z <- B %*% tildeZ
            B = cbind(U[, -1], Z)
            P <- diag(c(rep(0, ncol(U) - 1), rep(1, ncol(Z))))
        }
    }
    else if (type == "2dspline") {
        require(splines)
        B.deg = 2
        B.size = 20
        diff.size = 2
        x0 <- min(x[, 1]) - 0.001
        x1 <- max(x[, 1]) + 0.001
        dx = (x1 - x0)/(B.size - B.deg)
        Bx = splineDesign(knots = seq(x0 - dx * B.deg, x1 + dx * 
            B.deg, by = dx), x = newdata[, 1], ord = B.deg + 
            1)
        y0 <- min(x[, 2]) - 0.001
        y1 <- max(x[, 2]) + 0.001
        dy = (y1 - y0)/(B.size - B.deg)
        By = splineDesign(knots = seq(y0 - dy * B.deg, y1 + dy * 
            B.deg, by = dy), x = newdata[, 2], ord = B.deg + 
            1)
        B = matrix(NA, nrow = dim(x)[1], ncol = dim(Bx)[2] * 
            dim(By)[2])
        for (i in 1:dim(Bx)[1]) B[i, ] = as.vector(Bx[i, ] %o% 
            By[i, ])
        D = diff(diag(dim(Bx)[2]), diff = diff.size)
        P = t(D) %*% D
        P = diag(dim(Bx)[2]) %x% P + P %x% diag(dim(Bx)[2])
        if (center) {
            tildeX = matrix(1, nrow = dim(Bx)[2] * dim(By)[2], 
                ncol = 4)
            tildeX[, 2] = rep(1:dim(Bx)[2], times = dim(By)[2])
            tildeX[, 3] = rep(1:dim(By)[2], each = dim(Bx)[2])
            tildeX[, 4] = tildeX[, 2] * tildeX[, 3]
            ek = eigen(P)
            ek$values = ek$values[1:(length(ek$values) - 2 * 
                diff.size)]
            ek$vectors = ek$vectors[, 1:(dim(ek$vectors)[2] - 
                2 * diff.size)]
            L = ek$vectors %*% sqrt(diag(ek$values))
            tildeU = L %*% solve(t(L) %*% L)
            X = B %*% tildeX
            U = B %*% tildeU
            B = cbind(X[, -1], U)
            P <- diag(c(rep(0, ncol(X) - 1), rep(1, ncol(U))))
        }
    }
    else if (type == "radial") {
        x = x[order(x[, 1]), ]
        knots = x[seq(1, dim(x)[1], length = min(50, dim(x)[1])), 
            ]
        B = matrix(NA, nrow = dim(newdata)[1], ncol = dim(knots)[1])
        for (i in 1:dim(newdata)[1]) for (j in 1:dim(knots)[1]) {
            r = sqrt(sum((newdata[i, ] - knots[j, ])^2))
            if (r == 0) 
                B[i, j] = 0
            else B[i, j] = r^2 * log(r)
        }
        P = matrix(0, nrow = dim(B)[2], ncol = dim(B)[2])
        for (i in 1:dim(B)[2]) for (j in 1:dim(B)[2]) {
            r = sqrt(sum((knots[i, ] - knots[j, ])^2))
            if (r == 0) 
                P[i, j] = 0
            else P[i, j] = r^2 * log(r)
        }
    }
    else if (type == "krig") {
        c = 9.233
        x = x[order(x[, 1]), ]
        knots = x[seq(1, dim(x)[1], length = min(50, dim(x)[1])), 
            ]
        B = matrix(NA, nrow = dim(newdata)[1], ncol = dim(knots)[1])
        P = matrix(0, nrow = dim(B)[2], ncol = dim(B)[2])
        for (i in 1:dim(B)[2]) for (j in 1:dim(B)[2]) {
            P[i, j] = sqrt(sum((knots[i, ] - knots[j, ])^2))
        }
        phi = max(P)/c
        P = P/phi
        P = exp(-P) * (1 + P)
        for (i in 1:dim(newdata)[1]) for (j in 1:dim(knots)[1]) {
            r = sqrt(sum((newdata[i, ] - knots[j, ])^2))
            B[i, j] = exp(-r/phi) * (1 + r/phi)
        }
    }
    else if (type == "markov") {
        require(BayesX)
        if (any(!is.na(bnd)) && any(is.na(P))) 
            P = bnd2gra(bnd)
        if (all(is.na(P))) 
            stop("No Neighbourhood defined.")
        if (any(diag(P) < 0) || any(P - diag(P) > 0) || is.null(dimnames(P)[[1]])) 
            stop("Maldefined Neighbourhood.")
        districts = dimnames(P)[[1]]
        B = matrix(0, nrow = length(newdata), ncol = dim(P)[2])
        for (i in 1:length(newdata)) B[i, which(districts == 
            newdata[i])] = 1
        if (center) {
            e <- eigen(P)
            L <- e$vectors[, -dim(e$vectors)[2]] * sqrt(e$values[-length(e$values)])
            Zspathelp <- L %*% solve(t(L) %*% L)
            B <- B %*% Zspathelp
            P = diag(nrow = dim(P)[1] - 1)
        }
    }
    else if (type == "random") {
        districts = sort(unique(x))
        B = matrix(0, nrow = length(newdata), ncol = length(districts))
        for (i in 1:length(newdata)) B[i, which(districts == 
            newdata[i])] = 1
        P = diag(nrow = dim(B)[2])
    }
    else if (type == "ridge") {
        B = matrix(NA, nrow = dim(newdata)[1], ncol = dim(x)[2])
        for (i in 1:dim(x)[2]) B[, i] = x[, i]
        P = diag(nrow = dim(B)[2])
    }
    else if (type == "special") {
        if (is.na(B) || is.na(P)) 
            stop("In 'special' case: Base and Penalty matrix have to be specified.")
        if (is.vector(newdata)) 
            l = length(newdata)
        else if (is.matrix(newdata) || is.data.frame(newdata)) 
            l = nrow(newdata)
        B = as.matrix(object$B)[1:l, ]
        P = as.matrix(P)
    }
    else if (type == "parametric") {
        if (is.vector(x)) {
            B = matrix(newdata, ncol = 1)
        }
        else if (is.matrix(x)) {
            B = newdata
        }
        else if (is.data.frame(x)) {
            B = as.matrix(newdata)
        }
        if (center) 
            B = apply(B, 2, function(x) {
                x - sum(x)/length(x)
            })
        P = matrix(0, nrow = ncol(B), ncol = ncol(B))
    }
    B
}