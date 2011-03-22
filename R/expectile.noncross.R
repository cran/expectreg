expectile.noncross <-
function (formula, data = NULL, smooth = c("schall", "acv", "fixed"), 
    lambda = 0.1, expectiles = NA) 
{
    smooth = match.arg(smooth)
    if (any(is.na(expectiles)) || !is.vector(expectiles) || any(expectiles > 
        1) || any(expectiles < 0)) {
        p <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 
            0.98, 0.99)
        row.grid = 3
        col.grid = 4
    }
    else {
        p <- expectiles
        row.grid = floor(sqrt(length(p)))
        col.grid = ceiling(sqrt(length(p)))
        if (length(p) > row.grid * col.grid) 
            row.grid = row.grid + 1
    }
    y <- as.numeric(eval(parse(text = formula[2])))
    n <- length(y)
    mp <- length(p)
    types = list()
    design = list()
    if (labels(terms(formula))[1] == "1") {
        x <- c(rep(1, n))
        Bx <- as.matrix(rep(1, n))
        lambda <- 0
        P <- as.matrix(1)
        types[[1]] = "intercept"
    }
    else for (i in 1:length(labels(terms(formula)))) {
        types[[i]] = strsplit(labels(terms(formula))[i], "(", 
            fixed = TRUE)[[1]][1]
        if (types[[i]] == labels(terms(formula))[i]) {
            design[[i]] = base(matrix(eval(parse(text = labels(terms(formula))[i]), 
                envir = data, enclos = environment(formula)), 
                nrow = n), "parametric")
            formula = eval(substitute(update(formula, . ~ variable2 + 
                . - variable1), list(variable1 = as.name(types[[i]]), 
                variable2 = as.name(paste("base(", types[[i]], 
                  ",'parametric',center=FALSE)", sep = "")))))
            types[[i]] = "parametric"
            lambda = 0
        }
        else {
            if (sub("center=TRUE", "center=FALSE", labels(terms(formula))[i]) == 
                labels(terms(formula))[i]) {
                basetext = paste("base(center=FALSE,", strsplit(labels(terms(formula))[i], 
                  "base(", fixed = TRUE)[[1]][2], sep = "")
            }
            else basetext = sub("center=TRUE", "center=FALSE", 
                labels(terms(formula))[i])
            design[[i]] = eval(parse(text = basetext), envir = data, 
                enclos = environment(formula))
        }
    }
    Bx = design[[1]][[1]]
    types[[1]] = design[[1]][[4]]
    if (types[[1]] == "pspline") 
        P = t(design[[1]][[2]]) %*% design[[1]][[2]]
    else P = design[[1]][[2]]
    x = design[[1]][[3]]
    if (types[[1]] == "markov") 
        helper = list(design[[1]][[5]], design[[1]][[6]])
    else if (types[[1]] == "krig") 
        helper = list(design[[1]][[7]])
    else helper = list(NA)
    bdegp <- 1
    p.knots = p
    p.knots[1] <- p.knots[1] - 1e-07
    p.knots[mp] <- p.knots[mp] + 1e-07
    dx.left <- ((p.knots[1] + 1/(mp + 1))/bdegp)
    outer.left <- seq(p.knots[1] - bdegp * dx.left, p.knots[1] - 
        dx.left, by = dx.left)
    dx.right <- (((1 + 1/(mp + 1)) - p.knots[mp])/bdegp)
    outer.right <- seq(p.knots[mp] + dx.right, p.knots[mp] + 
        bdegp * dx.right, by = dx.right)
    knots <- c(outer.left, p.knots, outer.right)
    Bp <- splineDesign(knots, p.knots, bdegp + 1)
    B <- Bp %x% Bx
    yy <- rep(y, mp)
    W0 <- matrix(rep(0.5), ncol = ncol(B), nrow = nrow(B), byrow = FALSE)
    W0B <- W0 * B
    Dmat <- t(B) %*% W0B
    dvec <- t(yy) %*% W0B
    Amat <- matrix(0, ncol = ncol(B), nrow = ncol(Bx) * (nrow(Bp) - 
        1))
    bvec <- matrix(0, ncol = ncol(Bx) * (nrow(Bp) - 1), 1)
    a_p <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat), 
        bvec = bvec)
    a_p <- as.vector(unlist(a_p$solution))
    W_temp <- W0
    z <- B %*% a_p
    pw <- c()
    for (k in 1:mp) {
        pw <- c(pw, rep(p[k], n))
    }
    Amat <- diag(nrow = ncol(Bx) * (nrow(Bp) - 1))
    Amat <- cbind(matrix(0, ncol = ncol(Bx), nrow = ncol(Bx) * 
        (nrow(Bp) - 1)), Amat)
    diag(Amat) <- -1
    bvec <- matrix(0, ncol = ncol(Bx) * (nrow(Bp) - 1), 1)
    if (smooth == "acv") {
        acv.min = nlm(acv.noncross, p = lambda, yy = yy, Bx = Bx, 
            Bp = Bp, P = P, Amat = Amat, bvec = bvec, a_p = a_p, 
            pp = p, ndigit = 8, iterlim = 50, gradtol = 1e-04)
        lambda <- abs(acv.min$estimate)
    }
    difw <- 1
    difl <- 1
    it <- 0
    while ((difw > 0 || difl > 1e-05) && (it < 100)) {
        lambda0 <- lambda
        W_temp[, 1:(ncol(Bx) * ncol(Bp))] <- pw * (yy > z) + 
            (1 - pw) * (yy <= z)
        W_tempB <- W_temp * B
        Dmat <- t(B) %*% W_tempB + 2 * lambda * (diag(1, mp - 
            1 + bdegp, mp - 1 + bdegp) %x% P)
        dvec <- t(yy) %*% W_tempB
        a_pq <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat), 
            bvec = bvec)
        a_sol <- as.vector(unlist(a_pq$solution))
        z <- B %*% a_sol
        W1 <- W_temp
        W_temp[, 1:(ncol(Bx) * ncol(Bp))] <- pw * (yy > z) + 
            (1 - pw) * (yy <= z)
        difw <- sum(abs(W_temp - W1))
        if (smooth == "schall") {
            BW1B <- t(B) %*% (W1 * B)
            H <- solve(BW1B + 2 * lambda * (diag(1, mp - 1 + 
                bdegp, mp - 1 + bdegp) %x% P)) %*% BW1B
            df <- sum(diag(H))
            sig2 <- as.numeric((t(yy - z) %*% (W1[, 1] * (yy - 
                z)))/(length(yy) - df))
            lambda <- as.numeric(df/((t(a_sol) %*% (diag(1, mp - 
                1 + bdegp, mp - 1 + bdegp) %x% P) %*% a_sol) + 
                1e-06)) * sig2
        }
        it <- it + 1
        difl <- abs(lambda0 - lambda)
        print(it)
    }
    ncexpect <- list(lambda = list(lambda), intercepts = rep(0, 
        length(p)), values = list(matrix(z, ncol = length(p))), 
        coefficients = list(matrix(a_sol, ncol = length(p))), 
        response = y, formula = formula, expectiles = p, effects = types, 
        helper = helper, covariates = list(x), design = Bx, fitted = matrix(z, 
            ncol = length(p)))
    class(ncexpect) = c("expectreg", "noncross")
    ncexpect
}
