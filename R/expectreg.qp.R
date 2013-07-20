expectreg.qp <-
function (formula, data = NULL, smooth = c("schall", "acv", "fixed"), 
    lambda = 1, expectiles = NA) 
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
    y = eval(as.expression(formula[[2]]), envir = data, enclos = environment(formula))
    attr(y, "name") = deparse(formula[[2]])
    n <- length(y)
    mp <- length(p)
    i.ja <- 0
    if (attr(terms(formula), "intercept") == "1") 
        i.ja <- 1
    types = list()
    design = list()
    helper = list()
    x = list()
    if (attr(terms(formula), "intercept") == "1") {
        design[[1]] = rb(matrix(1, nrow = n, ncol = 1), "parametric", 
            center = F)
        types[[1]] = "intercept"
    }
    else if (formula[[3]] == ".") {
        design[[1]] = rb(data[, names(data) != all.vars(formula[[2]])], 
            "parametric")
    }
    if (length(labels(terms(formula))) > 0) 
        for (i in 1:length(labels(terms(formula)))) {
            types[[i.ja + i]] = strsplit(strsplit(labels(terms(formula))[i], 
                "(", fixed = TRUE)[[1]][2], ",", fixed = TRUE)[[1]][2]
            if (is.na(types[[i.ja + i]])) 
                types[[i.ja + i]] = "pspline"
            if (types[[i.ja + i]] == "pspline" || types[[i.ja + 
                i]] == " \"pspline\")") 
                types[[i.ja + i]] <- "pspline"
            if (types[[i.ja + i]] == " \"parametric\")" || types[[i.ja + 
                i]] == "parametric") {
                lambda[i.ja + i] <- 0
                P <- diag(0, 2, 2)
                types[[i.ja + i]] <- "parametric"
                x.mom <- eval(parse(text = labels(terms(formula))[i], 
                  srcfile = NULL), envir = data, enclos = environment(formula))$x
                knots <- c(min(x.mom) - 2, min(x.mom) - 1e-04, 
                  max(x.mom) + 1e-04, max(x.mom) + 2)
                B <- splineDesign(knots = knots, x = x.mom, ord = 2)
                design[[i.ja + i]] <- list(B = B, P = P, x = x.mom, 
                  type = types[[i.ja + i]], constraint = matrix(0, 
                    nrow = 2, ncol = ncol(B)))
                design[[i.ja + i]]$xname <- eval(parse(text = labels(terms(formula))[i], 
                  srcfile = NULL), envir = data, enclos = environment(formula))$xname
            }
            else if (strsplit(labels(terms(formula))[i], "(", 
                fixed = TRUE)[[1]][1] != "mono") {
                if (sub("center=TRUE", "center=FALSE", labels(terms(formula))[i]) == 
                  labels(terms(formula))[i]) {
                  basetext = paste("rb(center=FALSE,", strsplit(labels(terms(formula))[i], 
                    "rb(", fixed = TRUE)[[1]][2], sep = "")
                }
                else {
                  basetext = sub("center=TRUE", "center=FALSE", 
                    labels(terms(formula))[i])
                }
                design[[i.ja + i]] = eval(parse(text = basetext), 
                  envir = data, enclos = environment(formula))
            }
            else design[[i.ja + i]] = eval(parse(text = labels(terms(formula))[i]), 
                envir = data, enclos = environment(formula))
        }
    Bx = design[[1]][[1]]
    K <- vector(length = length(design))
    K[1] <- ncol(Bx)
    P = design[[1]][[2]]
    P = t(P) %*% P
    P.list <- list()
    P.list[[1]] <- P
    bigP <- P
    x[[1]] = design[[1]][[3]]
    names(x)[1] = design[[1]]$xname[1]
    if (types[[1]] == "markov") 
        helper[[1]] = list(design[[1]][[5]], design[[1]][[6]])
    else if (types[[1]] == "krig") 
        helper[[1]] = list(design[[1]][[7]])
    else helper[[1]] = list(NA)
    Cmat = list()
    Cmat[[1]] = as.matrix(design[[1]]$constraint)
    for (i in 2:length(p)) {
        Cmat[[1]] = rbind(cbind(Cmat[[1]], matrix(0, nrow = nrow(Cmat[[1]]), 
            ncol = ncol(as.matrix(design[[1]]$constraint)))), 
            cbind(matrix(0, nrow = nrow(as.matrix(design[[1]]$constraint)), 
                ncol = ncol(Cmat[[1]])), as.matrix(design[[1]]$constraint)))
    }
    if (length(design) > 1) {
        for (i in 2:length(design)) {
            Bx = cbind(Bx, design[[i]][[1]])
            K[i] <- ncol(design[[i]][[1]])
            types[[i]] = design[[i]][[4]]
            bigP = rbind(cbind(P, matrix(0, nrow = nrow(P), ncol = ncol(design[[i]][[2]]))), 
                cbind(matrix(0, nrow = nrow(design[[i]][[2]]), 
                  ncol = ncol(P)), design[[i]][[2]]))
            P.list[[i]] <- t(design[[i]][[2]]) %*% design[[i]][[2]]
            x[[i]] = design[[i]][[3]]
            names(x)[i] = design[[i]]$xname[1]
            if (types[[i]] == "markov") 
                helper[[i]] = list(design[[i]][[5]], design[[i]][[6]])
            else if (types[[i]] == "krig") 
                helper[[i]] = list(design[[i]][[7]])
            else helper[[i]] = list(NA)
            Cmat[[i]] = as.matrix(design[[i]]$constraint)
            for (k in 2:length(p)) {
                Cmat[[i]] = rbind(cbind(Cmat[[i]], matrix(0, 
                  nrow = nrow(Cmat[[i]]), ncol = ncol(as.matrix(design[[i]]$constraint)))), 
                  cbind(matrix(0, nrow = nrow(as.matrix(design[[i]]$constraint)), 
                    ncol = ncol(Cmat[[i]])), as.matrix(design[[i]]$constraint)))
            }
        }
    }
    if (length(lambda) < length(design)) {
        lambda <- rep(lambda[1], length(design))
    }
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
    dfparam <- vector(length = length(types))
    for (i in 1:length(types)) {
        dfparam[i] <- (types[[i]] == "parametric")
    }
    dfparam <- sum(dfparam) + i.ja
    bigBxnp <- matrix(0, ncol = 0, nrow = n)
    bigPnp <- matrix(0, ncol = 0, nrow = (21 - 2))
    lambdanp <- rep(0, times = length(lambda))
    indexnp <- vector(length = 0)
    for (i in 1:length(types)) {
        if (types[[i]] == "pspline") {
            bigBxnp <- cbind(bigBxnp, design[[i]][[1]])
            bigPnp <- rbind(cbind(bigPnp, matrix(0, nrow = nrow(bigPnp), 
                ncol = ncol(design[[i]][[2]]))), cbind(matrix(0, 
                nrow = nrow(design[[i]][[2]]), ncol = ncol(bigPnp)), 
                design[[i]][[2]]))
            lambdanp[i] <- lambda[i]
            indexnp <- c(indexnp, i)
        }
    }
    partbasis.list = list()
    for (k in 1:length(design)) {
        partbasis.list[[k]] = (sum(K[0:(k - 1)]) + 1):(sum(K[0:k]))
    }
    acv.noncross <- function(penalty, yy, Bx, Bp, pw, bdegp, 
        n, mp, P.list, design, partbasis.list) {
        aa <- noncross.backfit(yy, Bx, Bp, pw, bdegp, n, mp, 
            abs(penalty), P.list, design, partbasis.list, Cmat = Cmat)
        score = (aa$weights * (yy - aa$yyhat)^2)/((1 - sum(aa$dfi)/(n * 
            mp))^2)
        mean(score[which(is.finite(score))], na.rm = TRUE)
    }
    noncross.fit <- function(yy, Bx, Bp, bdegp, mp, lambda, P, 
        weights, Cmat, types, yywh, i) {
        B <- Bp %x% Bx
        W0 <- matrix(weights, nrow = nrow(B), ncol = ncol(B), 
            byrow = FALSE)
        W0B <- W0 * B
        Amat <- B[-(1:n), ] - B[-((nrow(B) - n + 1):nrow(B)), 
            ]
        if (types == "intercept") {
            Amat <- matrix(0, ncol = ncol(Amat), nrow = nrow(Amat))
            meq <- 0
        }
        else meq <- mp
        if (types != "intercept" && !i.ja && i == 1) 
            meq <- 0
        if (ncol(Amat) == ncol(Cmat)) {
            Amat <- rbind(Amat, Cmat)
            bvec <- matrix(0, ncol = n * (mp - 1) + meq + nrow(Cmat), 
                1)
        }
        else bvec <- matrix(0, ncol = n * (mp - 1) + meq, 1)
        if ((types != "intercept" && i.ja) || (types != "intercept" && 
            !i.ja && i != 1)) {
            Amat_pC_help <- matrix(0, ncol = n * mp, nrow = mp)
            for (m in 1:mp) Amat_pC_help[m, ((m - 1) * n + 1):(m * 
                n)] <- 1
            Amat_partCenter <- Amat_pC_help %*% B
            Amat <- rbind(Amat_partCenter, Amat)
            bvec_partNC <- yywh[-((nrow(B) - n + 1):nrow(B))] - 
                yywh[-(1:n)]
            bvec[, (meq + 1):(n * (mp - 1) + meq)] <- bvec_partNC
        }
        if (types != "intercept" && !i.ja && i == 1) {
            bvec_partNC <- yywh[-((nrow(B) - n + 1):nrow(B))] - 
                yywh[-(1:n)]
            bvec[, (meq + 1):(n * (mp - 1) + meq)] <- bvec_partNC
            print("hier!")
        }
        Dmat <- t(B) %*% W0B + 2 * lambda * (diag(1, mp - 1 + 
            bdegp, mp - 1 + bdegp) %x% P)
        dvec <- t(yy) %*% W0B
        a_pq <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat), 
            bvec = bvec, meq = meq)
        a_sol <- as.vector(unlist(a_pq$solution))
        z <- B %*% a_sol
        tau2dfi <- as.numeric((t(a_sol) %*% (diag(1, mp - 1 + 
            bdegp, mp - 1 + bdegp) %x% P) %*% a_sol))
        dfi <- sum(diag(solve(Dmat) %*% t(B) %*% W0B))
        ncexpect <- list(lambda = lambda, fitted = as.vector(z), 
            coefficients = a_sol, weights = W0[, 1], tau2dfi = tau2dfi, 
            Amat = Amat, bvec = bvec, dfi = dfi)
    }
    pw <- c()
    for (k in 1:mp) {
        pw <- c(pw, rep(p[k], n))
    }
    yy <- rep(y, mp)
    noncross.backfit <- function(yy, Bx, Bp, pw, bdegp, n, mp, 
        lambdanp, P.list, design, partbasis.list, Cmat) {
        weights <- rep(0.5, times = n * mp)
        f <- matrix(0, nrow = n * mp, ncol = length(design))
        alpha <- mean(yy)
        ok <- TRUE
        rss0 <- 0
        coefficients <- list()
        coeff.vec <- list()
        Z <- list()
        part.resid <- list()
        dfi <- vector(length = length(design))
        tau2dfi <- vector(length = length(design))
        while (ok) {
            for (i in 1:length(design)) {
                difw <- 1
                it <- 0
                while ((difw > 0) && (it < 150)) {
                  weights0 <- weights
                  epsilon <- yy - rowSums(as.matrix(f[, -i])) - 
                    alpha
                  yywi <- rowSums(as.matrix(f[, -i])) + alpha
                  nc.obj <- noncross.fit(yy = epsilon, Bx = as.matrix(Bx[, 
                    partbasis.list[[i]]]), Bp = Bp, mp = mp, 
                    bdegp = bdegp, lambda = lambdanp[i], P = as.matrix(P.list[[i]]), 
                    weights = weights, Cmat[[i]], types = types[[i]], 
                    yywh = yywi, i = i)
                  coeff.vec[[i]] <- nc.obj$coefficients
                  coefficients[[i]] <- matrix(nc.obj$coefficients, 
                    ncol = length(p))
                  f[, i] <- nc.obj$fitted
                  print(apply(f, 2, sum))
                  yyhat <- alpha + rowSums(f)
                  weights <- pw * (yy > yyhat) + (1 - pw) * (yy <= 
                    yyhat)
                  it <- it + 1
                  difw <- sum(abs(weights0 - weights))
                  if (it == 150) 
                    warning("Weights did not converge. Stopping after 150 iterations.")
                }
                part.resid[[i]] <- matrix(epsilon, nrow = n, 
                  ncol = length(p))
                Z[[i]] <- matrix(f[, i], nrow = n, ncol = length(p))
                dfi[i] <- nc.obj$dfi
                if (types[[i]] != "intercept") 
                  dfi[i] <- dfi[i] - 1
                tau2dfi[i] <- nc.obj$tau2dfi
            }
            rss <- sum((yy - yyhat)^2)
            if (abs(rss - rss0) < 1e-06 * rss) 
                ok <- FALSE
            rss0 <- rss
        }
        ncexpectbf <- list(lambdanp = lambdanp, fitted = f, coefficients = coefficients, 
            weights = weights, values = Z, alpha = alpha, yyhat = yyhat, 
            dfi = dfi, tau2dfi = tau2dfi, part.resid = part.resid, 
            coeff.vec = coeff.vec)
    }
    if (smooth == "schall") {
        dc <- 1
        it <- 1
        while (dc >= 1e-04 && it < 100) {
            nc.bf <- noncross.backfit(yy = yy, Bx = Bx, Bp = Bp, 
                pw = pw, bdegp = bdegp, n = n, mp = mp, lambdanp = lambdanp, 
                P.list = P.list, design = design, partbasis.list = partbasis.list, 
                Cmat = Cmat)
            tau2 <- (nc.bf$tau2dfi/nc.bf$dfi) + 1e-06
            sig2df <- as.numeric(t(yy - nc.bf$yyhat) %*% (nc.bf$weights * 
                (yy - nc.bf$yyhat)))
            df <- sum(nc.bf$dfi)
            sig2 <- as.numeric(sig2df/(length(yy) - df))
            lambdanp <- sig2/tau2
            for (i in 1:length(design)) {
                if ((types[[i]] == "parametric") || (types[[i]] == 
                  "intercept")) 
                  lambdanp[i] <- 0
            }
            dc <- sum((nc.bf$lambdanp - lambdanp)^2)
            it <- it + 1
            if (it == 100) 
                warning("Schall algorithm did not converge. Stopping after 100 iterations.")
        }
    }
    else if (smooth == "acv") {
        acv.min <- nlm(acv.noncross, p = lambdanp, yy = yy, Bx = Bx, 
            Bp = Bp, pw = pw, bdegp = bdegp, n = n, mp = mp, 
            P = P.list, design = design, partbasis.list = partbasis.list, 
            ndigit = 8, iterlim = 50, gradtol = 1e-04)
        lambdanp <- abs(acv.min$estimate)
        nc.bf <- noncross.backfit(yy = yy, Bx = Bx, Bp = Bp, 
            pw = pw, bdegp = bdegp, n = n, mp = mp, lambdanp = lambdanp, 
            P.list = P.list, design = design, partbasis.list = partbasis.list, 
            Cmat = Cmat)
    }
    else {
        nc.bf <- noncross.backfit(yy = yy, Bx = Bx, Bp = Bp, 
            pw = pw, bdegp = bdegp, n = n, mp = mp, lambdanp = lambdanp, 
            P.list = P.list, design = design, partbasis.list = partbasis.list, 
            Cmat = Cmat)
    }
    nc.bf$fitted <- nc.bf$fitted + nc.bf$alpha
    add_alpha <- function(x) x <- x + nc.bf$alpha
    Z <- lapply(nc.bf$values, add_alpha)
    alpha <- nc.bf$alpha
    if (attr(terms(formula), "intercept") == "1") {
        nc.bf$fitted[, -1] <- nc.bf$fitted[, -1] + nc.bf$fitted[, 
            1]
        add_intercept <- function(x) x <- x + nc.bf$fitted[, 
            1]
        Z <- lapply(nc.bf$values, add_intercept)
        Z[[1]] <- nc.bf$values[[1]]
    }
    coefficients <- nc.bf$coefficients
    coeff.vec <- nc.bf$coeff.vec
    dfi <- nc.bf$dfi
    final.lambdas <- nc.bf$lambdanp
    final.lambdas = list(final.lambdas)
    if (i.ja == 1) {
        intercepts <- Z[[1]][1, ]
        Z <- Z[-1]
        types <- types[-1]
        x <- x[-1]
    }
    else intercepts <- rep(0, length(p))
    fitted <- matrix(NA, ncol = length(p), nrow = n * (length(design) - 
        i.ja))
    for (j in 1:(length(design) - i.ja)) {
        fitted[((j - 1) * n + 1):(j * n), ] <- Z[[j]]
        names(Z)[j] = design[[j]]$xname[1]
        names(coefficients)[j] = design[[j]]$xname[1]
    }
    yyhat = matrix(nc.bf$yyhat, nrow = n, ncol = length(p))
    ncexpect <- list(lambda = final.lambdas, intercepts = intercepts, 
        values = Z, coefficients = coefficients, response = y, 
        formula = formula, asymmetries = p, effects = types, 
        helper = helper, covariates = x, design = Bx, bases = design, 
        fitted = yyhat, weights = nc.bf$weights, part.resid = nc.bf$part.resid, 
        dfi = dfi, coeff.vec = coeff.vec, alpha = alpha)
    ncexpect$predict <- function(newdata = NULL) {
        BB = list()
        values = list()
        bmat = NULL
        betavec = NULL
        fitted <- matrix(NA, ncol = length(p), nrow = length(newdata))
        if (attr(terms(formula), "intercept") == "1") 
            it = 2
        else it = 1
        for (k in it:length(coefficients)) {
            BB[[k]] = predict(design[[k]], newdata)
            values[[k]] <- BB[[k]] %*% coefficients[[k]]
            bmat = cbind(bmat, BB[[k]])
            betavec = rbind(betavec, coefficients[[k]])
        }
        add_alpha <- function(x) x <- x + nc.bf$alpha
        values <- lapply(values, add_alpha)
        names(values) = names(coefficients)
        if (attr(terms(formula), "intercept") == "1") {
            bmat = cbind(1, bmat)
            betavec = rbind(intercepts, betavec)
        }
        fitted = bmat %*% betavec + alpha
        names(values) = names(coefficients)
        list(fitted = fitted, values = values)
    }
    class(ncexpect) = c("expectreg", "noncross")
    ncexpect
}
