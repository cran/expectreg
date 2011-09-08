expectile.sheets <-
function (formula, data = NULL, smooth = c("acv", "fixed"), lambda = 0.1, 
    lambdap = 5, expectiles = NA, density = FALSE) 
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
        if (length(expectiles) == 1) 
            warning("Method best applied to more than one expectile.")
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
    attr(yy, "name") = deparse(formula[[2]])
    m = length(yy)
    design = list()
    x = list()
    types = list()
    bnd = list()
    Zspathelp = list()
    nb = vector()
    krig.phi = list()
    B = list()
    DD = list()
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
    B[[1]] = design[[1]][[1]]
    DD[[1]] = design[[1]][[2]]
    x[[1]] = design[[1]][[3]]
    names(x)[1] = design[[1]]$xname
    types[[1]] = design[[1]][[4]]
    bnd[[1]] = design[[1]][[5]]
    Zspathelp[[1]] = design[[1]][[6]]
    nb[1] = ncol(design[[1]][[1]])
    krig.phi[[1]] = design[[1]][[7]]
    center = center && design[[1]][[8]]
    varying[[1]] = design[[1]][[9]]
    if (length(design) > 1) 
        for (i in 2:length(labels(terms(formula)))) {
            B[[i]] = design[[i]][[1]]
            DD[[i]] = design[[i]][[2]]
            x[[i]] = design[[i]][[3]]
            names(x)[i] = design[[i]]$xname
            types[[i]] = design[[i]][[4]]
            bnd[[i]] = design[[i]][[5]]
            Zspathelp[[i]] = design[[i]][[6]]
            nb[i] = ncol(design[[i]][[1]])
            krig.phi[[i]] = design[[i]][[7]]
            center = center && design[[i]][[8]]
            varying[[i]] = design[[i]][[9]]
        }
    lala <- matrix(c(rep(lambda, nterms), rep(lambdap, nterms)), 
        nrow = nterms, ncol = 2, dimnames = list(1:nterms, c("curve", 
            "sheet")))
    med = which(pp == 0.5)
    if (smooth == "acv") {
        acv.min = nlm(acv.sheets, p = lala, yy = yy, B = B, pp = pp, 
            DD = DD, nb = nb, center = center, by = varying, 
            ndigit = 8, iterlim = 50, gradtol = 1e-04)
        min.lambda = matrix(abs(acv.min$estimate), ncol = 2)
        lala[, 1] <- min.lambda[, 1]
        lala[, 2] <- min.lambda[, 2]
        ynp <- rep(yy, np)
        ps <- rep(pp, each = m)
        vv <- log(ps/(1 - ps))
        w <- runif(m * np)
        p2f.new <- pspfit2d.new(B, DD, ps, ynp, w, lala[, 1], 
            lala[, 2], center, varying)
    }
    else {
        ynp <- rep(yy, np)
        ps <- rep(pp, each = m)
        vv <- log(ps/(1 - ps))
        w <- runif(m * np)
        p2f.new <- pspfit2d.new(B, DD, ps, ynp, w, lala[, 1], 
            lala[, 2], center, varying)
    }
    curves = p2f.new$curves
    Z <- list()
    coefficients <- p2f.new$coef
    final.lambdas <- list()
    helper <- list()
    intercept = p2f.new$intercept
    desmat = 1
    for (k in 1:length(design)) {
        final.lambdas[[k]] = lala[k, ]
        names(final.lambdas)[k] = design[[k]]$xname
        desmat = cbind(desmat, B[[k]])
        if (types[[k]] == "markov") {
            coefficients[[k]] = matrix(NA, nrow = nb[k] + 1 * 
                center, ncol = np)
            helper[[k]] = list(bnd[[k]], Zspathelp[[k]])
        }
        else if (types[[k]] == "krig") {
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = krig.phi[[k]]
        }
        else helper[[k]] = NA
        names(curves)[k] = design[[k]]$xname
        names(coefficients)[k] = design[[k]]$xname
    }
    result = list(lambda = final.lambdas, intercepts = intercept, 
        coefficients = coefficients, values = curves, response = yy, 
        covariates = x, formula = formula, expectiles = pp, effects = types, 
        helper = helper, design = desmat, fitted = p2f.new$fit)
    result$predict <- function(newdata = NULL) {
        BB = list()
        values = list()
        for (k in 1:length(types)) {
            BB[[k]] = predict(design[[k]], newdata)
        }
        vv <- log(ps/(1 - ps))
        x0 <- min(vv, -5) - 0.501
        x1 <- max(vv, 5) + 0.501
        B.size = 10
        B.deg = 2
        dx = (x1 - x0)/(B.size - 1)
        By = splineDesign(knots = seq(x0 - dx * B.deg, x1 + dx * 
            B.deg, by = dx), x = vv, ord = B.deg + 1)
        basisX = NULL
        Bx = NULL
        for (k in 1:(length(ynp)/nrow(BB[[1]]))) {
            Bx = rbind(Bx, BB[[1]])
        }
        if (center) {
            Bx = cbind(1, Bx)
        }
        nx <- ncol(Bx)
        ny <- ncol(By)
        B1 <- kronecker(Bx, t(rep(1, ny)))
        B2 <- kronecker(t(rep(1, nx)), By)
        basisX <- B1 * B2
        nb = ncol(basisX) - 1 * center
        if (length(BB) > 1) 
            for (i in 2:length(BB)) {
                Bx = NULL
                for (k in 1:(length(ynp)/nrow(BB[[i]]))) {
                  Bx = rbind(Bx, BB[[i]])
                }
                nx <- ncol(Bx)
                ny <- ncol(By)
                B1 <- kronecker(Bx, t(rep(1, ny)))
                B2 <- kronecker(t(rep(1, nx)), By)
                nb = c(nb, ncol(B1))
                basisX <- cbind(basisX, B1 * B2)
            }
        fitted <- basisX %*% as.vector(coefficients)
        my = nrow(B[[1]])
        np = length(ps)/my
        Bg = basisX[, 1:nb[1]]
        ag = coefficients[1:nb[1]]
        values[[1]] = Bg %*% ag
        dim(values[[1]]) = c(my, np)
        if (length(B) > 1) 
            for (i in 2:length(B)) {
                Bg = NULL
                ag = NULL
                if (center) {
                  Bx = NULL
                  for (k in 1:(length(yy)/nrow(B[[i]]))) {
                    if (any(!is.na(by[[i]]))) 
                      Bx = rbind(Bx, B[[i]] * by[[i]])
                    else Bx = rbind(Bx, B[[i]])
                  }
                  Bx = cbind(1, Bx)
                  nx <- ncol(Bx)
                  ny <- ncol(By)
                  B1 <- kronecker(Bx, t(rep(1, ny)))
                  B2 <- kronecker(t(rep(1, nx)), By)
                  Bg = B1 * B2
                  ag = matrix(coefficients[(sum(nb[1:(i - 1)]) + 
                    1):(sum(nb[1:i]))], nrow = B.size + 1)
                  ag = cbind(intercept, ag)
                  ag = as.vector(ag)
                }
                else {
                  Bg = basisX[, (sum(nb[1:(i - 1)]) + 1):(sum(nb[1:i]))]
                  ag = coefficients[(sum(nb[1:(i - 1)]) + 1):(sum(nb[1:i]))]
                }
                values[[i]] = Bg %*% ag
                dim(values[[i]]) = c(my, np)
            }
        list(fitted = fitted, values = values)
    }
    class(result) = c("expectreg", "sheets")
    result
}
