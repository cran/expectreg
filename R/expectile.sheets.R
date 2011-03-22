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
    B = list()
    DD = list()
    center = TRUE
    varying = list()
    for (i in 1:length(labels(terms(formula)))) {
        types[[i]] = strsplit(labels(terms(formula))[i], "(", 
            fixed = TRUE)[[1]][1]
        if (types[[i]] == labels(terms(formula))[i]) {
            design[[i]] = base(types[[i]], "parametric")
            types[[i]] = "parametric"
        }
        else design[[i]] = eval(parse(text = labels(terms(formula))[i]), 
            envir = data, enclos = environment(formula))
    }
    nterms = length(design)
    B[[1]] = design[[1]][[1]]
    DD[[1]] = design[[1]][[2]]
    x[[1]] = design[[1]][[3]]
    types[[1]] = design[[1]][[4]]
    bnd[[1]] = design[[1]][[5]]
    Zspathelp[[1]] = design[[1]][[6]]
    nb[1] = ncol(design[[1]][[1]])
    krig.phi[[1]] = design[[1]][[7]]
    center = center && design[[1]][[8]]
    varying[[1]] = design[[1]][[9]]
    if (length(labels(terms(formula))) > 1) 
        for (i in 2:length(labels(terms(formula)))) {
            B[[i]] = design[[i]][[1]]
            DD[[i]] = design[[i]][[2]]
            x[[i]] = design[[i]][[3]]
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
    vector.a.ma.schall <- matrix(NA, nrow = sum(nb) + 1, ncol = np)
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
        vector.a.ma.schall <- matrix(p2f.new$coef, ncol = np, 
            byrow = T)
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
    coefficients <- NA
    final.lambdas <- list()
    helper <- list()
    intercept = p2f.new$intercept
    desmat = 1
    for (k in 1:length(design)) {
        final.lambdas[[k]] = lala[k, ]
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
    }
    result = list(lambda = final.lambdas, intercepts = intercept, 
        coefficients = NULL, values = curves, response = yy, 
        covariates = x, formula = formula, expectiles = pp, effects = types, 
        helper = helper, design = desmat, fitted = p2f.new$fit)
    class(result) = c("expectreg", "sheets")
    result
}
