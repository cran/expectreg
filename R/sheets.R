sheets <-
function (B, DD, yy, pp, lambda, smooth, nb, center) 
{
    nterms = length(nb)
    m = length(yy)
    np = length(pp)
    lala <- matrix(c(rep(lambda, nterms), rep(lambda, nterms)), 
        nrow = nterms, ncol = 2, dimnames = list(1:nterms, c("curve", 
            "sheet")))
    med = which(pp == 0.5)
    if (smooth == "acv") {
        acv.min = nlm(acv.sheets, p = lala, yy = yy, B = B, pp = pp, 
            DD = DD, nb = nb, center = center, ndigit = 8, iterlim = 50, 
            gradtol = 1e-04)
        min.lambda = matrix(abs(acv.min$estimate), ncol = 2)
        lala[, 1] <- min.lambda[, 1]
        lala[, 2] <- min.lambda[, 2]
        ynp <- rep(yy, np)
        ps <- rep(pp, each = m)
        w <- runif(m * np)
        p2f.new <- pspfit2d.new(B, DD, ps, ynp, w, lala[, 1], 
            lala[, 2], center)
    }
    else {
        ynp <- rep(yy, np)
        ps <- rep(pp, each = m)
        w <- runif(m * np)
        p2f.new <- pspfit2d.new(B, DD, ps, ynp, w, lala[, 1], 
            lala[, 2], center)
    }
    curves = p2f.new$curves
    coefficients <- p2f.new$coef
    diag.hat = matrix(p2f.new$hatma, ncol = np)
    result = list(coefficients, lala, diag.hat)
    result
}
