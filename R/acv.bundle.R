acv.bundle <-
function (penalty, yy, B, pp, DD, nb) 
{
    penalty = matrix(abs(penalty), ncol = 2)
    aa <- asyregpen.lsfit(yy, B, 0.5, penalty[, 1], DD, nb)
    residuals = yy - B %*% aa$a
    b <- rep(1, ncol(B))
    cc <- pp - 0.5
    for (i in 1:20) {
        mo <- fitampllsfit(residuals, B, b, pp, cc, DD, penalty[, 
            2], nb)
        b <- mo$b/((sum(mo$b^2)/length(mo$b))^(1/2))
        c0 <- cc
        cc <- fitasy(residuals, B, b, pp, cc)
        dc <- max(abs(cc - c0))
        if (dc < 1e-06) 
            break
    }
    coef.vector = matrix(NA, nrow = sum(nb) + 1, ncol = length(pp))
    for (q in 1:length(pp)) {
        coef.vector[, q] = aa$a + cc[q] * b
    }
    score = 0
    for (i in 1:length(pp)) {
        score = score + mo$weight[, i] * (yy - B %*% coef.vector[, 
            i])^2/(1 - aa$diag.hat.ma)^2
    }
    mean(score[which(is.finite(score))], na.rm = TRUE)
}
