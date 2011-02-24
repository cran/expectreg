asyregpen.lsfit <-
function (y, B, p, lambda, DD, nb) 
{
    w1 <- 0 * y + 0.5
    n <- ncol(B)
    lambda = c(rep(0, times = n - sum(nb)), rep(lambda, times = nb))
    P <- sqrt(lambda) * DD
    augm <- rep(0, dim(P)[1])
    dw1 = 1
    it = 1
    while (dw1 != 0 && it < 50) {
        model <- lsfit(x = rbind(B, P), y = c(y, augm), wt = c(w1, 
            (augm + 1)), intercept = FALSE)
        a1 <- model$coefficients
        z1 <- B %*% a1
        w01 <- w1
        w1 <- as.vector(ifelse(y > z1, p, 1 - p))
        dw1 <- sum(w1 != w01, na.rm = TRUE)
        it = it + 1
    }
    diag.hat.ma1 <- hat(model$qr)[1:length(y)]
    list(a = a1, diag.hat.ma = diag.hat.ma1, weight = w1, fitted = z1)
}
