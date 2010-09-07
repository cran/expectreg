enorm <-
function (e, m = 0, sd = 1, log.p = FALSE) 
{
    mapply(function(y) optimize(function(z) abs((m * pnorm((z - 
        m)/sd, log.p = log.p) - sd * dnorm((z - m)/sd, log = log.p) - 
        z * pnorm(z, m = m, sd = sd, log.p = log.p))/(2 * (m * 
        pnorm((z - m)/sd, log.p = log.p) - sd * dnorm((z - m)/sd, 
        log = log.p) - z * pnorm(z, m = m, sd = sd, log.p = log.p)) + 
        z - m) - y), c(-10, 10) * sd + m)$minimum, e)
}
