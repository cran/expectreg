elnorm <-
function (e, meanlog = 0, sdlog = 1, log.p = FALSE) 
{
    mapply(function(y) optimize(function(z) abs((exp(meanlog + 
        0.5 * sdlog^2) * (1 - pnorm((meanlog + sdlog^2 - log(z))/sdlog, 
        log.p = log.p)) - z * plnorm(z, meanlog, sdlog, log.p = log.p))/(2 * 
        (exp(meanlog + 0.5 * sdlog^2) * (1 - pnorm((meanlog + 
            sdlog^2 - log(z))/sdlog, log.p = log.p)) - z * plnorm(z, 
            meanlog, sdlog, log.p = log.p)) + z - exp(meanlog + 
        0.5 * sdlog^2)) - y), c(0, 1.2 * qlnorm(0.995, meanlog, 
        sdlog)))$minimum, e)
}
