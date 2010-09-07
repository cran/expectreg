elnorm <-
function (e, meanlog = 0, sdlog = 1, log.p = FALSE) 
{
    mapply(function(y) optimize(function(z) abs((exp(meanlog + 
        0.5 * sdlog^2) * pnorm((log(z) - meanlog - sdlog^2)/sdlog, 
        log.p = log.p) - z * plnorm(z, meanlog, sdlog, log.p))/(2 * 
        (exp(meanlog + 0.5 * sdlog^2) * pnorm((log(z) - meanlog - 
            sdlog^2)/sdlog, log.p = log.p) - z * plnorm(z, meanlog, 
            sdlog, log.p)) + z - exp(meanlog + 0.5 * sdlog^2)) - 
        y), c(0, 15 * sdlog + meanlog))$minimum, e)
}
