ebeta <-
function (e, a = 1, b = 1, log.p = FALSE) 
{
    mapply(function(y) optimize(function(z) abs((gamma(1 + a) * 
        gamma(a + b)/gamma(a)/gamma(1 + a + b) * pbeta(z, 1 + 
        a, b, log.p = log.p) - z * pbeta(z, a, b, log.p = log.p))/(2 * 
        (gamma(1 + a) * gamma(a + b)/gamma(a)/gamma(1 + a + b) * 
            pbeta(z, 1 + a, b, log.p = log.p) - z * pbeta(z, 
            a, b, log.p = log.p)) + z - a/(a + b)) - y), c(0, 
        1))$minimum, e)
}
