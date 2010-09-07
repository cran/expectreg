egamma <-
function (e, shape, rate = 1, scale = 1/rate, log.p = FALSE) 
{
    mapply(function(y) optimize(function(z) abs((gamma(1 + rate)/gamma(rate)/shape * 
        pgamma(z, shape, 1 + rate) - z * pgamma(z, shape, rate, 
        scale, log.p))/(2 * (gamma(1 + rate)/gamma(rate)/shape * 
        pgamma(z, shape, 1 + rate) - z * pgamma(z, shape, rate, 
        scale, log.p)) + z - shape/rate) - y), c(0, 1.2 * qgamma(0.995, 
        shape, rate, scale)))$minimum, e)
}
