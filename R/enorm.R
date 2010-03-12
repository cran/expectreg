enorm <-
function (e, m = 0, sd = 1) 
{
    mapply(function(y) optimize(function(z) abs((integrate(function(x) x * 
        dnorm(x, m = m, sd = sd), -Inf, z)$value - z * pnorm(z, 
        m = m, sd = sd))/(2 * (integrate(function(x) x * dnorm(x, 
        m = m, sd = sd), -Inf, z)$value - z * pnorm(z, m = m, 
        sd = sd)) + z - m) - y), c(-10, 10) * sd)$minimum, e)
}
