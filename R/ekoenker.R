ekoenker <-
function (q) 
{
    mapply(function(y) {
        if (y == 0) 
            return(-Inf)
        else if (y == 1) 
            return(Inf)
        else nlm(function(z) abs((integrate(function(x) x * dkoenker(x), 
            -Inf, z)$value - z * pkoenker(z))/(2 * (integrate(function(x) x * 
            dkoenker(x), -Inf, z)$value - z * pkoenker(z)) + 
            z) - y), p = 0)$estimate
    }, q)
}
