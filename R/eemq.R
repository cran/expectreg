eemq <-
function (q, ncp = 0, s = 1) 
{
    if (s <= 0) 
        stop("scaling parameter must be strictly >0.")
    mapply(function(y) {
        if (y == 0) 
            return(-Inf)
        else if (y == 1) 
            return(Inf)
        else nlm(function(z) abs((integrate(function(x) x * demq(x, 
            ncp, s), -Inf, z)$value - z * pemq(z, ncp, s))/(2 * 
            (integrate(function(x) x * demq(x, ncp, s), -Inf, 
                z)$value - z * pemq(z, ncp, s)) + z - ncp) - 
            y), p = 0)$estimate
    }, q)
}
