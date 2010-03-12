eunif <-
function (e, min = 0, max = 1, lower.tail = TRUE, log.p = FALSE) 
{
    mapply(function(y) optimize(function(z) abs((integrate(function(x) x * 
        dunif(x, min = min, max = max, log = log.p), min, z)$value - 
        z * punif(z, min = min, max = max, lower.tail = lower.tail, 
            log.p = log.p))/(2 * (integrate(function(x) x * dunif(x, 
        min = min, max = max, log = log.p), min, z)$value - z * 
        punif(z, min = min, max = max, lower.tail = lower.tail, 
            log.p = log.p)) + z - (min + max)/2) - y), c(min, 
        max))$minimum, e)
}
