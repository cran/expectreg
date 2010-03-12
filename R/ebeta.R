ebeta <-
function (e, a = 1, b = 1) 
{
    mapply(function(y) optimize(function(z) abs((integrate(function(x) x * 
        dbeta(x, a, b), 0, z)$value - z * pbeta(z, a, b))/(2 * 
        (integrate(function(x) x * dbeta(x, a, b), 0, z)$value - 
            z * pbeta(z, a, b)) + z - a/(a + b)) - y), c(0, 1))$minimum, 
        e)
}
