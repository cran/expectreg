eunif <-
function (e, min = 0, max = 1) 
{
    mapply(function(y) optimize(function(z) abs(((z^2 - min^2)/2/(max - 
        min) - z * punif(z, min = min, max = max))/(2 * ((z^2 - 
        min^2)/2/(max - min) - z * punif(z, min = min, max = max)) + 
        z - (min + max)/2) - y), c(min, max))$minimum, e)
}
