dkoenker <-
function (z) 
{
    mapply(function(x) ifelse(x == 0, 0.25, 2 * abs(x)/sqrt(1 - 
        4/(4 + x^2))/(4 + x^2)^2), z)
}
