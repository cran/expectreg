pkoenker <-
function (z) 
{
    mapply(function(x) 0.5 + 0.5 * sign(x) * sqrt(1 - 4/(4 + 
        x^2)), z)
}
