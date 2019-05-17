penorm <-
function (q, m = 1, sd = 2) 
{
    z = (q - m)/sd
    p = pnorm(z)
    d = dnorm(z)
    u = -d - z * p
    asy = u/(3 * u + z)
    return(asy)
}
