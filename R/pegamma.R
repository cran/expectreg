pegamma <-
function (e, shape, rate = 1, scale = 1/rate) 
{
    u = gamma(1 + shape)/gamma(shape)/rate * pgamma(e, 1 + shape, 
        rate) - e * pgamma(e, shape, rate, scale)
    asy = u/(2 * u + e - shape/rate)
    return(asy)
}
