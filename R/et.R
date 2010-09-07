et <-
function (e, df, log.p = FALSE) 
{
    mapply(function(y) optimize(function(z) abs(((df + z^2) * 
        dt(z, df, log = log.p)/(1 - df) - z * pt(z, df, log.p = log.p))/(2 * 
        ((df + z^2) * dt(z, df, log = log.p)/(1 - df) - z * pt(z, 
            df, log.p = log.p)) + z - 0) - y), c(-10, 10))$minimum, 
        e)
}
