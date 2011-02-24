echisq <-
function (e, df, log.p = FALSE) 
{
    egamma(e, shape = df/2, rate = 0.5, log.p = log.p)
}
