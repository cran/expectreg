echisq <-
function (e, df, log.p = FALSE) 
{
    egamma(e, 0.5, df/2, log.p = log.p)
}
