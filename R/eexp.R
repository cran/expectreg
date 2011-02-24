eexp <-
function (e, rate = 1, log.p = FALSE) 
{
    egamma(e, shape = 1, rate = rate, log.p = log.p)
}
