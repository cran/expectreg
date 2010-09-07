eexp <-
function (e, rate = 1, log.p = FALSE) 
{
    egamma(e, rate, 1, 1, log.p)
}
