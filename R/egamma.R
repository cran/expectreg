egamma <-
function (asy, shape, rate = 1, scale = 1/rate) 
{
    if (scale != 1/rate) 
        rate = 1/scale
    zz = 0 * asy
    for (k in 1:length(asy)) {
        root = function(z) pegamma(z, shape, rate) - asy[k]
        z = uniroot(root, interval = c(0, 1.5 * qgamma(0.995, 
            shape, rate)), tol = 1e-06)
        zz[k] = z$root
    }
    return(zz)
}
