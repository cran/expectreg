enorm <-
function (asy, m = 0, sd = 1) 
{
    zz = 0 * asy
    for (k in 1:length(asy)) {
        root = function(z) penorm(z) - asy[k]
        z = uniroot(root, interval = c(-10, 10), tol = 1e-06)
        zz[k] = z$root * sd + m
    }
    return(zz)
}
