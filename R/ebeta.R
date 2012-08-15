ebeta <-
function (asy, a = 1, b = 1) 
{
    zz = 0 * asy
    for (k in 1:length(asy)) {
        root = function(z) pebeta(z, a, b) - asy[k]
        z = uniroot(root, interval = c(0, 1), tol = 1e-06)
        zz[k] = z$root
    }
    return(zz)
}
