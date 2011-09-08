eunif <-
function (asy, min = 0, max = 1) 
{
    zz = 0 * asy
    for (k in 1:length(asy)) {
        root = function(z) peunif(z, min, max) - asy[k]
        z = uniroot(root, interval = c(min, max), tol = 1e-06)
        zz[k] = z$root
    }
    return(zz)
}
