et <-
function (asy, df) 
{
    zz = 0 * asy
    for (k in 1:length(asy)) {
        root = function(z) pet(z, df) - asy[k]
        z = uniroot(root, interval = c(-10, 10), tol = 1e-06)
        zz[k] = z$root
    }
    return(zz)
}
