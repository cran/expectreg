elnorm <-
function (asy, meanlog = 0, sdlog = 1) 
{
    zz = 0 * asy
    for (k in 1:length(asy)) {
        root = function(z) pelnorm(z, meanlog, sdlog) - asy[k]
        z = uniroot(root, interval = c(0, 1.5 * qlnorm(0.995, 
            meanlog, sdlog)), tol = 1e-06)
        zz[k] = z$root
    }
    return(zz)
}
