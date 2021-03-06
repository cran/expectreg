elnorm <-
function (asy, meanlog = 0, sdlog = 1) 
{
    asy[asy > 1 | asy < 0] = NA
    zz = 0 * asy
    lower = rep(0, length(asy))
    upper = rep(1.5 * qlnorm(0.995, meanlog, sdlog), length(asy))
    diff = 1
    index = 1
    while (diff > 1e-10 && index < 1000) {
        root = pelnorm(zz, meanlog, sdlog) - asy
        root[is.na(root)] = 0
        lower[root < 0] = zz[root < 0]
        upper[root > 0] = zz[root > 0]
        zz = (upper + lower)/2
        diff = max(abs(root), na.rm = T)
        index = index + 1
    }
    zz[is.na(asy)] = NA
    return(zz)
}
