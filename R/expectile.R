expectile <-
function (x, probs = seq(0, 1, 0.25), dec = 4) 
{
    if (!is.vector(x)) 
        stop("observations are needed in vector form.")
    if (min(probs) < 0 || max(probs) > 1) 
        stop("only asymmetries between 0 and 1 allowed.")
    res = NULL
    for (i in 1:length(probs)) {
        if (probs[i] == 1) 
            res[i] = max(x, na.rm = T)
        else res[i] = asyregpen.lsfit(x, matrix(1, nrow = length(x), 
            ncol = 1), probs[i], NULL, matrix(0, nrow = 1, ncol = 1), 
            0)$a
    }
    names(res) = probs
    round(res, dec)
}
