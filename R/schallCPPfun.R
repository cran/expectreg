schallCPPfun <-
function(glatterms, y, B, tau, lambdashort_in, DD, NB, NBP, NBPC, center) {
    .Call('_expectreg_schallCPPfun', PACKAGE = 'expectreg', glatterms, y, B, tau, lambdashort_in, DD, NB, NBP, NBPC, center)
}
