print.expectreg <-
function (x, ...) 
{
    cat("\nCall:\n", deparse(x$formula), "\n", sep = "")
    cat("\nHead of Data:\n", sep = "")
    print(cbind(response = x$response, matrix(unlist(x$covariates), 
        nrow = length(x$response)))[1:min(6, length(x$response)), 
        ])
    cat("\n")
    cat("\nFitted Expectiles:\n", sep = "")
    print(x$expectiles)
    cat("\n")
    if (!inherits(x, "boost")) {
        cat("\nSmoothing Parameters:\n", sep = "")
        print(x$lambda)
        cat("\n")
        cat("\nIntercepts:\n", sep = "")
        print(x$inter)
        cat("\n")
    }
    cat("\nRegression Coefficients:\n", sep = "")
    y = list()
    for (i in 1:length(x$coef)) y[[i]] = x$coef[[i]][1:min(20, 
        nrow(x$coef[[i]])), ]
    print(y)
    cat("\n")
    invisible(x)
}
