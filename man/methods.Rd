\name{methods}
\Rdversion{1.1}
\alias{plot.expectreg}
\alias{print.expectreg}
\alias{predict.expectreg}
\alias{residuals.expectreg}
\alias{resid.expectreg}
\alias{coef.expectreg}
\alias{coefficients.expectreg}
\alias{fitted.expectreg}
\alias{fitted.values.expectreg}
\alias{effects.expectreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for expectile regression objects
}
\description{
Methods for objects returned by expectile regression functions.
}
\usage{
\method{plot}{expectreg}(x, ...)

\method{print}{expectreg}(x, ...)

\method{predict}{expectreg}(object, newdata = NULL, ...)

\method{residuals}{expectreg}(object, ...)
\method{resid}{expectreg}(object, ...)

\method{fitted}{expectreg}(object, ...)
\method{fitted.values}{expectreg}(object, ...)

\method{effects}{expectreg}(object, ...)

\method{coef}{expectreg}(object, ...)
\method{coefficients}{expectreg}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{expectreg} as returned e.g. by the function \code{\link{expectile.laws}}.
}
  \item{object}{
An object of class \code{expectreg} as returned e.g. by the function \code{\link{expectile.laws}}.
}
  \item{newdata}{ 
Optionally, a data frame in which to look for variables with which to predict.
}
  \item{\dots}{ additional arguments passed over. }
}
\details{
These functions can be used to extract details from fitted models.
  \code{print} shows a dense representation of the model fit.

The \code{plot} function gives a visual representation of the fitted expectiles
separately for each covariate.

\code{resid} returns the residuals in order of the response.

\code{fitted} returns the overall fitted values \eqn{\hat{y}} while \code{effects} returns the values
for each covariate in a list.

The function \code{coef} extracts the regression coefficients for each covariate listed separately.
For the functions \code{\link{expectile.boost}} and \code{\link{expectile.sheets}} this is not possible.
}
\references{
Schnabel S and Eilers P (2009)
\emph{ Optimal expectile smoothing }
Computational Statistics and Data Analysis, 53:4168-4177

Sobotka F and Kneib T (2010)
\emph{ Geoadditive Expectile Regression }
Computational Statistics and Data Analysis,
doi: 10.1016/j.csda.2010.11.015.

}
\author{
Fabian Sobotka, Thomas Kneib \cr
University Oldenburg \cr
\url{http://www.uni-oldenburg.de} \cr

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{expectile.laws}}, \code{\link{expectile.boost}}, \code{\link{expectile.restricted}}
}
\examples{
data(dutchboys)

expreg <- expectile.laws(hgt ~ base(age,"pspline"),data=dutchboys,smooth="schall",expectiles=c(0.05,0.2,0.8,0.95))
plot(expreg)

print(expreg)

coef(expreg)

new.d = dutchboys[1:10,]
new.d[,2] = 1:10

predict(expreg,newdata=new.d)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
