\name{expectreg-package}
\Rdversion{1.1}
\alias{expectreg-package}
\alias{expectreg}
\docType{package}
\title{
Expectile Regression
}
\description{
Expectile and quantile regression of models with nonlinear effects
  e.g. spatial, random, ridge using least asymmetric weighed squares / absolutes
  as well as boosting; also supplies expectiles for common distributions.
}
\details{
\tabular{ll}{
Package: \tab expectreg\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2010-03-12\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

\itemize{
\item This package requires the packages \code{\link[BayesX]{BayesX}}, \code{\link[mboost]{mboost}} and \code{\link[splines]{splines}}.
}

}
\author{
Fabian Sobotka, Thomas Kneib \cr
University Oldenburg \cr
\url{http://www.uni-oldenburg.de} \cr

Sabine Schnabel \cr
Wageningen University and Research Centre \cr
\url{http://www.wur.nl}

Paul Eilers \cr
Erasmus Medical Center Rotterdam \cr
\url{http://www.erasmusmc.nl}

Maintainer: Fabian Sobotka <fabian.sobotka@uni-oldenburg.de>

}
\references{
Koenker R (2005)
\emph{ Quantile Regression }
Cambridge University Press, New York

Schnabel S and Eilers P (2009)
\emph{ Optimal expectile smoothing }
Computational Statistics and Data Analysis, 53:4168-4177

Fenske N and Kneib T and Hothorn T (2009)
\emph{ Identifying Risk Factors for Severe Childhood Malnutrition
by Boosting Additive Quantile Regression}
Technical Report 052, University of Munich

He X (1997)
\emph{ Quantile Curves without Crossing }
The American Statistician, 51(2):186-192
}

\keyword{ package }
\keyword{ nonparametric }
\keyword{ smooth }% __ONLY ONE__ keyword per line
\keyword{ multivariate }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ models }
\seealso{
\code{\link[mboost]{mboost}}, \code{\link[BayesX]{BayesX}}
}
\examples{
data(dutchboys)
## Expectile Regression using the restricted approach
expreg <- expectile.restricted(dutchboys[,3] ~ base(dutchboys[,2],"pspline"),smooth="schall")
## The calculation of expectiles for given distributions
enorm(0.1)
}
