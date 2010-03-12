\name{expectile.boost}
\Rdversion{1.1}
\alias{expectile.boost}
\alias{quant.boost}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantile and expectile regression using boosting
}
\description{
Generalized additive models are fitted with gradient boosting for optimizing arbitrary 
loss functions to obtain the graphs of 11 different expectiles
for continuous, spatial or random effects.
}
\usage{
expectile.boost(formula, data, mstop = NA, parallel = FALSE)

quant.boost(formula, data, mstop = NA, parallel = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
An R formula object consisting of the response variable, '~'
and the sum of all effects that should be taken into consideration (see \code{\link[mboost]{gamboost}}).
Each effect can be linear or represented through a nonlinear or spatial base
(see \code{\link[mboost]{bbs}}). Each variable has to be named consistently with \code{data}.
}
  \item{data}{
data frame or matrix.
}
  \item{mstop}{
vector, number of bootstrap iterations for each of the 11 quantiles/expectiles that are fitted.
Default is 10 000.
}
  \item{parallel}{
On Unix machines with the package \code{\link[multicore]{multicore}} installed the different expectiles
can be calculated simultaneously, if the computer has multiple CPUs. 
}
}
\details{
A (generalized) additive model is fitted using a boosting algorithm based on component-wise univariate base learners. 
The base learner can be specified via the formula object. After fitting the model a cross-validation is done using
\code{\link[mboost]{cvrisk}} to determine the optimal stopping point for the boosting which results in the best fit.
}
\value{
An object of class 'expectreg', which is basically a list consisting of:
    \item{values}{ The fitted values for each observation and all expectiles,
                   separately in a list for each effect in the model,
                   sorted in order of ascending covariate values. }
    \item{response}{ Vector of the response variable. }
    \item{formula}{ The formula object that was given to the function. }
}
\references{
Fenske N and Kneib T and Hothorn T (2009)
\emph{ Identifying Risk Factors for Severe Childhood Malnutrition
by Boosting Additive Quantile Regression}
Technical Report 052, University of Munich
}
\author{
Fabian Sobotka, Thomas Kneib \cr
University Oldenburg \cr
\url{http://www.uni-oldenburg.de} \cr
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{expectile.laws}}, \code{\link[mboost]{gamboost}}, \code{\link[mboost]{bbs}}, \code{\link[mboost]{cvrisk}}
}
\examples{
data(dutchboys)

expreg <- expectile.boost(hgt ~ bbs(age,df=5,degree=2),dutchboys,mstop=rep(600,11))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nonparametric }
\keyword{ smooth }% __ONLY ONE__ keyword per line
\keyword{ multivariate }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ models }
