\name{expectile.laws}
\Rdversion{1.1}
\alias{expectile.laws}
\alias{expectile.sheets}
\alias{expectile.noncross}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expectile regression of additive models
}
\description{
Generalized additive models are fitted with least asymmetrically weighted squares
or restricted with respect to the mean to obtain the graphs of 11 different expectiles
for continuous, spatial or random effects.
}
\usage{
expectile.laws(formula, data = NULL, smooth = c("schall", "acv", "fixed"), lambda = 0.1, expectiles = NA, parallel = FALSE)

expectile.sheets(formula, data = NULL, smooth = c("acv", "fixed"), lambda = 0.1, lambdap = 5, expectiles = NA, density = FALSE)

expectile.noncross(formula, data  = NULL, smooth = c("schall", "acv", "fixed"), lambda = 0.1, expectiles = NA) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
An R formula object consisting of the response variable, '~'
and the sum of all effects that should be taken into consideration.
Each effect has to be given through the function \code{\link{base}}.
}
  \item{data}{
Optional data frame containing the variables used in the model, if the data is not explicitely given in the
formula.  
}
  \item{smooth}{
There are different smoothing algorithms that should prevent overfitting.
The 'schall' algorithm iterates the smoothing penalty \code{lambda} until it converges,
the asymmetric cross-validation 'acv' minimizes a score-function using \code{\link[stats]{nlm}}
or the function uses a fixed penalty.
}
  \item{lambda}{
The fixed penalty can be adjusted. Also serves as starting value for
the smoothing algorithms.
}
  \item{lambdap}{
A penalty for smoothing across the expectiles in the sheets function. Used if no smoothing algorithm is applied.
}
  \item{expectiles}{
In default setting, the expectiles (0.01,0.02,0.05,0.1,0.2,0.5,0.8,0.9,0.95,0.98,0.99) are calculated.
You may specify your own set of expectiles in a vector. Is overridden by \code{density}.
}
  \item{parallel}{
On Unix machines with the package \code{\link[multicore]{multicore}} installed the different expectiles
can be calculated simultaneously, if the computer has multiple CPU cores. 
}
  \item{density}{
If \code{TRUE}, 99 expectiles from 1\% to 99\% are fitted to allow for a density estimation afterwards.
}
}
\details{
In least asymmetrically weighted squares (LAWS) each expectile is fitted independently from the others.
LAWS minimizes:

\eqn{ S = \sum_{i=1}^{n}{ w_i(p)(y_i - \mu_i(p))^2} }

with

\eqn{ w_i(p) = p 1_{(y_i > \mu_i(p))} + (1-p) 1_{(y_i < \mu_i(p))} }.

The expectile sheets construct a p-spline basis for the expectiles and perform a continuous fit over all expectiles
by fitting the tensor product of the expectile spline basis and the basis of the covariates.
In consequence there will be most likely no crossing of expectiles but also a good fit in heteroscedastic scenarios.
"schall" smoothing does not yet work for sheets.

The noncross function also fits a sheet over all expectiles, but it uses quadratic programming with constraints,
so crossing of expectiles will definitely not happen. So far the function is implemented for one nonlinear or spatial covariate
and further parametric covariates. It works with all smoothing methods.
}
\value{
An object of class 'expectreg', which is basically a list consisting of:
    \item{lambda }{The final smoothing parameters for all expectiles and for all effects in a list.
                   For the restricted and the bundle regression there are only the mean and the residual lambda.}
    \item{intercepts }{The intercept for each expectile.}
    \item{coefficients}{ A matrix of all the coefficients, for each base element
                         a row and for each expectile a column. }
    \item{values}{ The fitted values for each observation and all expectiles,
                   separately in a list for each effect in the model,
                   sorted in order of ascending covariate values. }
    \item{response}{ Vector of the response variable. }
    \item{covariates}{ List with the values of the covariates. }
    \item{formula}{ The formula object that was given to the function. }
    \item{expectiles}{ Vector of fitted expectile asymmetries as given by argument \code{expectiles}. }
    \item{effects}{ List of characters giving the types of covariates. }
    \item{helper}{ List of additional parameters like neighbourhood structure for spatial effects or 'phi' for kriging. }
    \item{design}{ Complete design matrix. }
    \item{fitted}{ Fitted values \eqn{ \hat{y} }. }
\code{\link[=plot.expectreg]{plot}}, \code{\link[=predict.expectreg]{predict}}, \code{\link[=resid.expectreg]{resid}}, \code{\link[=fitted.expectreg]{fitted}} and \code{\link[=effects.expectreg]{effects}}
methods are available for class 'expectreg'.
}
\references{
Schnabel S and Eilers P (2009)
\emph{ Optimal expectile smoothing }
Computational Statistics and Data Analysis, 53:4168-4177

Sobotka F and Kneib T (2010)
\emph{ Geoadditive Expectile Regression }
Computational Statistics and Data Analysis,
doi: 10.1016/j.csda.2010.11.015.

Schnabel S and Eilers P (2011)
\emph{ Expectile sheets for joint estimation of expectile curves }
(under review at Statistical Modelling)
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

Linda Schulze Waltrup \cr
University Bielefeld \cr
\url{http://www.uni-bielefeld.de} \cr
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{base}}, \code{\link{expectile.boost}}, \code{\link{expectile.restricted}}
}
\examples{
data(dutchboys)

expreg <- expectile.laws(dutchboys[,3] ~ base(dutchboys[,2],"pspline"),smooth="schall",expectiles=c(0.05,0.2,0.8,0.95))
plot(expreg)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nonparametric }
\keyword{ smooth }% __ONLY ONE__ keyword per line
\keyword{ multivariate }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ models }
