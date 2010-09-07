\name{expectile.laws}
\Rdversion{1.1}
\alias{expectile.laws}
\alias{expectile.restricted}
\alias{expectile.bundle}
\alias{expectile.sheets}
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
expectile.laws(formula, data = NULL, smooth = c("schall", "acv", "none"), lambda = 0.1, expectiles = NA, parallel = FALSE)

expectile.restricted(formula, data = NULL, smooth = c("schall", "acv", "none"), lambda = 0.1, expectiles = NA, density = FALSE)

expectile.bundle(formula, data = NULL, smooth = c("schall", "acv", "none"), lambda = 0.1, expectiles = NA, density = FALSE)

expectile.sheets(formula, data = NULL, smooth = c("acv", "none"), lambda = 0.1, lambdap = 5, expectiles = NA, density = FALSE)
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

The restricted version fits the 0.5 expectile at first and then the residuals.
Afterwards the other expectiles are fitted as deviation by a factor of the residuals from the mean expectile.
This algorithm is based on He(1997). The advantage is that expectile crossing cannot occur,
the disadvantage is a suboptimal fit in certain heteroscedastic settings.
Also, since the number of fits is significantly decreased, the restricted version is much faster.

The expectile bundle has a resemblence to the restricted regression. At first, a trend curve is fitted
and then an iteration is performed between fitting the residuals and calculating the deviation factors
for all the expectiles until the results are stable. Therefore this function shares the (dis)advantages
of the restricted.

The expectile sheets construct a p-spline basis for the expectiles and perform a continuous fit over all expectiles
by fitting the tensor product of the expectile spline basis and the basis of the covariates.
In consequence there will be most likely no crossing of expectiles but also a good fit in heteroscedastic scenarios.
"schall" smoothing does not yet work for sheets.
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
The output of the bundle and restricted function consists of two more elements:
    \item{trend.coef}{ Coefficients of the trend function. }
    \item{residual.coef}{ Vector of the coefficients the residual curve was fitted with. }
    \item{asymmetry}{ Vector of the asymmetry factors for all expectiles. }
A \code{\link{plot}} method is available.
}
\references{
Schnabel S and Eilers P (2009)
\emph{ Optimal expectile smoothing }
Computational Statistics and Data Analysis, 53:4168-4177

He X (1997)
\emph{ Quantile Curves without Crossing }
The American Statistician, 51(2):186-192

Schnabel S and Eilers P (2010)
\emph{ Non crossing expectiles and quantiles}
submitted to Journal of Computational and Graphical Statistics

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
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{base}}, \code{\link{expectile.boost}}
}
\examples{
data(dutchboys)

expreg <- expectile.laws(dutchboys[,3] ~ base(dutchboys[,2],"pspline"),smooth="schall",expectiles=c(0.05,0.2,0.8,0.95))

exprest <- expectile.restricted(dutchboys[,3] ~ base(dutchboys[,2],"pspline"),smooth="acv")

expbund <- expectile.bundle(dutchboys[,3] ~ base(dutchboys[,2],"pspline"),smooth="none")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nonparametric }
\keyword{ smooth }% __ONLY ONE__ keyword per line
\keyword{ multivariate }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ models }
