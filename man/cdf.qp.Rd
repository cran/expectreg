\name{cdf.qp}
\alias{cdf.qp}
\alias{cdf.bundle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of the conditional CDF based on expectile curves
}
\description{
Estimating the CDF of the response for a given value of covariate. Additionally quantiles are computed
from the distribution function which allows for the calculation of regression quantiles.
}
\usage{
cdf.qp(expectreg, x = NA, qout = NA, extrap = FALSE)

cdf.bundle(bundle, qout = NA, extrap = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expectreg, bundle}{
An object of class expectreg or subclass bundle respectively. The number of expectiles should be high enough to ensure accurate estimation. 
One approach would be to take as many expectiles as data points.
}
  \item{x}{
The covariate value where the CDF is estimated. By default the first covariate value.
}
  \item{qout}{
Vector of quantiles that will be computed from the CDF.
}
  \item{extrap}{
If TRUE, extreme quantiles will be extrapolated linearly, otherwise the maximum of the CDF is used.
}
}
\details{
Expectile curves can describe very well the spread and location of a scatterplot. With
a set of curves they give good impression about the nature of the data. This information
can be used to estimate the conditional density from the expectile curves.
The results of the bundle model are especially suited in this case
as only one density will be estimated which can then be modulated to
over the independent variable x. The 
density estimation can be formulated as penalized least squares problem that results in a smooth non-negative
density.
The theoretical values of a quantile regression at this covariate value
are also returned for adjustable probabilities \code{qout}.
}
\value{
A list consisting of
  \item{x}{ vector of expectiles where the CDF is computed. }
  \item{cdf }{ vector of values of the CDF at the expectiles \code{x}. }
  \item{quantiles }{ vector of quantile values estimated from the CDF. }
  \item{qout }{ vector of probabilities for the calculated quantiles. }
%% ... return(list( x = e, cdf = F, quantiles = quant, qout = qout))
}
\references{
Schnabel SK and Eilers PHC (2010)
\emph{ A location scale model for non-crossing expectile curves }
(working paper)

Schulze Waltrup L, Sobotka F, Kauermann G and Kneib T (2011)
\emph{ Comparing Expectiles and Quantiles Regarding Efficiency }
Working Paper.
}
\author{
Goeran Kauermann, Linda Schulze Waltrup \cr
Ludwig Maximilian University Munich \cr
\url{http://www.lmu.de} \cr

Fabian Sobotka \cr
University Oldenburg \cr
\url{http://www.uni-oldenburg.de} \cr

Sabine Schnabel \cr
Wageningen University and Research Centre \cr
\url{http://www.wur.nl} \cr

Paul Eilers \cr
Erasmus Medical Center Rotterdam \cr
\url{http://www.erasmusmc.nl} \cr
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{expectreg.ls}}, \code{\link{expectreg.qp}}
}
\examples{
data(dutchboys)
d <- expectreg.ls(dutchboys[,3] ~ rb(dutchboys[,2]),estimate="restricted")
e = cdf.qp(d,15,extrap=TRUE)
e 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ multivariate }
\keyword{ nonlinear }
