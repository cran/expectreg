\name{expectile.cdf}
\alias{expectile.cdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of the conditional CDF based on expectile curves
}
\description{
Estimating the CDF of the response for a given value of covariate. Additionally quantiles are computed
from the distribution function which allows for the calculation of regression quantiles.
}
\usage{
expectile.cdf(expectreg, x = NA, qout = NA, extrap = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expectreg}{
An object of class expectreg. The number of expectiles should be high enough to ensure accurate estimation. 
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
For a given fixed covariate value, the conditional CDF of the response is estimated from the
calculated expectiles. The theoretical values of a quantile regression at this covariate value
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
Schulze Waltrup L, Sobotka F, Kauermann G and Kneib T (2011)
\emph{ Comparing Expectiles and Quantiles Regarding Efficiency }
Working Paper.
}
\author{
Goeran Kauermann, Linda Schulze Waltrup \cr
University Bielefeld \cr
\url{http://www.uni-bielefeld.de} \cr

Fabian Sobotka \cr
University Oldenburg \cr
\url{http://www.uni-oldenburg.de} \cr
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bundle.density}}
}
\examples{
data(dutchboys)
d <- expectile.restricted(dutchboys[,3] ~ base(dutchboys[,2]))
e = expectile.cdf(d,15,extrap=TRUE)
e 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ multivariate }
\keyword{ nonlinear }
