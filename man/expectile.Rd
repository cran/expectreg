\name{expectile}
\Rdversion{1.1}
\alias{expectile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expectile regression of additive models
}
\description{
Expectiles are fitted to univariate samples with least asymmetrically weighted squares for asymmetries between 0 and 1.
}
\usage{
expectile(x, probs = seq(0, 1, 0.25), dec = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numeric vector of univariate observations.
}
  \item{probs}{
Numeric vector of asymmetries between 0 and 1 where 0.5 corresponds to the mean.  
}
  \item{dec}{
Number of decimals remaining after rounding the results.  
}
}
\details{
In least asymmetrically weighted squares (LAWS) each expectile is fitted independently from the others.
LAWS minimizes:

\eqn{ S = \sum_{i=1}^{n}{ w_i(p)(x_i - \mu(p))^2} }

with

\eqn{ w_i(p) = p 1_{(x_i > \mu(p))} + (1-p) 1_{(x_i < \mu(p))} }.

\eqn{\mu(p)} is determined by iteration process with recomputed weights \eqn{w_i(p)}.
}
\value{
Numeric vector with the fitted expectiles.
}
\references{

Sobotka F and Kneib T (2010)
\emph{ Geoadditive Expectile Regression }
Computational Statistics and Data Analysis,
doi: 10.1016/j.csda.2010.11.015.

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
\code{\link{expectile.laws}}
}
\examples{
data(dutchboys)

expectile(dutchboys[,3])

x = rnorm(1000)

expectile(x,probs=c(0.01,0.02,0.05,0.1,0.2,0.5,0.8,0.9,0.95,0.98,0.99))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nonparametric }
\keyword{ univariate }% __ONLY ONE__ keyword per line
