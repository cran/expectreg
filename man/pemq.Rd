\name{pkoenker}
\Rdversion{1.1}
\alias{pemq}
\alias{demq}
\alias{remq}
\alias{qemq}
\alias{eemq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The "expectiles-meet-quantiles" distribution family.
}
\description{
Density, distribution function, quantile function, random generation and expectile function
for a family of distributions for which expectiles and quantiles coincide.
}
\usage{
pemq(z,ncp=0,s=1)
demq(z,ncp=0,s=1)
qemq(q,ncp=0,s=1)
remq(n,ncp=0,s=1)
eemq(q,ncp=0,s=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ncp}{
non centrality parameter and mean of the distribution.
}
  \item{s}{
scaling parameter, has to be positive.
}
  \item{z}{
vector of quantiles.
}
  \item{q}{
vector of asymmetries.
}
  \item{n}{
number of observations. If length(n) > 1, the length is taken to be the number required.
}
}
\details{
This distribution has the cumulative distribution function:
\eqn{ F(x;ncp,s) = \frac{1}{2}(1 + sgn(\frac{x-ncp}{s}) \sqrt{1 - \frac{2}{2 + (\frac{x-ncp}{s})^2}}) }

and the density:
\eqn{ f(x;ncp,s) = \frac{1}{s}( \frac{1}{2 + (\frac{x-ncp}{s})^2} )^\frac{3}{2}
}

It has infinite variance, still can be scaled by the parameter \code{s}.
It has mean \code{ncp}.
In the canonical parameters it is equal to a students-t distribution with 2 degrees of freedom.
For \eqn{ s = \sqrt{2} } it is equal to a distribution introduced by Koenker(2005).
}

\value{
\code{demq} gives the density, \code{pemq} gives the distribution function, 
\code{qemq} gives the quantile function, \code{eemq} computes the expectiles numerically and is only provided for completeness,
since the quantiles = expectiles can be determined analytically using \code{qemq}, 
and \code{remq} generates random deviates.
}
\references{
Koenker R (2005)
\emph{ Quantile Regression }
Cambridge University Press, New York
}
\author{
Fabian Sobotka, Thomas Kneib \cr
University Oldenburg \cr
\url{http://www.uni-oldenburg.de} \cr
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{enorm}}
}
\examples{
x <- seq(-5,5,length=100)
plot(x,demq(x))
plot(x,pemq(x,ncp=1))

z <- remq(100,s=sqrt(2))

y <- seq(0,1,0.2)
qemq(y)
eemq(y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
