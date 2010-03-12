\name{pkoenker}
\Rdversion{1.1}
\alias{pkoenker}
\alias{dkoenker}
\alias{rkoenker}
\alias{qkoenker}
\alias{ekoenker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A special distribution developed by Roger Koenker
}
\description{
Density, distribution function, quantile function, random generation and expectile function
for a special distribution for which expectiles and quantiles coincide.
}
\usage{
pkoenker(z)
dkoenker(z)
qkoenker(q)
rkoenker(n)
ekoenker(q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
vector of quantiles.
}
  \item{q}{
vector of probabilities.
}
  \item{n}{
number of observations. If length(n) > 1, the length is taken to be the number required.
}
}
\details{
This distribution has the distribution function:
\eqn{ F(x) = \frac{1}{2}(1 + sgn(x) \sqrt{1 - \frac{4}{4 + x^2}}) }

and the density:
\eqn{ f(x) = \frac{2|x|}{(4+x^2)^2 \sqrt{1 - \frac{4}{4 + x^2}}}}
}
\value{
dkoenker gives the density, pkoenker gives the distribution function, 
qkoenker gives the quantile function, ekoenker gives the expectile function, 
and rkoenker generates random deviates.
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
x <- seq(-10,10,length=100)
plot(x,dkoenker(x))
plot(x,pkoenker(x))

z <- rkoenker(100)

y <- seq(0,1,0.2)
qkoenker(y)
ekoenker(y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
