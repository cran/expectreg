\name{enorm}
\Rdversion{1.1}
\alias{enorm}
\alias{ebeta}
\alias{eunif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expectiles of distributions
}
\description{
Much like the 0.5 quantile of a distribution is the median,
the 0.5 expectile is the mean / expected value.
These functions add the possibility of calculating expectiles
of known distributions.
}
\usage{
enorm(e, m = 0, sd = 1)

ebeta(e, a = 1, b = 1)

eunif(e, min = 0, max = 1, lower.tail = TRUE, log.p = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{e}{
vector of probabilities.
}
  \item{m, sd}{
mean and standard deviation of the Normal distribution.
}
  \item{a, b}{
positive parameters of the Beta distribution.
}
  \item{min, max, lower.tail, log.p}{
minimum, maximum of the uniform distribution, and if
probabilities are given from the lower tail and in log.
}
}
\details{
An expectile of a distribution cannot be determined explicitely,
but instead is given by an equation.
The expectile z for a probability p is:
\eqn{ p = \frac{G(z) - z F(z)}{2(G(z) - z F(z)) + z - m}}
where m is the mean, F the cdf and G the partial moment function
\eqn{ G(z) = \int\limits_{-\infty}^{z} uf(u) \mbox{d}u }.
}
\value{
Vector of the expectiles for the desired distribution.
}
\references{
Newey W and Powell J (1987)
\emph{ Asymmetric least squares estimation and testing }
Econometrica, 55:819-847
}
\author{
Fabian Sobotka, Thomas Kneib \cr
University Oldenburg \cr
\url{http://www.uni-oldenburg.de} \cr
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ekoenker}}
}
\examples{
x <- seq(0,1,0.2)

qnorm(x)
enorm(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
