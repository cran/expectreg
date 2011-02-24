\name{bundle.density}
\Rdversion{1.1}
\alias{bundle.density}
\title{
Calculation of the conditional density based on expectile curves}
\description{
Based on the results of an expectile bundle model a smooth conditional denisty underlying the
data is estimated.
}
\usage{bundle.density(bundle)}

\arguments{
    \item{bundle}{A bundle object is the output of the function \code{\link{expectile.bundle}}.}
  }

\details{Expectile curves can describe very well the spread and location of a scatterplot. With
a set of curves they give good impression about the nature of the data. This information
can be used to estimate the conditional density from the expectile curves.
The results of the bundle model are especially suited in this case
as only one density will be estimated which can then be modulated to
over the independent variable x. The 
density estimation can be formulated as penalized least squares problem that results in a smooth non-negative
density.
}

\value{
\describe{
\item{\code{random}}{standardized residual for the values of the initial data set}
\item{\code{density}}{discrete density values at x}
\item{\code{x}}{regular grid over the range of the standardized residuals}
}
A histogram of \code{random} as well as the estimated density is also produced.
}

\references{
Schnabel SK and Eilers PHC (2010)
\emph{ A location scale model for non-crossing expectile curves }
Under review at Statistics and Computing
}
\examples{
data(dutchboys)
expbund <- expectile.bundle(dutchboys[,3] ~ base(dutchboys[,2],"pspline"),smooth="fixed")
bundle.density(expbund)
}
