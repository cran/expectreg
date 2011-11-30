\name{rb}
\Rdversion{1.1}
\alias{rb}
\alias{mono}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates base for a regression based on covariates
}
\description{
Based on given observations a matrix is created that creates a basis e.g.
of splines or a markov random field that is evaluated for each observation.
Additionally a penalty matrix is generated. Shape constraint p-spline bases can
also be specified.
}
\usage{
rb(x, type = c("pspline", "2dspline", "markov", "radial", "krig", 
                 "random", "ridge", "special", "parametric"),
                 B = NA, P = NA, bnd = NA, center = TRUE, by = NA)
                 
mono(x, constraint = c("increase", "decrease", "convex", "concave", "flatend"), by = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data vector, matrix or data frame. In case of '2dspline', 'radial' or 'krig' \code{type}
number of variables of \code{x} has to be 2, more dimensions are allowed in 'ridge' and 'special' \code{type}.
'markov' and 'random' \code{type} require a vector of a factor.
}
  \item{type}{
Character string defining the type of base that is generated for the given variable(s) \code{x}.
Further description of the possible options is given below in \code{details}.
}
  \item{B}{
For the 'special' \code{type} the base \code{B} and penalization matrix \code{P}
are entered manually. The data frame or matrix needs as many rows as observations in \code{x}
and as many columns as \code{P}.
}
  \item{P}{
Square matrix that has to be provided in 'special' case and with 'markov' \code{type}
if no \code{bnd} is given.
}
  \item{bnd}{
Object of class \code{bnd}, required with 'markov' type if \code{P} is not given.
See \code{\link[BayesX]{read.bnd}}.
}
  \item{center}{
Logical to state whether the basis shall be centered in order to fit additive models
with one central intercept.  
}
  \item{by}{
An optional variable defining varying coefficients, either a factor or numeric variable.
Per default treatment coding is used. Note that the main effect needs to be specified in a separate basis.
}
  \item{constraint}{
Character string defining the type of shape constraint that is imposed on the spline curve. The last option
'flatend' results in constant functions at the covariate edges.
}
}
\details{

Possible \code{types} of bases:
\tabular{ll}{
pspline \tab Penalized splines made upon 20 equidistant knots
             and with degree 2. The penalization matrix consists 
             of differences of the second order, 
             see \code{\link[base]{diff}}. \cr
2dspline \tab Tensor product of 2 p-spline bases 
              with the same properties as above. \cr
markov \tab Gaussian markov random field with a neighbourhood 
            structure given by \code{P} or \code{bnd}. \cr
radial \tab 2-dimensional base, knots are subset of observations, 
            base is calculated as r^2*log(r) with r equalling 
            the euclidean distance between the current observation 
            and the knots. The base has size 50 if possible.\cr
krig \tab 'kriging' produces basically the same base as 'radial', 
          but the base is calculated as exp(-r/phi)*(1+r/phi)
          where \code{phi} is the maximum euclidean distance 
          between two knots divided by a constant. \cr
random \tab A 'random' effect is like the 'markov' random field 
            based on a categorial variable, 
            and since there is no neighbourhood structure, P = I. \cr
ridge \tab In a 'ridge' regression, the base is made 
           from the independent variables while the goal is 
           to determine significant variables from the coefficients. 
           Therefore no penalization is used (P = I). \cr
special \tab In the 'special' case, \code{B} and \code{P} are user defined.\cr
parametric \tab A parametric effect. \cr
}

}
\value{
List consisting of:
  \item{B }{ Matrix of the evaluated base, one row for each observation,
             one column for each base element. }
  \item{P }{ Penalty square matrix, needed for the smoothing in the regression. }
  \item{x }{ The observations \code{x} given to the function. }
  \item{type}{ The \code{type} as given to the function. }
  \item{bnd}{ The \code{bnd} as given to the function, only needed with 'markov' \code{type}. }
  \item{Zspathelp}{ Matrix that is also only needed with 'markov' \code{type}
                    for calculation of the fitted values. }
  \item{phi}{ Constant only needed with 'kriging' \code{type},
              otherwise 'NA'. }
  \item{center}{ The boolean value of the argument \code{center}.
         }
  \item{by}{ The variable included in the \code{by} argument if available. }
  \item{xname}{ Name of the variable \code{x} given to the function. }
}
\references{
    Fahrmeir L and Kneib T and Lang S (2009)
    \emph{ Regression }
    Springer, New York
}
\author{
Fabian Sobotka, Thomas Kneib \cr
University Oldenburg \cr
\url{http://www.uni-oldenburg.de}

Sabine Schnabel \cr
Wageningen University and Research Centre \cr
\url{http://www.wur.nl}

Paul Eilers \cr
Erasmus Medical Center Rotterdam \cr
\url{http://www.erasmusmc.nl}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{quant.bundle}}, \code{\link{expectreg.ls}}
}
\examples{
x <- rnorm(100)

bx <- rb(x,"pspline")

y <- sample(10,100,replace=TRUE)

by <- rb(y,"random")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ models }
\keyword{ spatial }
