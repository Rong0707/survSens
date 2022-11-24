\name{plotsens}
\alias{plotsens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A contour plot of sensitivity analysis results.
}
\description{
A contour plot of sensitivity analysis results.
}
\usage{
plotsens(tau.res, coeff0, partialRsq = FALSE)
}
\arguments{
  \item{tau.res}{a data.frame that can be generated from either \code{survSensitivity} or \code{comprSensitivity.}}
  \item{coeff0}{the value of estimated treatment effect ignoring any confounding.}
  \item{partialRsq}{whether to use partial R^2 instead of coefficients in the contour plot.}
}

\details{This function gives a contour plot in order to visualize results from either \code{survSensitivity} or \code{comprSensitivity}. The name of sensitivity parameter in the treatment model needs to be "zetaz", the name of sensitivity parameter in the response model needs to be "zetat1", and the name of estimated treatment effect needs to be "tau1".}
\value{
A contour plot corresponding to the output from either \code{survSensitivity} or \code{comprSensitivity.}
}

\author{
Rong Huang
}
\examples{
data(tau.res) #an example output
plotsens(tau.res, coeff0 = 1.131)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{contour plot}% use one of  RShowDoc("KEYWORDS")
