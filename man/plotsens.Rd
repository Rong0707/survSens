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
plotsens(tau.res, zetaz, zetat, tau1, coeff0)
}
\arguments{
  \item{tau.res}{a data.frame that can be generated from either \code{survSensitivity} or \code{comprSensitivity.}}
  \item{zetaz}{the name of sensitivity parameter in the treatment model.}
  \item{zetat}{the name of sensitivity parameter in the response model.}
  \item{tau1}{the name of estimated treatment effect.}
  \item{coeff0}{the value of estimated treatment effect ignoring any confounding.}
}

\details{This function gives a contour plot in order to visualize results from either \code{survSensitivity} or \code{comprSensitivity}.}
\value{
A contour plot corresponding to the output from either \code{survSensitivity} or \code{comprSensitivity.}
}

\author{
Rong Huang
}
\examples{
data(tau.res) #an example output
plotsens(tau.res, "zetaz", "zetat", "tau1", coeff0 = 1.131)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{contour plot}% use one of  RShowDoc("KEYWORDS")
