\name{tau.res}
\alias{tau.res}
\docType{data}
\title{
Sensitivity analysis output example
}
\description{
An example output from \code{survSensitivity}.
}
\usage{data("tau.res")}
\format{
  A data frame with 81 observations on the following 7 variables.
  \describe{
    \item{\code{zetaz}}{a numeric vector, corresponding to the sensitivity parameter in the treatment model.}
    \item{\code{zetat1}}{a numeric vector, corresponding to the sensitivity parameter in the response model.}
    \item{\code{tau1}}{a numeric vector, corresponding to the estimated treatment effect.}
    \item{\code{tau1.se}}{a numeric vector, corresponding to the standard error of the estimated treatment effect.}
    \item{\code{pR2z}}{a numeric vector, corresponding to the Rsquared in the treatment model.}
    \item{\code{pR2t1}}{a numeric vector, corresponding to the Rsquared in the response model.}
    \item{\code{t}}{a numeric vector, corresponding to the \code{t} statistic.}
  }
}
\examples{
data(tau.res)
}
\keyword{datasets}
