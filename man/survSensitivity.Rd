\name{survSensitivity}
\alias{survSensitivity}
\title{
Sensitivity analysis of treatment effect to unmeasured confounding with survival outcomes.
}
\description{
\code{survSensitivity} performs dual-parameter sensitivity analysis of treatment effect to unmeasured confounding in observational studies with survival outcomes
}
\usage{
survSensitivity(t, d, Z, X, method, zetaT = seq(-2,2,by=0.5),
zetaZ = seq(-2,2,by=0.5), theta = 0.5, B = 100, Bem = 1000)
}
\arguments{
  \item{t}{survival outcomes.}
  \item{d}{indicator of occurence of event, d == 0 means censoring.}
  \item{Z}{indicator of treatment.}
  \item{X}{pre-treatment covariates.}
  \item{method}{needs to be one of \code{"stoEM_reg"}, \code{"stoEM_IPW"}, \code{"EM_reg"}.}
}
\details{
This function performs dual-parameter sensitivity analysis of treatment effect to unmeasured confounding by either drawing simulated potential confounders U from the conditional distribution of U given observed response, treatment and covariates or the Expecteation-Maximization method.
}
\value{
\item{tau1}{a matrix with zetaz, zetat, tau1, tau1.se and t statistic.}
\item{Sensplot}{a contour plot corresponding to the output tau1.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Rong Huang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
\keyword{sensitivity analysis}
\keyword{survival outcome}
