\name{survSensitivity}
\alias{survSensitivity}
\title{
Sensitivity analysis of treatment effect to unmeasured confounding with survival outcomes.
}
\description{
\code{survSensitivity} performs a dual-parameter sensitivity analysis of treatment effect to unmeasured confounding in observational studies with survival outcomes.
}
\usage{
survSensitivity(t, d, Z, X, method, zetaT = seq(-2,2,by=0.5),
zetaZ = seq(-2,2,by=0.5), theta = 0.5, B = 50, Bem = 200, contour = TRUE)
}
\arguments{
  \item{t}{survival outcomes.}
  \item{d}{indicator of occurrence of event, with \code{d == 0} denotes right censoring.}
  \item{Z}{indicator of treatment.}
  \item{X}{pre-treatment covariates that will be included in the model as measured confounders.}
  \item{method}{needs to be one of \code{"stoEM_reg"}, \code{"stoEM_IPW"}, and \code{"EM_reg"}.}
  \item{zetaT}{range of coefficient of \eqn{U} in the response model.}
  \item{zetaZ}{range of coefficient of \eqn{U} in the treatment model.}
  \item{theta}{marginal probability of \eqn{U=1}.}
  \item{B}{iteration in the stochastic EM algorithm.}
  \item{Bem}{iteration used to estimate the variance-covariance matrix in the EM algorithm.}
  \item{contour}{a contour plot will be generated, default \code{TRUE}}
}
\details{
This function performs a dual-parameter sensitivity analysis of treatment effect to unmeasured confounding by either drawing simulated potential confounders \eqn{U} from the conditional distribution of \eqn{U} given observed response, treatment and covariates or the Expectation-Maximization algorithm. We assume \eqn{U} is following \eqn{Bernoulli(\pi)} (default 0.5). Given \eqn{Z}, \eqn{X} and \eqn{U}, the hazard rate is modeled using the Cox proportional hazards (PH) regression:
\deqn{\lambda (t | Z, X, U) = \lambda_{0} (t) exp(\tau Z + X ' \beta + \zeta U).}
Given \eqn{X} and \eqn{U}, \eqn{Z} follows a generalized linear model:
\deqn{P( Z=1 | X,U ) = \Phi(X' \beta_z + \zeta_z U).}
}
\value{
\item{tau}{a data.frame with zetaz, zetat, tau1, tau1.se and t statistic.}
\item{Sensplot}{if \code{contour = TRUE}, a contour plot corresponding to the output tau, which will be stored in a separate pdf file.}
}
\references{
Huang, R., Xu, R., & Dulai, P. S. (2019). Sensitivity Analysis of Treatment Effect to Unmeasured Confounding in Observational Studies with Survival and Competing Risks Outcomes. arXiv preprint arXiv:1908.01444.
}
\author{
Rong Huang
}
\examples{
#load the dataset included in the package.
data(survdata)
#stochastic EM with regression
tau.sto = survSensitivity(survdata$t, survdata$d, survdata$Z, survdata$X,
"stoEM_reg", zetaT = 0.5, zetaZ = 0.5, B = 3, contour = FALSE)

#EM with regression
tau.em = survSensitivity(survdata$t, survdata$d, survdata$Z, survdata$X,
"EM_reg", zetaT = 0.5, zetaZ = 0.5, Bem = 50, contour = FALSE)
}
\keyword{sensitivity analysis}
\keyword{survival outcomes}
