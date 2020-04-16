\name{comprSensitivity}
\alias{comprSensitivity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sensitivity analysis of treatment effect to unmeasured confounding with competing risks outcomes.
}
\description{
\code{comprSensitivity} performs a dual-parameter sensitivity analysis of treatment effect to unmeasured confounding in observational studies with competing risks outcomes.
}
\usage{
comprSensitivity(t, d, Z, X, method, zetaT = seq(-2,2,by=0.5),
zetat2 = 0, zetaZ = seq(-2,2,by=0.5), theta = 0.5, B = 100, Bem = 1000)
}
\arguments{
  \item{t}{survival outcomes with competing risks.}
  \item{d}{indicator of occurrence of event, with \code{d == 0} denotes right censoring, \code{d==1} denotes event of interest, \code{d==2} denotes competing risk.}
  \item{Z}{indicator of treatment.}
  \item{X}{pre-treatment covariates, will be included in the model as measured confounders.}
  \item{method}{needs to be one of \code{"stoEM_reg"}, \code{"stoEM_IPW"} and \code{"EM_reg"}.}
  \item{zetaT}{range of coefficient of \eqn{U} in the event of interest response model.}
  \item{zetat2}{value of coefficient of \eqn{U} in the competing risk response model}
  \item{zetaZ}{range of coefficient of \eqn{U} in the treatment model.}
  \item{theta}{marginal probability of \eqn{U=1}.}
  \item{B}{iteration in the stochastic EM algorithm.}
  \item{Bem}{iteration used to estimate variance-covariance matrix in the EM algorithm.}
}
\details{
This function performs a dual-parameter sensitivity analysis of treatment effect to unmeasured confounding by either drawing simulated potential confounders \eqn{U} from the conditional distribution of \eqn{U} given observed response, treatment and covariates or the Expectation-Maximization algorithm. We assume \eqn{U} is following \eqn{Bernoulli(\pi)} (default 0.5). Given \eqn{Z}, \eqn{X} and \eqn{U}, the hazard rate of jth type of failure of is modeled using the Cox proportional hazards (PH) regression:
\deqn{\lambda_j (t|Z, X, U) = \lambda_{j0} (t) exp( \tau_j Z + X' \beta_j + \zeta_j U).}
Given \eqn{X} and \eqn{U}, \eqn{Z} follows a generalized linear model:
\deqn{P(Z=1|X, U) = \Phi(X' \beta_z + \zeta_z U).}
}

\value{
\item{tau1}{a data.frame with zetaz, zetat1, zetat2, tau1, tau1.se and t statistic in the event of interest response model.}
\item{tau2}{a data.frame with zetaz, zetat, zetat2, tau2, tau2.se and t statistic in the competing risks response model.}
\item{Sensplot}{a contour plot corresponding to the output tau1, which will be stored in a seperate pdf file.}
}
\references{
Huang, R., Xu, R., & Dulai, P. S. (2019). Sensitivity Analysis of Treatment Effect to Unmeasured Confounding in Observational Studies with Survival and Competing Risks Outcomes. arXiv preprint arXiv:1908.01444.
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
#load the dataset included in the package
data(comprdata)
#stochastic EM with regression
tau = comprSensitivity(comprdata$t, comprdata$d, comprdata$Z, comprdata$X, "stoEM_reg", B = 10)
#EM with regression
tau = comprSensitivity(comprdata$t, comprdata$d, comprdata$Z, comprdata$X, "EM_reg", Bem = 100)
}
\keyword{sensitivity analysis}
\keyword{competing risks outcomes}
