compr_stoEM_ipw <- function(data, zetat, zetaz, theta = 0.5, B = 100, offset = TRUE){
  t = data$t
  d1 = (data$d==1)
  d2 = (data$d==2)
  Z = data$Z
  X = data.matrix(data$X)
  nx = dim(X)[2]
  n = length(t)

  #coefficients with simulated U
  Coeft1 = Coeft1.se = numeric(B)
  Coeft2 = Coeft2.se = numeric(B)

  for (j in 1:B){
    Usim = SimulateU_compr(t, cbind(d1,d2), Z, X, zetat = zetat, zetaz = zetaz, theta = theta, offset = TRUE)

    Z.fit = glm(Z ~ X, offset = zetaz * Usim$U, family=binomial(link="probit"))
    #ps = pnorm(cbind(1,X,Usim$U) %*% c(Z.fit$coefficients, zetaz))
    ps = Z.fit$fitted.values
    ipw = (sum(Z)/n) * Z/ps + (1-(sum(Z)/n))*(1-Z)/(1-ps)
    ipw = pmin(ipw, 10)
    ipw = pmax(ipw, 0.1)

    t1.ipw = coxph(Surv(t, d1) ~ Z, weights = ipw, robust = TRUE)
    Coeft1[j] = t1.ipw$coefficients
    Coeft1.se[j] = summary(t1.ipw)$coefficients[,'robust se']

    t2.ipw = coxph(Surv(t, d2) ~ Z, weights = ipw, robust = TRUE)
    Coeft2[j] = t2.ipw$coefficients
    Coeft2.se[j] = summary(t2.ipw)$coefficients[,'robust se']
  }

  names(Coeft1) = names(Coeft1.se) = names(t1.ipw$coefficients)
  names(Coeft2) = names(Coeft2.se) = names(t2.ipw$coefficients)

  tau1 = mean(Coeft1)
  tau1.se = sqrt(mean((Coeft1.se)^2) + (1+1/B) * var(Coeft1))

  tau2 = mean(Coeft2)
  tau2.se = sqrt(mean((Coeft2.se)^2) + (1+1/B) * var(Coeft2))

  return (list(tau1 = tau1, tau1.se = tau1.se, tau2 = tau2, tau2.se = tau2.se))
}
