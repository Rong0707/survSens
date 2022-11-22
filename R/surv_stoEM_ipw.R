surv_stoEM_ipw <- function(data, zetat, zetaz, B = 100, theta = 0.5){
  t = data$t
  d = data$d
  Z = data$Z
  X = data.matrix(data$X)
  nx = dim(X)[2]
  n = length(t)

  #Record coefficients with simulated U
  Coeft1 = Coeft1.se = numeric(B)
  partialR2z = c()
  partialR2t1 = c()

  for (j in 1:B){
      Usim = SimulateU_surv(t, d, Z, X, zetat = zetat, zetaz = zetaz, theta = theta, offset=TRUE)

      Z.fit = glm(Z ~ X, offset = zetaz * Usim$U, family=binomial(link="probit"))
      # Calculate partial R-sq of z ~ u | x
      Z.fit_reduced = glm(Z ~ X, family=binomial(link="probit"))
      if (zetaz >= 0)
        partialR2z[j] = 1 - Z.fit$deviance/Z.fit_reduced$deviance
      else
        partialR2z[j] = - (1 - Z.fit$deviance/Z.fit_reduced$deviance)

      #ps = pnorm(cbind(1,X,Usim$U) %*% c(Z.fit$coefficients, zetaz))
      ps = Z.fit$fitted.values
      ipw = (sum(Z)/n) * Z/ps + (1-(sum(Z)/n))*(1-Z)/(1-ps)
      ipw = pmin(ipw, 10)
      ipw = pmax(ipw, 0.1)

      t1.ipw = coxph(Surv(t, d) ~ Z, weights = ipw, robust = TRUE)
      Coeft1[j] = t1.ipw$coefficients
      Coeft1.se[j] = summary(t1.ipw)$coefficients[,'robust se']

      # Calculate partial R-sq of (t, d) ~ u | x, z
      t1.fit = coxph(Surv(t, d) ~ X + Z + offset(zetat * Usim$U))
      t1.fit_reduced = coxph(Surv(t, d) ~ X + Z)
      logtest <- -2 * (t1.fit_reduced$loglik[2] - t1.fit$loglik[2])
      if (zetat >= 0)
        partialR2t1[j] = (1 - exp(-logtest/t1.fit$nevent))
      else
        partialR2t1[j] = - (1 - exp(-logtest/t1.fit$nevent))
  }

  names(Coeft1) = names(Coeft1.se) = names(t1.ipw$coefficients)

  tau1 = mean(Coeft1)
  tau1.se = sqrt(mean((Coeft1.se)^2) + (1+1/B) * var(Coeft1))
  pR2z = mean(partialR2z)
  pR2t1 = mean(partialR2t1)

  return (list(tau1 = tau1, tau1.se = tau1.se, pR2z = pR2z, pR2t1 = pR2t1))
}
