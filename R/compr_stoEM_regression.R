compr_stoEM_regression <- function(data, zetat, zetaz, theta = 0.5, B = 100, offset = TRUE){
  t = data$t
  d1 = (data$d==1)
  d2 = (data$d==2)
  Z = data$Z
  X = data.matrix(data$X)
  n = length(t)
  nx = dim(X)[2]

  if(offset){
    # Record coefficients with simulated U
    Coefz = Coefz.se = matrix(0, nrow = B, ncol = nx+1)
    Coeft1 = Coeft1.se = matrix(0, nrow = B, ncol = nx+1)
    Coeft2 = Coeft2.se = matrix(0, nrow = B, ncol = nx+1)
    partialR2z = c()
    partialR2t1 = c()
    partialR2t2 = c()

    for (j in 1:B){
      Usim = SimulateU_compr(t, cbind(d1,d2), Z, X, zetat = zetat, zetaz = zetaz, theta = theta, offset = TRUE)

      Z.fit = glm(Z ~ X, offset = zetaz * Usim$U, family=binomial(link="probit"))
      Coefz[j,] = Z.fit$coefficients
      Coefz.se[j,] = summary(Z.fit)$coefficients[,'Std. Error']

      # Calculate partial R-sq of z ~ u | x
      Z.fit_reduced = glm(Z ~ X, family=binomial(link="probit"))
      if (zetaz >= 0)
        partialR2z[j] = 1 - exp((Z.fit$deviance-Z.fit_reduced$deviance)/n)
      else
        partialR2z[j] = - (1 - exp((Z.fit$deviance-Z.fit_reduced$deviance)/n))

      t1.fit = coxph(Surv(t, d1) ~ X + Z + offset(zetat[1] * Usim$U))
      Coeft1[j,] = t1.fit$coefficients
      Coeft1.se[j,] = summary(t1.fit)$coefficients[,'se(coef)']

      # Calculate partial R-sq of (t, d1) ~ u | x, z
      t1.fit_reduced = coxph(Surv(t, d1) ~ X + Z)
      logtest <- -2 * (t1.fit_reduced$loglik[2] - t1.fit$loglik[2])
      if (zetat[1] >= 0)
        partialR2t1[j] = (1 - exp(-logtest/t1.fit$nevent))
      else
        partialR2t1[j] = - (1 - exp(-logtest/t1.fit$nevent))

      t2.fit = coxph(Surv(t, d2) ~ X + Z + offset(zetat[2] * Usim$U))
      Coeft2[j,] = t2.fit$coefficients
      Coeft2.se[j,] = summary(t2.fit)$coefficients[,'se(coef)']

      # Calculate partial R-sq of (t, d2) ~ u | x, z
      t2.fit_reduced = coxph(Surv(t, d2) ~ X + Z)
      logtest <- -2 * (t2.fit_reduced$loglik[2] - t2.fit$loglik[2])
      if (zetat[2] >= 0)
        partialR2t2[j] = (1 - exp(-logtest/t2.fit$nevent))
      else
        partialR2t2[j] = - (1 - exp(-logtest/t2.fit$nevent))
    }
    colnames(Coefz) = colnames(Coefz.se) = names(Z.fit$coefficients)
    colnames(Coeft1) =  colnames(Coeft1.se) =  names(t1.fit$coefficients)
    colnames(Coeft2) = colnames(Coeft2.se) = names(t2.fit$coefficients)
  }
  else{
    # Record coefficients with simulated U
    Coefz = Coefz.se =  matrix(0, nrow = B, ncol = nx+2)
    Coeft1 = Coeft1.se = matrix(0, nrow = B, ncol = nx+2)
    Coeft2 = Coeft2.se = matrix(0, nrow = B, ncol = nx+2)
    partialR2z = c()
    partialR2t1 = c()
    partialR2t2 = c()

    for (j in 1:B){
      Usim = SimulateU_compr(t, cbind(d1,d2), Z, X, zetat = zetat, zetaz = zetaz, theta = theta, offset = FALSE)

      Z.fit = glm(Z ~ X + Usim$U, family=binomial(link="probit"))
      Coefz[j,] = Z.fit$coefficients
      Coefz.se[j,] = summary(Z.fit)$coefficients[,'Std. Error']

      # Calculate partial R-sq of z ~ u | x
      Z.fit_reduced = glm(Z ~ X, family=binomial(link="probit"))
      if (zetaz >= 0)
        partialR2z[j] = 1 - exp((Z.fit$deviance-Z.fit_reduced$deviance)/n)
      else
        partialR2z[j] = - (1 - exp((Z.fit$deviance-Z.fit_reduced$deviance)/n))

      t1.fit = coxph(Surv(t, d1) ~ X + Z + Usim$U)
      Coeft1[j,] = t1.fit$coefficients
      Coeft1.se[j,] = summary(t1.fit)$coefficients[,'se(coef)']

      # Calculate partial R-sq of (t, d1) ~ u | x, z
      t1.fit_reduced = coxph(Surv(t, d1) ~ X + Z)
      logtest <- -2 * (t1.fit_reduced$loglik[2] - t1.fit$loglik[2])
      if (zetat[1] >= 0)
        partialR2t1[j] = (1 - exp(-logtest/t1.fit$nevent))
      else
        partialR2t1[j] = - (1 - exp(-logtest/t1.fit$nevent))

      t2.fit = coxph(Surv(t, d2) ~ X + Z + Usim$U)
      Coeft2[j,] = t2.fit$coefficients
      Coeft2.se[j,] = summary(t2.fit)$coefficients[,'se(coef)']

      # Calculate partial R-sq of (t, d2) ~ u | x, z
      t2.fit_reduced = coxph(Surv(t, d2) ~ X + Z)
      logtest <- -2 * (t2.fit_reduced$loglik[2] - t2.fit$loglik[2])
      if (zetat[2] >= 0)
        partialR2t2[j] = (1 - exp(-logtest/t2.fit$nevent))
      else
        partialR2t2[j] = - (1 - exp(-logtest/t2.fit$nevent))
    }
    colnames(Coefz) = colnames(Coefz.se) = names(Z.fit$coefficients)
    colnames(Coeft1) =  colnames(Coeft1.se) =  names(t1.fit$coefficients)
    colnames(Coeft2) = colnames(Coeft2.se) = names(t2.fit$coefficients)
  }

  tau1 = mean(Coeft1[,"Z"])
  tau1.se = sqrt(mean((Coeft1.se[,"Z"])^2) + (1+1/B) * var(Coeft1[,"Z"]))

  tau2 = mean(Coeft2[,"Z"])
  tau2.se = sqrt(mean((Coeft2.se[,"Z"])^2) + (1+1/B) * var(Coeft2[,"Z"]))

  pR2z = mean(partialR2z)
  pR2t1 = mean(partialR2t1)
  pR2t2 = mean(partialR2t2)

  return (list(tau1 = tau1, tau1.se = tau1.se, tau2 = tau2, tau2.se = tau2.se,
               pR2z = pR2z, pR2t1 = pR2t1, pR2t2 = pR2t2))
}
