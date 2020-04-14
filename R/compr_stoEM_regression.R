compr_stoEM_regression <- function(data, zetat, zetaz, theta = 0.5, B = 100, offset = TRUE){
  t = data$t
  d1 = (data$d==1)
  d2 = (data$d==2)
  Z = data$Z
  X = data.matrix(data$X)
  n = length(t)
  nx = dim(X)[2]

  if(offset){
    #coefficients with simulated U
    Coefz = Coefz.se = matrix(0, nrow = B, ncol = nx+1)
    Coeft1 = Coeft1.se = matrix(0, nrow = B, ncol = nx+1)
    Coeft2 = Coeft2.se = matrix(0, nrow = B, ncol = nx+1)

    for (j in 1:B){
      Usim = SimulateU(t, cbind(d1,d2), Z, X, zetat = zetat, zetaz = zetaz, theta = theta, offset = TRUE)

      Z.fit = glm(Z ~ X, offset = zetaz * Usim$U, family=binomial(link="probit"))
      Coefz[j,] = Z.fit$coefficients
      Coefz.se[j,] = summary(Z.fit)$coefficients[,'Std. Error']

      t1.fit = coxph(Surv(t, d1) ~ X + Z + offset(zetat[1] * Usim$U))
      Coeft1[j,] = t1.fit$coefficients
      Coeft1.se[j,] = summary(t1.fit)$coefficients[,'se(coef)']

      t2.fit = coxph(Surv(t, d2) ~ X + Z + offset(zetat[2] * Usim$U))
      Coeft2[j,] = t2.fit$coefficients
      Coeft2.se[j,] = summary(t2.fit)$coefficients[,'se(coef)']
    }
    colnames(Coefz) = colnames(Coefz.se) = names(Z.fit$coefficients)
    colnames(Coeft1) =  colnames(Coeft1.se) =  names(t1.fit$coefficients)
    colnames(Coeft2) = colnames(Coeft2.se) = names(t2.fit$coefficients)
  }
  else{
    #coefficients with simulated U
    Coefz = matrix(0, nrow = B, ncol = nx+2)
    Coeft1 = matrix(0, nrow = B, ncol = nx+2)
    Coeft2 = matrix(0, nrow = B, ncol = nx+2)

    for (j in 1:B){
      Usim = SimulateU(t, cbind(d1,d2), Z, X, zetat = zetat, zetaz = zetaz, theta = theta, offset = FALSE)

      Z.fit = glm(Z ~ X + Usim$U, family=binomial(link="probit"))
      Coefz[j,] = Z.fit$coefficients
      Coefz.se[j,] = summary(Z.fit)$coefficients[,'Std. Error']

      t1.fit = coxph(Surv(t, d1) ~ X + Z + Usim$U)
      Coeft1[j,] = t1.fit$coefficients
      Coeft1.se[j,] = summary(t1.fit)$coefficients[,'se(coef)']

      t2.fit = coxph(Surv(t, d2) ~ X + Z + Usim$U)
      Coeft2[j,] = t2.fit$coefficients
      Coeft2.se[j,] = summary(t2.fit)$coefficients[,'se(coef)']
    }
    colnames(Coefz) = colnames(Coefz.se) = names(Z.fit$coefficients)
    colnames(Coeft1) =  colnames(Coeft1.se) =  names(t1.fit$coefficients)
    colnames(Coeft2) = colnames(Coeft2.se) = names(t2.fit$coefficients)
  }

  tau1 = mean(Coeft1[,"Z"])
  tau1.se = sqrt(mean((Coeft1.se[,"Z"])^2) + (1+1/B) * var(Coeft1[,"Z"]))

  tau2 = mean(Coeft2[,"Z"])
  tau2.se = sqrt(mean((Coeft2.se[,"Z"])^2) + (1+1/B) * var(Coeft2[,"Z"]))

  return (list(tau1 = tau1, tau1.se = tau1.se, tau2 = tau2, tau2.se = tau2.se))
}
