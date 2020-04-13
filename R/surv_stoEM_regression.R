surv_stoEM_regression <- function(data, zetat, zetaz, B = 100, theta = 0.5, offset = TRUE){
  t = data$t
  d = data$d
  Z = data$Z
  X = data.matrix(data$X)
  nx = dim(X)[2]

  if(offset){
    #coefficients with simulated U
    Coefz = Coefz.se = matrix(0, nrow = B, ncol = nx+1)
    Coeft1 = Coeft1.se = matrix(0, nrow = B, ncol = nx+1)

    for (j in 1:B){
      Usim = SimulateU_surv(t, d, Z, X, zetat = zetat, zetaz = zetaz, theta = theta, offset = TRUE)

      Z.fit = glm(Z ~ X, offset = zetaz * Usim$U, family=binomial(link="probit"))
      Coefz[j,] = Z.fit$coefficients
      Coefz.se[j,] = summary(Z.fit)$coefficients[,'Std. Error']

      t1.fit = coxph(Surv(t, d) ~ X + Z + offset(zetat * Usim$U))
      Coeft1[j,] = t1.fit$coefficients
      Coeft1.se[j,] = summary(t1.fit)$coefficients[,'se(coef)']
    }

    colnames(Coefz) = colnames(Coefz.se) = names(Z.fit$coefficients)
    colnames(Coeft1) = colnames(Coeft1.se) = names(t1.fit$coefficients)
  }
  else{
    #coefficients with simulated U
    Coefz = Coefz.se = matrix(0, nrow = B, ncol = nx + 2)
    Coeft1 = Coeft1.se = matrix(0, nrow = B, ncol = nx + 2)

    for (j in 1:B){
      Usim = SimulateU_surv(t, d, Z, X, zetat = zetat, zetaz = zetaz, theta = theta, offset = FALSE)

      Z.fit = glm(Z ~ X + Usim$U, family=binomial(link="probit"))
      Coefz[j,] = Z.fit$coefficients
      Coefz.se[j,] = summary(Z.fit)$coefficients[,'Std. Error']

      t1.fit = coxph(Surv(t, d) ~ X + Z + Usim$U)
      Coeft1[j,] = t1.fit$coefficients
      Coeft1.se[j,] = summary(t1.fit)$coefficients[,'se(coef)']
    }
    colnames(Coefz) = colnames(Coefz.se) = names(Z.fit$coefficients)
    colnames(Coeft1) = colnames(Coeft1.se) = names(t1.fit$coefficients)
  }

  tau1 = mean(Coeft1[,"Z"])
  tau1.se = sqrt(mean((Coeft1.se[,"Z"])^2) + (1+1/B) * var(Coeft1[,"Z"]))

  return (list(tau1 = tau1, tau1.se = tau1.se))
}
