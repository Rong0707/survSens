comprSensitivity <- function(t, d, Z, X, method, zetaT = seq(-2,2,by=0.5), zetat2 = 0, zetaZ = seq(-2,2,by=0.5), theta = 0.5, B = 50, Bem = 200){
  n = length(t)
  data = list(t = t, d = d, Z = Z, X = X)
  data1 = data.frame(data)

  tau1.res = data.frame()
  tau2.res = data.frame()

  Z.fit = glm(Z ~ X, family = binomial(link = "probit"))
  Coefz.noU = Z.fit$coefficients

  if(method == "stoEM_reg" | method == "EM_reg"){
    t1.fit = coxph(Surv(t, d==1) ~ X + Z)
    Coeft1.noU = t1.fit$coefficients
    t2.fit = coxph(Surv(t, d==2) ~ X + Z)
    Coeft2.noU = t2.fit$coefficients
  }
  else if(method == "stoEM_IPW"){
    ps = Z.fit$fitted.values
    ipw = (sum(Z)/n) * Z/ps + (1-(sum(Z)/n))*(1-Z)/(1-ps)
    ipw = pmin(ipw, 10)
    ipw = pmax(ipw, 0.1)
    t1.ipw = coxph(Surv(t, d==1) ~ Z, weights = ipw, robust = TRUE)
    Coeft1.noU = t1.ipw$coefficients
    t2.ipw = coxph(Surv(t, d==2) ~ Z, weights = ipw, robust = TRUE)
    Coeft2.noU = t2.ipw$coefficients
  }
  else{
    print("No method found.")
    return(NA)
  }

  for(i in 1:length(zetaZ)){
    for(j in 1:length(zetaT)){
      if(method == "stoEM_reg"){
        #stochastic EM/regression
        tau.sto = compr_stoEM_regression(data, zetat = c(zetaT[j],zetat2), zetaz = zetaZ[i], theta = theta, B = B)
        tau1.res = rbind(tau1.res,
                         data.frame(zetaz = zetaZ[i], zetat1 = zetaT[j], zetat2 = zetat2,
                                    tau1 = tau.sto$tau1, tau1.se = tau.sto$tau1.se,
                                    pR2z = tau.sto$pR2z, pR2t1 = tau.sto$pR2t1, pR2t2 = tau.sto$pR2t2))
        tau2.res = rbind(tau2.res,
                         data.frame(zetaz = zetaZ[i], zetat1 = zetaT[j], zetat2 = zetat2,
                                    tau2 = tau.sto$tau2, tau2.se = tau.sto$tau2.se,
                                    pR2z = tau.sto$pR2z, pR2t1 = tau.sto$pR2t1, pR2t2 = tau.sto$pR2t2))
      }
      else if(method == "stoEM_IPW"){
        #stochastic EM/IPW
        tau.ipw = compr_stoEM_ipw(data, zetat = c(zetaT[j],zetat2), zetaz = zetaZ[i], theta = theta, B = B)
        tau1.res = rbind(tau1.res,
                         data.frame(zetaz = zetaZ[i], zetat1 = zetaT[j], zetat2 = zetat2,
                                    tau1 = tau.ipw$tau1, tau1.se = tau.ipw$tau1.se,
                                    pR2z = tau.ipw$pR2z, pR2t1 = tau.ipw$pR2t1, pR2t2 = tau.ipw$pR2t2))
        tau2.res = rbind(tau2.res,
                         data.frame(zetaz = zetaZ[i], zetat1 = zetaT[j], zetat2 = zetat2,
                                    tau2 = tau.ipw$tau2, tau2.se = tau.ipw$tau2.se,
                                    pR2z = tau.ipw$pR2z, pR2t1 = tau.ipw$pR2t1, pR2t2 = tau.ipw$pR2t2))
      }
      else if(method == "EM_reg"){
        #EM/regression
        tau.em = emU_compr(t, cbind(d==1, d==2), Z, X, zetat = c(zetaT[j],zetat2), zetaz = zetaZ[i], theta = theta)
        data1$p = tau.em$p

        # Calculate partial R-sq of z ~ u | x
        U = rbinom(n,1,tau.em$p)
        partialR2z = c()
        for (k in 1:B){
          Z.fit = glm(Z ~ X, offset = zetaZ[i] * U, family=binomial(link="probit"))
          Z.fit_reduced = glm(Z ~ X, family=binomial(link="probit"))
          if (zetaZ[i] >= 0)
            partialR2z[k] = 1 - Z.fit$deviance/Z.fit_reduced$deviance
          else
            partialR2z[k] = - (1 - Z.fit$deviance/Z.fit_reduced$deviance)
        }

        nx = length(tau.em$z.coef) - 2
        tau.em.final = compr_EM_variance(data1, zetat = c(zetaT[j],zetat2), zetaz = zetaZ[i], z.coef = tau.em$z.coef[1:(nx+1)], B = Bem)
        tau1.res = rbind(tau1.res,
                         data.frame(zetaz = zetaZ[i], zetat1 = zetaT[j], zetat2 = zetat2,
                                    tau1 = tau.em.final$t1.coef[nx+1], tau1.se = tau.em.final$t1.coef.se[nx+1],
                                    pR2z = mean(partialR2z), pR2t1 = tau.em.final$pR2t1, pR2t2 = tau.em.final$pR2t2))
        tau2.res = rbind(tau2.res,
                         data.frame(zetaz = zetaZ[i], zetat1 = zetaT[j], zetat2 = zetat2,
                                    tau2 = tau.em.final$t2.coef[nx+1], tau2.se = tau.em.final$t2.coef.se[nx+1],
                                    pR2z = mean(partialR2z), pR2t1 = tau.em.final$pR2t1, pR2t2 = tau.em.final$pR2t2))
      }
      else{
        print("No method found.")
        return(NA)
      }
    }
  }

  tau1.res$t = tau1.res$tau1/tau1.res$tau1.se
  tau2.res$t = tau2.res$tau2/tau2.res$tau2.se

  return (list(tau1 = tau1.res, tau2 = tau2.res))
}
