survSensitivity <- function(t, d, Z, X, method, zetaT = seq(-2,2,by=0.5), zetaZ = seq(-2,2,by=0.5), theta = 0.5, B = 50, Bem = 200){
  n = length(t)
  data = list(t = t, d = d, Z = Z, X = X)
  data1 = data.frame(data)

  tau1.res = data.frame()

  Z.fit = glm(Z ~ X, family = binomial(link="probit"))
  Coefz.noU = Z.fit$coefficients

  if(method == "stoEM_reg" | method == "EM_reg"){
    t1.fit = coxph(Surv(t, d) ~ X + Z)
    Coeft1.noU = t1.fit$coefficients
    }
  else if(method == "stoEM_IPW"){
    ps = Z.fit$fitted.values
    ipw = (sum(Z)/n) * Z/ps + (1-(sum(Z)/n))*(1-Z)/(1-ps)
    ipw = pmin(ipw, 10)
    ipw = pmax(ipw, 0.1)
    t1.ipw = coxph(Surv(t, d) ~ Z, weights = ipw, robust = TRUE)
    Coeft1.noU = t1.ipw$coefficients
  }
  else{
    print("No method found.")
    return(NA)
  }

  for(i in 1:length(zetaZ)){
    for(j in 1:length(zetaT)){
      if(method == "stoEM_reg"){
        #stochastic EM/regression
        tau.sto = surv_stoEM_regression(data, zetat = zetaT[j], zetaz = zetaZ[i], theta = theta, B = B)
        tau1.res = rbind(tau1.res, data.frame(zetaz = zetaZ[i], zetat = zetaT[j], tau1 = tau.sto$tau1, tau1.se = tau.sto$tau1.se))}
      else if(method == "stoEM_IPW"){
        #stochastic EM/IPW
        tau.ipw = surv_stoEM_ipw(data, zetat = zetaT[j], zetaz = zetaZ[i], theta = theta, B = B)
        tau1.res = rbind(tau1.res, data.frame(zetaz = zetaZ[i], zetat = zetaT[j], tau1 = tau.ipw$tau1, tau1.se = tau.ipw$tau1.se))}
      else if(method == "EM_reg"){
        #EM/regression
        tau.em = emU_surv(t, d, Z, X, zetat = zetaT[j], zetaz = zetaZ[i], theta = theta)
        data1$p = tau.em$p
        nx = length(tau.em$z.coef) - 2
        tau.em.final = surv_EM_variance(data1, zetat = zetaT[j], zetaz = zetaZ[i], z.coef = tau.em$z.coef[1:(nx+1)], B = Bem)
        tau1.res = rbind(tau1.res, data.frame(zetaz = zetaZ[i], zetat = zetaT[j], tau1 = tau.em.final$coef[nx+1], tau1.se = tau.em.final$coef.se[nx+1]))}
      else{
        print("No method found.")
        return(NA)
      }
    }
  }

  tau1.res$t = tau1.res$tau1/tau1.res$tau1.se

  return (tau1.res)
}
