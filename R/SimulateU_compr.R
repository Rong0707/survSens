SimulateU_compr <- function(t, d, z, x, zetat, zetaz, theta = 0.5, iter = 20, weights = NULL, offset = TRUE){
  #t is a vector of n, observed time
  #d is matrix n*2, d[,1] for event of interest, d[,2] for competing risk
  #z is a vector of n, indicator of treatment
  #x is matrix, observed covariates
  #zetat is a vector of 2, zetaz is a scaler, sensitivity parameters

  x = data.matrix(x)
  n = length(t) #number of observations
  nx = dim(x)[2] #number of observed covariates

  p = theta
  U = rbinom(n,1,p)
  #Upath = U

  for(j in 1:iter){
    if(offset){
      #fit time to event of interest
      U.fit1 = coxph(Surv(t,d[,1]) ~ z + x + offset(U*zetat[1]), weights = weights)
      t.coef1 = c(U.fit1$coef, zetat[1])
      #fit time to competing risk
      U.fit2 = coxph(Surv(t,d[,2]) ~ z + x + offset(U*zetat[2]), weights = weights)
      t.coef2 = c(U.fit2$coef, zetat[2])
      #fit treatment
      z.coef = c(glm(z ~ x, family = binomial(link = "probit"), offset = zetaz*U, control = glm.control(epsilon = 1e-6, maxit = 50))$coef, zetaz)
    }
    else{
      #fit time to event of interest
      U.fit1 = coxph(Surv(t,d[,1]) ~ z + x + U, weights = weights)
      t.coef1 = U.fit1$coef
      t.coef1[length(t.coef1)]  = zetat[1]
      #fit time to competing risk
      U.fit2 = coxph(Surv(t,d[,2]) ~ z + x + U, weights = weights)
      t.coef2 = U.fit2$coef
      t.coef2[length(t.coef2)]  = zetat[2]
      #fit treatment
      z.coef = glm(z ~ x + U, family = binomial(link = "probit"), control = glm.control(epsilon = 1e-6, maxit = 50))$coef
      z.coef[length(z.coef)] = zetaz
    }

    t.coef1[is.na(t.coef1)] = 0
    t.coef2[is.na(t.coef2)] = 0
    z.coef[is.na(z.coef)] = 0

    #cumulative baseline hazard for event of interest
    bh1 = basehaz(U.fit1, centered = F)
    index1 = match(t, bh1$time)

    #cumulative baseline hazard for competing risk
    bh2 = basehaz(U.fit2, centered = F)
    index2 = match(t, bh2$time)

    if(offset){
      ptzu1 = (1-pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1)))^(1-z)*
        pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1))^z*theta*
        exp(cbind(z,x,1)%*%matrix(t.coef1, ncol = 1))^d[,1] * exp(-bh1[index1,1]/exp(mean(U*zetat[1])) * exp(cbind(z,x,1)%*%matrix(t.coef1, ncol = 1)))*
        exp(cbind(z,x,1)%*%matrix(t.coef2, ncol = 1))^d[,2] * exp(-bh2[index2,1]/exp(mean(U*zetat[2])) * exp(cbind(z,x,1)%*%matrix(t.coef2, ncol = 1)))

      ptzu0 = (1-pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1)))^(1-z)*
        pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1))^z*(1-theta)*
        exp(cbind(z,x,0)%*%matrix(t.coef1, ncol = 1))^d[,1] * exp(-bh1[index1,1]/exp(mean(U*zetat[1])) * exp(cbind(z,x,0)%*%matrix(t.coef1, ncol = 1)))*
        exp(cbind(z,x,0)%*%matrix(t.coef2, ncol = 1))^d[,2] * exp(-bh2[index2,1]/exp(mean(U*zetat[2])) * exp(cbind(z,x,0)%*%matrix(t.coef2, ncol = 1)))
      }
    else{
      ptzu1 = (1-pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1)))^(1-z)*
        pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1))^z*theta*
        exp(cbind(z,x,1)%*%matrix(t.coef1, ncol = 1))^d[,1] * exp(-bh1[index1,1] * exp(cbind(z,x,1)%*%matrix(t.coef1, ncol = 1)))*
        exp(cbind(z,x,1)%*%matrix(t.coef2, ncol = 1))^d[,2] * exp(-bh2[index2,1] * exp(cbind(z,x,1)%*%matrix(t.coef2, ncol = 1)))

      ptzu0 = (1-pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1)))^(1-z)*
        pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1))^z*(1-theta)*
        exp(cbind(z,x,0)%*%matrix(t.coef1, ncol = 1))^d[,1] * exp(-bh1[index1,1] * exp(cbind(z,x,0)%*%matrix(t.coef1, ncol = 1)))*
        exp(cbind(z,x,0)%*%matrix(t.coef2, ncol = 1))^d[,2] * exp(-bh2[index2,1] * exp(cbind(z,x,0)%*%matrix(t.coef2, ncol = 1)))
    }

    p = ptzu1/(ptzu1 + ptzu0)
    p[ptzu1==0 & ptzu0==0] = 0

    U = rbinom(n,1,p)
    #Upath = cbind(Upath, U)
  }

  return(list(U = U,p = p))
}
