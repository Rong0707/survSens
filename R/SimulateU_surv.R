SimulateU_surv <- function(t, d, z, x, zetat, zetaz, theta = 0.5, iter = 20, weights = NULL, offset = TRUE){
  #t is a vector of n, time to event
  #d is a vector of n, indicator of event
  #z is a vector of n, treatment
  #x is a matrix, covariates
  #zetat is a scaler, zetaz is a scaler, sensitivity parameters

  n = length(t) #number of observations
  x = data.matrix(x)
  nx = dim(x)[2] #number of observed covariates

  p = theta
  U = rbinom(n,1,p)
  #Upath = U (useful for checking convergence)

  for(j in 1:iter){
    if(offset){
      #fit time to event
      U.fit1 = coxph(Surv(t,d) ~ z + x + offset(U*zetat), weights=weights)
      t.coef1 = c(U.fit1$coef, zetat)
      #fit treatment
      z.coef = c(glm(z ~ x, family=binomial(link="probit"), offset=zetaz*U, control = glm.control(epsilon = 1e-6, maxit = 50))$coef, zetaz)
    }
    else{
      #fit time to event
      U.fit1 = coxph(Surv(t,d) ~ z + x + U, weights=weights)
      t.coef1 = U.fit1$coef
      t.coef1[length(t.coef1)]  = zetat
      #fit treatment
      z.coef = glm(z ~ x + U, family = binomial(link="probit"), control = glm.control(epsilon = 1e-6, maxit = 50))$coef
      z.coef[length(z.coef)] = zetaz
    }

    t.coef1[is.na(t.coef1)] = 0
    z.coef[is.na(z.coef)] = 0

    bh1 = basehaz(U.fit1, centered=F) #cumulative baseline hazard for remission
    index1 = match(t,bh1$time)

    if(offset){
      ptzu1 = (1-pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1)))^(1-z)*
        pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1))^z*theta*
        exp(cbind(z,x,1)%*%matrix(t.coef1, ncol = 1))^d * exp(-bh1[index1,1]/exp(mean(U*zetat)) * exp(cbind(z,x,1)%*%matrix(t.coef1, ncol = 1)))

      ptzu0 = (1-pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1)))^(1-z)*
        pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1))^z*(1-theta)*
        exp(cbind(z,x,0)%*%matrix(t.coef1, ncol = 1))^d * exp(-bh1[index1,1]/exp(mean(U*zetat)) * exp(cbind(z,x,0)%*%matrix(t.coef1, ncol = 1)))
    }
    else{
      ptzu1 = (1-pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1)))^(1-z)*
        pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1))^z*theta*
        exp(cbind(z,x,1)%*%matrix(t.coef1, ncol = 1))^d * exp(-bh1[index1,1] * exp(cbind(z,x,1)%*%matrix(t.coef1, ncol = 1)))

      ptzu0 = (1-pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1)))^(1-z)*
        pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1))^z*(1-theta)*
        exp(cbind(z,x,0)%*%matrix(t.coef1, ncol = 1))^d * exp(-bh1[index1,1] * exp(cbind(z,x,0)%*%matrix(t.coef1, ncol = 1)))
    }

    p = ptzu1/(ptzu1 + ptzu0)
    p[ptzu1==0 & ptzu0==0] = 0

    U = rbinom(n,1,p)
    #Upath = cbind(Upath, U)
  }

  return(list(U = U, p = p))
}
