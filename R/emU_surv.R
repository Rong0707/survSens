emU_surv <- function(t, d, z, x, zetat, zetaz, theta = 0.5, iter = 20, weights = NULL){
  #t is a vector of n, time to event
  #d is a vector of n, indicator of event
  #z is a vector of n, treatment
  #x is a matrix, covariates
  #zetat is a scaler, zetaz is a scaler, sensitivity parameters

  fn <- function(beta, x, z, p, zetaz) {
    -sum(z * (log(pnorm(x %*% beta + zetaz))*p + log(pnorm(x %*% beta))*(1-p)) +
           (1-z) * (log(1-pnorm(x %*% beta + zetaz))*p + log(1-pnorm(x %*% beta))*(1-p)))
  }

  x = data.matrix(x)
  n = length(t) #number of observations
  nx = dim(x)[2] #number of covariates

  #initialize parameters
  #fit time to event
  U.fit1 = coxph(Surv(t,d) ~ z + x , weights=weights)
  t.coef1 = c(U.fit1$coef, zetat)
  #fit treatment
  z.coef = c(glm(z ~ x, family=binomial(link="probit"), control = glm.control(epsilon = 1e-6, maxit = 50))$coef, zetaz)

  t.coef1[is.na(t.coef1)] = 0
  z.coef[is.na(z.coef)] = 0

  p = rep(0,n)

  for(j in 1:iter){
    bh1 = basehaz(U.fit1, centered=F) #cumulative baseline hazard for event with mean(offset)
    index1 = match(t,bh1$time)

    ptzu1 = (1-pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1)))^(1-z)*
      pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1))^z*theta*
      exp(cbind(z,x,1)%*%matrix(t.coef1, ncol = 1))^d * exp(-bh1[index1,1]/exp(mean(log(exp(zetat)*p+(1-p)))) * exp(cbind(z,x,1)%*%matrix(t.coef1, ncol = 1)))

    ptzu0 = (1-pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1)))^(1-z)*
      pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1))^z*(1-theta)*
      exp(cbind(z,x,0)%*%matrix(t.coef1, ncol = 1))^d * exp(-bh1[index1,1]/exp(mean(log(exp(zetat)*p+(1-p)))) * exp(cbind(z,x,0)%*%matrix(t.coef1, ncol = 1)))

    p = ptzu1/(ptzu1 + ptzu0) #posterior dist of U
    p[ptzu1==0 & ptzu0==0] = 0

    #fit time to event
    U.fit1 = coxph(Surv(t,d) ~ z + x + offset(log(exp(zetat)*p+(1-p))), weights=weights)
    t.coef1 = c(U.fit1$coef, zetat[1])
    #fit treatment
    z.fit = optim(par = z.coef[1:(nx+1)], fn, x=cbind(1,x), z=z, p=p, zetaz=zetaz)
    z.coef = c(z.fit$par, zetaz)

    t.coef1[is.na(t.coef1)] = 0
    z.coef[is.na(z.coef)] = 0
  }

  return(list(p = p, t.coef1 = t.coef1, z.coef = z.coef))
}
