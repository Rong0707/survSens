emU_compr <- function(t, d, z, x, zetat, zetaz, theta = 0.5, iter = 20, weights = NULL){
  #t is a vector of n, observed time
  #d is matrix n*2, d[,1] for event of interest, d[,2] for competing risk
  #z is a vector of n, indicator of treatment
  #x is matrix, observed covariates
  #zetat is a vector of 2, zetaz is a scaler

  fn <- function(beta, x, z, p, zetaz) {
    -sum(z * (log(pnorm(x %*% beta + zetaz))*p + log(pnorm(x %*% beta))*(1-p)) +
          (1-z) * (log(1-pnorm(x %*% beta + zetaz))*p + log(1-pnorm(x %*% beta))*(1-p)))
  }

  x = data.matrix(x)
  n = length(t) #number of observations
  nx = dim(x)[2] #number of covariates

  #initialize parameters
  #fit time to event of interest
  U.fit1 = coxph(Surv(t,d[,1]) ~ z + x , weights = weights)
  t.coef1 = c(U.fit1$coef, zetat[1])
  #fit time to competing risk
  U.fit2 = coxph(Surv(t,d[,2]) ~ z + x, weights = weights)
  t.coef2 = c(U.fit2$coef, zetat[2])
  #fit treatment
  z.coef = c(glm(z ~ x, family=binomial(link="probit"), control = glm.control(epsilon = 1e-6, maxit = 50))$coef, zetaz)

  t.coef1[is.na(t.coef1)] = 0
  t.coef2[is.na(t.coef2)] = 0
  z.coef[is.na(z.coef)] = 0

  #used to check convergence
  #t.coef1.path =  t.coef1
  #t.coef2.path =  t.coef2

  p = rep(0,n)

  for(j in 1:iter){
    #cumulative baseline hazard for event of interest with mean(offset)
    bh1 = basehaz(U.fit1, centered=F)
    index1 = match(t,bh1$time)

    #cumulative baseline hazard for competing risk with mean(offset)
    bh2 = basehaz(U.fit2, centered=F)
    index2 = match(t,bh2$time)

    ptzu1 = (1-pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1)))^(1-z)*
      pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1))^z*theta*
      exp(cbind(z,x,1)%*%matrix(t.coef1, ncol = 1))^d[,1] * exp(-bh1[index1,1]/exp(mean(log(exp(zetat[1])*p+(1-p)))) * exp(cbind(z,x,1)%*%matrix(t.coef1, ncol = 1)))*
      exp(cbind(z,x,1)%*%matrix(t.coef2, ncol = 1))^d[,2] * exp(-bh2[index2,1]/exp(mean(log(exp(zetat[2])*p+(1-p)))) * exp(cbind(z,x,1)%*%matrix(t.coef2, ncol = 1)))

    ptzu0 = (1-pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1)))^(1-z)*
      pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1))^z*(1-theta)*
      exp(cbind(z,x,0)%*%matrix(t.coef1, ncol = 1))^d[,1] * exp(-bh1[index1,1]/exp(mean(log(exp(zetat[1])*p+(1-p)))) * exp(cbind(z,x,0)%*%matrix(t.coef1, ncol = 1)))*
      exp(cbind(z,x,0)%*%matrix(t.coef2, ncol = 1))^d[,2] * exp(-bh2[index2,1]/exp(mean(log(exp(zetat[2])*p+(1-p)))) * exp(cbind(z,x,0)%*%matrix(t.coef2, ncol = 1)))

    p = ptzu1/(ptzu1 + ptzu0)
    p[ptzu1==0 & ptzu0==0] = 0

    #fit time to event of interest
    U.fit1 = coxph(Surv(t,d[,1]) ~ z + x + offset(log(exp(zetat[1])*p+(1-p))), weights = weights)
    t.coef1 = c(U.fit1$coef, zetat[1])
    #fit time to competing risk
    U.fit2 = coxph(Surv(t,d[,2]) ~ z + x + offset(log(exp(zetat[2])*p+(1-p))), weights = weights)
    t.coef2 = c(U.fit2$coef, zetat[2])
    #fit treatment
    #z.coef = c(glm(z ~ x, family=binomial(link="probit"), offset=zetaz*p, control = glm.control(epsilon = 1e-6, maxit = 50))$coef, zetaz)
    #z.coef = c(optimize(f, lower=-2, upper=2, x=x, z=z, p=p, zetaz=zetaz, maximum = TRUE)$maximum, zetaz)
    z.fit = optim(par = z.coef[1:(nx+1)], fn, x=cbind(1,x), z=z, p=p, zetaz=zetaz)
    z.coef = c(z.fit$par, zetaz)

    t.coef1[is.na(t.coef1)] = 0
    t.coef2[is.na(t.coef2)] = 0
    z.coef[is.na(z.coef)] = 0

    #t.coef1.path = rbind(t.coef1.path, t.coef1)
    #t.coef2.path = rbind(t.coef2.path, t.coef2)
  }

  return(list(p = p, t.coef1 = t.coef1, t.coef2 = t.coef2, z.coef = z.coef))
}
