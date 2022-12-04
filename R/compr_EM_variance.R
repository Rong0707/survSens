compr_EM_variance <- function(data, zetat, zetaz, z.coef, B = 1000){
  #data has to be (t,d,z,x,p)
  n = dim(data)[1]
  nx = dim(data)[2] - 4 #number of observed covariates
  data$d1 = as.numeric(data$d == 1)
  data$d2 = as.numeric(data$d == 2)

  data = data[order(data$t),]
  data$u1 = log(exp(zetat[1])*data$p+(1-data$p)) #u1 is the offset for t1
  data$u2 = log(exp(zetat[2])*data$p+(1-data$p)) #u2 is the offset for t2

  X1 = cbind(1,data.matrix(data[,4:(3+nx)]))
  XZ = data.matrix(data[,c(4:(3+nx),3)]) #put X and Z together

  # CoxPH model for t1
  U.fit1 = coxph(Surv(t,d1) ~ XZ + offset(u1), data = data)
  bh1 = basehaz(U.fit1, centered=F)
  K1 = dim(bh1)[1] #number of distinct time for t1
  bh1[,1] = bh1[,1]/exp(mean(data$u1))
  bh1$lambda = bh1[,1] - c(0,bh1[1:(K1-1),1])

  # CoxPH model for t2
  U.fit2 = coxph(Surv(t,d2) ~ XZ + offset(u2), data = data)
  bh2 = basehaz(U.fit2, centered=F)
  K2 = dim(bh2)[1] #number of distinct time for t2
  bh2[,1] = bh2[,1]/exp(mean(data$u2))
  bh2$lambda = bh2[,1] - c(0,bh2[1:(K2-1),1])

  data$t1.lin = XZ %*% U.fit1$coefficients + data$u1
  data$t2.lin = XZ %*% U.fit2$coefficients + data$u2
  data$z.lin0 = X1 %*% z.coef
  data$z.lin1 = X1 %*% z.coef + zetaz

  #codes for lambda when there are no ties
  #bh1$lambda1 = numeric(K1)
  bh1$event = numeric(K1)
  for(i in 1:K1){
    #bh1$lambda1[i] = 1/sum((data$t>=bh1$time[i]) * exp(data$t1.lin)) * sum(data$d1[data$t==bh1$time[i]])
    bh1$event[i] = sum(data$d1[data$t==bh1$time[i]])
  }

  #bh2$lambda1 = numeric(K2)
  bh2$event = numeric(K2)
  for(i in 1:K2){
    #bh2$lambda1[i] = 1/sum((data$t>=bh2$time[i]) * exp(data$t2.lin)) * sum(data$d2[data$t==bh2$time[i]])
    bh2$event[i] = sum(data$d2[data$t==bh2$time[i]])
  }

  index1 = match(data$t, bh1[,2])
  data$Lambda1 = bh1[index1,1]

  index2 = match(data$t, bh2[,2])
  data$Lambda2 = bh2[index2,1]

  dim_t1 = nx + 1 + sum(bh1$lambda!=0)
  dim_t2 = nx + 1 + sum(bh2$lambda!=0)
  dim_z = nx + 1
  dim_total =  dim_t1 + dim_t2 + dim_z
  #total dim is the sum of number of coefs in t1, t2 and z.
  SS = matrix(0, ncol = dim_total, nrow = dim_total)
  for (b in 1:B){
    SS = SS + ss_compr(data, zetat, zetaz, U.fit1, U.fit2, z.coef, nx, XZ, X1, bh1, bh2)
  }
  SS = SS/B

  #matrix of second derivative consists of 3 blocks for t1, t2, z
  l2 = matrix(0, ncol = dim_total, nrow = dim_total)

  #block for t1 is from 1:dim_t1
  XZw1 = XZ * sqrt(data$Lambda1 * exp(data$t1.lin))[,1]
  l2[1:(nx+1),1:(nx+1)] = t(XZw1) %*% XZw1 #block for (beta, zeta) for t1

  l2[(nx+2):dim_t1,(nx+2):dim_t1] = diag((bh1$event/bh1$lambda^2)[bh1$lambda!=0]) #the block for lambda

  l2_lam1 = matrix(0, ncol = nx+1, nrow = K1)
  for (i in 1:K1){
    if(bh1$lambda[i]!=0) l2_lam1[i,] = t(XZ*(data$t>=bh1$time[i])) %*% (exp(data$t1.lin))
  }

  l2[(nx+2):dim_t1, 1:(nx+1)] = l2_lam1[bh1$lambda!=0,]
  l2[1:(nx+1), (nx+2):dim_t1] = t(l2_lam1[bh1$lambda!=0,])

  #block for t2 (dim_t1 + 1): (dim_t1+dim_t2)
  XZw2 = XZ * sqrt(data$Lambda2 * exp(data$t2.lin))[,1]
  l2[(dim_t1+1):(dim_t1+nx+1),(dim_t1+1):(dim_t1+nx+1)] = t(XZw2) %*% XZw2 #block for (beta, zeta) for t2

  l2[(dim_t1+nx+2):(dim_t1+dim_t2),(dim_t1+nx+2):(dim_t1+dim_t2)] = diag((bh2$event/bh2$lambda^2)[bh2$lambda!=0]) #the block for lambda

  l2_lam2 = matrix(0, ncol = nx+1, nrow = K2)
  for (i in 1:K2){
    if(bh2$lambda[i]!=0) l2_lam2[i,] = t(XZ*(data$t>=bh2$time[i])) %*% (exp(data$t2.lin))
  }

  l2[(dim_t1+nx+2):(dim_t1+dim_t2),(dim_t1+1):(dim_t1+nx+1)] = l2_lam2[bh2$lambda!=0,]
  l2[(dim_t1+1):(dim_t1+nx+1), (dim_t1+nx+2):(dim_t1+dim_t2)] = t(l2_lam2[bh2$lambda!=0,])

  #block for z
  X1w = X1 * (sqrt(dnorm(data$z.lin0) * (data$Z * (dnorm(data$z.lin0) + data$z.lin0 * pnorm(data$z.lin0))/pnorm(data$z.lin0)^2 +
                  (1-data$Z)*(dnorm(data$z.lin0) - data$z.lin0 *(1- pnorm(data$z.lin0)))/(1- pnorm(data$z.lin0))^2) * (1-data$p) +
              dnorm(data$z.lin1) * (data$Z * (dnorm(data$z.lin1) + data$z.lin1 * pnorm(data$z.lin1))/pnorm(data$z.lin1)^2 +
                  (1-data$Z)*(dnorm(data$z.lin1) - data$z.lin1 *(1- pnorm(data$z.lin1)))/(1- pnorm(data$z.lin1))^2) * data$p))[,1]

  l2[(dim_t1 + dim_t2 + 1):dim_total,(dim_t1 + dim_t2 + 1):dim_total] = t(X1w) %*% X1w

  #used to check l/beta
  #test = matrix(0, ncol=3, nrow=3)
  #for(i in 1:n){
  #  test = test + X1[i,] %*% matrix(X1[i,], nrow=1)  * (dnorm(data$z.lin0[i]) * (data$Z[i] * (dnorm(data$z.lin0[i]) + data$z.lin0[i] * pnorm(data$z.lin0[i]))/pnorm(data$z.lin0[i])^2 +
  #                                                 (1-data$Z[i])*(dnorm(data$z.lin0[i]) - data$z.lin0[i] *(1- pnorm(data$z.lin0[i])))/(1- pnorm(data$z.lin0[i]))^2) * (1-data$p[i]) +
  #    dnorm(data$z.lin1[i]) * (data$Z[i] * (dnorm(data$z.lin1[i]) + data$z.lin1[i] * pnorm(data$z.lin1[i]))/pnorm(data$z.lin1[i])^2 +
  #                            (1-data$Z[i])*(dnorm(data$z.lin1[i]) - data$z.lin1[i] *(1- pnorm(data$z.lin1[i])))/(1- pnorm(data$z.lin1[i]))^2) * data$p[i])
  #}

  I_theta = l2 - SS
  coef.se = sqrt(diag(solve(I_theta)))
  return (list(t1.coef = U.fit1$coeff, t1.coef.se=coef.se[1:(nx + 1)],
               t2.coef = U.fit2$coeff, t2.coef.se=coef.se[(dim_t1 + 1):(dim_t1 + nx + 1)]))
}

ss_compr <- function(data, zetat, zetaz, U.fit1, U.fit2, z.coef, nx, XZ, X1, bh1, bh2){
  n = dim(data)[1]
  U = rbinom(n,1,data$p) #Sample U

  data$t1.lin = XZ %*% U.fit1$coefficients + zetat[1] * U
  data$t2.lin = XZ %*% U.fit2$coefficients + zetat[2] * U
  data$z.lin = X1 %*% z.coef + zetaz * U

  #vector of first derivative of t1
  s1 = numeric(nx + 1 + sum(bh1$lambda!=0))
  s1[1:(nx+1)] = t(XZ) %*% (data$d1 - data$Lambda1 * exp(data$t1.lin))

  K1 = dim(bh1)[1]
  s1_lambda = numeric(K1)
  for(i in 1:K1){
    if(bh1$lambda[i]!=0) s1_lambda[i] = bh1$event[i]/bh1$lambda[i] - sum((data$t>=bh1$time[i]) * exp(data$t1.lin))
  }
  s1[(nx+2):(nx + 1 + sum(bh1$lambda!=0))] = s1_lambda[bh1$lambda!=0]

  #vector of first derivative of t2
  s2 = numeric(nx + 1 + sum(bh2$lambda!=0))
  s2[1:(nx+1)] = t(XZ) %*% (data$d2 - data$Lambda2 * exp(data$t2.lin))
  K2 = dim(bh2)[1]
  s2_lambda = numeric(K2)
  for(i in 1:K2){
    if(bh2$lambda[i]!=0) s2_lambda[i] = bh2$event[i]/bh2$lambda[i] - sum((data$t>=bh2$time[i]) * exp(data$t2.lin))
  }
  s2[(nx+2):(nx + 1 + sum(bh2$lambda!=0))] = s2_lambda[bh2$lambda!=0]

  #vector of first derivative of z
  sz = t(X1) %*% (data$Z * dnorm(data$z.lin)/pnorm(data$z.lin) - (1-data$Z) * dnorm(data$z.lin)/(1-pnorm(data$z.lin)))

  #sz = c(0,0,0)
  #for(i in 1:n){
  #  sz = sz + (data$Z[i] * dnorm(data$z.lin[i])/pnorm(data$z.lin[i]) - (1-data$Z[i]) * dnorm(data$z.lin[i])/(1-pnorm(data$z.lin[i]))) * X1[i,]
  #}

  s = c(s1,s2,sz)
  ss = matrix(s, ncol=1) %*% matrix(s, nrow=1)
  return (ss)
}
