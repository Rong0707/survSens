surv_EM_variance <- function(data, zetat, zetaz, z.coef, B = 1000){
  #data has to be (t,d,z,x,p)
  n = dim(data)[1]
  nx = dim(data)[2] - 4 #number of observed covariates

  data = data[order(data$t),]
  data$u = log(exp(zetat)*data$p+(1-data$p)) #u is the offset

  #XZ = cbind(data$X.1, data$X.2, data$Z)
  X1 = cbind(1,data.matrix(data[,4:(3+nx)]))
  XZ = data.matrix(data[,c(4:(3+nx),3)]) #put X and Z together

  U.fit1 = coxph(Surv(t,d) ~ XZ + offset(u), data = data)
  # Calculate partial R-sq of (t, d) ~ u | x, z
  U.fit1_reduced = coxph(Surv(t,d) ~ XZ, data = data)
  logtest <- -2 * (U.fit1_reduced$loglik[2] - U.fit1$loglik[2])
  if (zetat >= 0)
    partialR2t1 = (1 - exp(-logtest/U.fit1$nevent))
  else
    partialR2t1 = - (1 - exp(-logtest/U.fit1$nevent))

  bh1 = basehaz(U.fit1, centered=F)
  K = dim(bh1)[1] #number of distinct time
  bh1[,1] = bh1[,1]/exp(mean(data$u))
  bh1$lambda = bh1[,1] - c(0,bh1[1:(K-1),1])

  data$t.lin = XZ %*% U.fit1$coefficients + data$u
  data$z.lin0 = X1 %*% z.coef
  data$z.lin1 = X1 %*% z.coef + zetaz

  #codes for checking lambda when there are no ties
  #bh1$lambda1 = numeric(K)
  bh1$event = numeric(K)
  for(i in 1:K){
    #bh1$lambda1[i] = 1/sum((data$t>=bh1$time[i]) * exp(data$t.lin)) * sum(data$d[data$t==bh1$time[i]])
    bh1$event[i] = sum(data$d[data$t==bh1$time[i]])
  }

  index = match(data$t, bh1[,2])
  data$Lambda = bh1[index,1] #cumulative baseline hazard

  dim_total = nx + 1 + sum(bh1$lambda!=0) #dim for t part
  SS = matrix(0, ncol = dim_total+(1+nx), nrow = dim_total+(1+nx))
  for (b in 1:B){
    SS = SS + ss(data, zetat, zetaz, U.fit1, z.coef, nx, XZ, X1, bh1)
  }
  SS = SS/B

  #codes for checking vectorized operation.
  #l2_beta = matrix(0, nrow = nx+1, ncol = nx+1)
  #for(i in 1:n){
  #  l2_beta = l2_beta + XZ[i,] %*% matrix(XZ[i,], nrow=1) * data$Lambda[i] * exp(data$t.lin[i])
  #}

  l2 = matrix(0, ncol = dim_total+(1+nx), nrow = dim_total+(1+nx)) #matrix of second derivative
  XZw = XZ * sqrt(data$Lambda * exp(data$t.lin))[,1]
  l2[1:(nx+1),1:(nx+1)] = t(XZw) %*% XZw #block for (beta, zeta)

  l2[(nx+2):dim_total,(nx+2):dim_total] = diag((bh1$event/bh1$lambda^2)[bh1$lambda!=0]) #the block for lambda

  l2_lam = matrix(0, ncol = nx+1, nrow = K)
  for (i in 1:K){
    if(bh1$lambda[i]!=0) l2_lam[i,] = t(XZ*(data$t>=bh1$time[i])) %*% (exp(data$t.lin))
  }

  l2[(nx+2):dim_total, 1:(nx+1)] = l2_lam[bh1$lambda!=0,]
  l2[1:(nx+1), (nx+2):dim_total] = t(l2_lam[bh1$lambda!=0,])

  X1w = X1 * (sqrt(dnorm(data$z.lin0) * (data$Z * (dnorm(data$z.lin0) + data$z.lin0 * pnorm(data$z.lin0))/pnorm(data$z.lin0)^2 +
                  (1-data$Z)*(dnorm(data$z.lin0) - data$z.lin0 *(1- pnorm(data$z.lin0)))/(1- pnorm(data$z.lin0))^2) * (1-data$p) +
              dnorm(data$z.lin1) * (data$Z * (dnorm(data$z.lin1) + data$z.lin1 * pnorm(data$z.lin1))/pnorm(data$z.lin1)^2 +
                  (1-data$Z)*(dnorm(data$z.lin1) - data$z.lin1 *(1- pnorm(data$z.lin1)))/(1- pnorm(data$z.lin1))^2) * data$p))[,1]

  l2[(dim_total+1):(dim_total+1+nx),(dim_total+1):(dim_total+1+nx)] = t(X1w) %*% X1w

  #non-vectorized codes used to check l/beta
  #test = matrix(0, ncol=nx+1, nrow=nx+1)
  #for(i in 1:n){
  #  test = test + X1[i,] %*% matrix(X1[i,], nrow=1)  * (dnorm(data$z.lin0[i]) * (data$Z[i] * (dnorm(data$z.lin0[i]) + data$z.lin0[i] * pnorm(data$z.lin0[i]))/pnorm(data$z.lin0[i])^2 +
  #                                                 (1-data$Z[i])*(dnorm(data$z.lin0[i]) - data$z.lin0[i] *(1- pnorm(data$z.lin0[i])))/(1- pnorm(data$z.lin0[i]))^2) * (1-data$p[i]) +
  #    dnorm(data$z.lin1[i]) * (data$Z[i] * (dnorm(data$z.lin1[i]) + data$z.lin1[i] * pnorm(data$z.lin1[i]))/pnorm(data$z.lin1[i])^2 +
  #                            (1-data$Z[i])*(dnorm(data$z.lin1[i]) - data$z.lin1[i] *(1- pnorm(data$z.lin1[i])))/(1- pnorm(data$z.lin1[i]))^2) * data$p[i])
  #}

  I_theta = l2 - SS
  coef.se = sqrt(diag(solve(I_theta)))
  return (list(coef = U.fit1$coeff, coef.se = coef.se[1:(nx+1)],
               pR2t1 = partialR2t1))
}

ss <- function(data, zetat, zetaz, U.fit1, z.coef, nx, XZ, X1, bh1){
  n = dim(data)[1]
  U = rbinom(n,1,data$p) #Sample U
  dim_total = nx+1+sum(bh1$lambda!=0)

  s1 = numeric(dim_total) #vector of first derivative
  data$t.lin = XZ %*% U.fit1$coefficients + zetat * U
  data$z.lin = X1 %*% z.coef + zetaz * U

  s1[1:(nx+1)] = t(XZ) %*% (data$d - data$Lambda * exp(data$t.lin))

  K = dim(bh1)[1]
  s_lambda = numeric(K)
  for(i in 1:K){
    if(bh1$lambda[i]!=0) s_lambda[i] = bh1$event[i]/bh1$lambda[i] - sum((data$t>=bh1$time[i]) * exp(data$t.lin))
  }
  s1[(nx+2):dim_total] = s_lambda[bh1$lambda!=0]

  sz = t(X1) %*% (data$Z * dnorm(data$z.lin)/pnorm(data$z.lin) - (1-data$Z) * dnorm(data$z.lin)/(1-pnorm(data$z.lin)))

  #sz = c(0,0,0)
  #for(i in 1:n){
  #  sz = sz + (data$Z[i] * dnorm(data$z.lin[i])/pnorm(data$z.lin[i]) - (1-data$Z[i]) * dnorm(data$z.lin[i])/(1-pnorm(data$z.lin[i]))) * X1[i,]
  #}

  s1 = c(s1,sz)
  ss = matrix(s1, ncol=1) %*% matrix(s1, nrow=1)
  return (ss)
}
