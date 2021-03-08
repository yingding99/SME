#-----------------------------------------#
# Weibull_Median_RatDiff_2grp.R includes 4 functions: 
#  	1. weib: weibull fitting for 2 group case
#   2. median.surv: calculate median survival for each subgroup and combined group
#  	3. var.median.surv.dt.dim6
#  	4. log.ratio.median.surv.dt.dim3
#   5. log.ratio.simu.CI: CI computation is based on log-scale ratio 
#	  6. dif.median.surv.dt.dim3
#   7. dif.median.simu.CI
#
#-----------------------------------------#

weib <- function(formula, data){
  # data contain the following columns: Y, Delta, Trt, M
  # Trt=1: treatment; Trt=0: control
  fit.weib <- weibreg(formula=formula, data=data)
  #fit.weib <- weibreg(Surv(Y,Delta) ~ relevel(as.factor(Trt),"1")+as.factor(M)+relevel(as.factor(Trt),"1")*as.factor(M), data=data)
  fit.LogCoeff <- fit.weib$coef
  fit.Coeff <- exp(fit.LogCoeff) # the coefficient in the "original" scale (exponentiated from the log scale)
  # change the order of the parameters (so that they correspond to "lambda","k","theta1","theta2","theta3")
  fit.Coeff <- fit.Coeff[c(4,5,1:3)]
  names(fit.Coeff)[1:2]<-c("Scale","Shape")
  # change the order of the var parameters accordingly 
  fit.LogVar <- fit.weib$var
  tmp1 <- fit.LogVar[1:3,1:3]
  tmp2 <- fit.LogVar[4:5,4:5]
  tmp3 <- fit.LogVar[1:3,4:5]
  
  fit.LogVar.2 <- rbind(cbind(tmp2,t(tmp3)),cbind(tmp3,tmp1))
  # delta method to compute the Variance matrix for the parameter estimates in the "original" scale
  fit.Var <- diag(fit.Coeff) %*% fit.LogVar.2 %*% diag(fit.Coeff) 
  
  return(list(log.coef=fit.LogCoeff,log.var=fit.LogVar,coef=fit.Coeff,var=fit.Var))
}



# calculate the median survival time for each subgroup and the combined groups
median.surv <- function(lambda,k,theta,p.neg){
  
  # lambda: Weibull distribution's scale parameter 
  # k: Weibull distribution's shape parameter 
  # theta=(theta1,theta2,theta3): =exp(beta) where beta is the slope parameter in the Cox PH model
  # p.neg: Prob(g-) (assume Prob(g-|C)=Prob(g-|Rx))
  
  median.C.neg <- as.numeric(lambda*(log(2))^(1/k))
  median.C.pos <- as.numeric(lambda*(log(2)/theta[2])^(1/k))
  median.Rx.neg <- as.numeric(lambda*(log(2)/theta[1])^(1/k))  
  median.Rx.pos <- as.numeric(lambda*(log(2)/(theta[1]*theta[2]*theta[3]))^(1/k))
  
  f.C <- function(t){
    f=p.neg*exp(-(t/lambda)^k)+(1-p.neg)*exp(-theta[2]*(t/lambda)^k)-0.5
    return(f)
  }
  f.Rx <- function(t){
    f=p.neg*exp(-theta[1]*(t/lambda)^k)+(1-p.neg)*exp(-theta[1]*theta[2]*theta[3]*(t/lambda)^k)-0.5
    return(f)
  }
  median.C <- uniroot(f.C,c(max(median.C.neg,median.C.pos),min(median.C.neg,median.C.pos)))$root
  median.Rx <- uniroot(f.Rx,c(max(median.Rx.neg,median.Rx.pos),min(median.Rx.neg,median.Rx.pos)))$root
  
  return(list(median.C.neg=median.C.neg,median.C.pos=median.C.pos,
              median.Rx.neg=median.Rx.neg,median.Rx.pos=median.Rx.pos,
              median.C=median.C,median.Rx=median.Rx))
  
}

# use the delta method for implicitly defined RVs
# reference paper: Jacques Benichou and Mitchell Gail (The American Statistician, 1989, Vol.43, No.1) 

var.median.surv.dt.dim6 <- function(med.C.neg, med.C.pos, med.Rx.neg, med.Rx.pos, med.C, med.Rx,
                                    p.neg,lambda,k,theta,param.var){
  

  # calculate the variance estimate for the median survival time estimate of all group
  # using delta method for implicitly defined random variables
  
  J.11 <- (-k*exp(-(med.C.neg/lambda)^k)*(med.C.neg/lambda)^(k-1))/lambda 
  J.22 <- (-k*theta[2]*exp(-theta[2]*(med.C.pos/lambda)^k)*(med.C.pos/lambda)^(k-1))/lambda 
  J.33 <- (-k*theta[1]*exp(-theta[1]*(med.Rx.neg/lambda)^k)*(med.Rx.neg/lambda)^(k-1))/lambda 
  J.44 <- (-k*theta[1]*theta[2]*theta[3]*exp(-theta[1]*theta[2]*theta[3]*(med.Rx.pos/lambda)^k)*(med.Rx.pos/lambda)^(k-1))/lambda 
  J.55 <- (-p.neg*k*exp(-(med.C/lambda)^k)*(med.C/lambda)^(k-1)
           -(1-p.neg)*k*theta[2]*exp(-theta[2]*(med.C/lambda)^k)*(med.C/lambda)^(k-1))/lambda 
  J.66 <- (-p.neg*k*theta[1]*exp(-theta[1]*(med.Rx/lambda)^k)*(med.Rx/lambda)^(k-1)
           -(1-p.neg)*k*theta[1]*theta[2]*theta[3]*exp(-theta[1]*theta[2]*theta[3]*(med.Rx/lambda)^k)*(med.Rx/lambda)^(k-1))/lambda 
  
  J.mat <- diag(c(J.11, J.22, J.33, J.44, J.55, J.66))  
  J.mat.inv <- solve(J.mat)
  
  H.mat <- mat.or.vec(6,5)
  H.mat[1,] <- c(exp(-(med.C.neg/lambda)^k)*k*(med.C.neg/lambda)^(k-1)*(med.C.neg/lambda^2), 
                 -exp(-(med.C.neg/lambda)^k)*(med.C.neg/lambda)^k*log(med.C.neg/lambda),
                 0,
                 0,
                 0)
  H.mat[2,] <- c(exp(-theta[2]*(med.C.pos/lambda)^k)*theta[2]*k*(med.C.pos/lambda)^(k-1)*(med.C.pos/lambda^2), 
                 -exp(-theta[2]*(med.C.pos/lambda)^k)*theta[2]*(med.C.pos/lambda)^k*log(med.C.pos/lambda),
                 0,
                 -exp(-theta[2]*(med.C.pos/lambda)^k)*(med.C.pos/lambda)^k,
                 0)
  H.mat[3,] <- c(exp(-theta[1]*(med.Rx.neg/lambda)^k)*theta[1]*k*(med.Rx.neg/lambda)^(k-1)*(med.Rx.neg/lambda^2),
                 -exp(-theta[1]*(med.Rx.neg/lambda)^k)*theta[1]*(med.Rx.neg/lambda)^k*log(med.Rx.neg/lambda),   
                 -exp(-theta[1]*(med.Rx.neg/lambda)^k)*(med.Rx.neg/lambda)^k,
                 0,
                 0)
  H.mat[4,] <- c(exp(-theta[1]*theta[2]*theta[3]*(med.Rx.pos/lambda)^k)*theta[1]*theta[2]*theta[3]*k*(med.Rx.pos/lambda)^(k-1)*(med.Rx.pos/lambda^2),
                 -exp(-theta[1]*theta[2]*theta[3]*(med.Rx.pos/lambda)^k)*theta[1]*theta[2]*theta[3]*(med.Rx.pos/lambda)^k*log(med.Rx.pos/lambda),   
                 -exp(-theta[1]*theta[2]*theta[3]*(med.Rx.pos/lambda)^k)*theta[2]*theta[3]*(med.Rx.pos/lambda)^k,
                 -exp(-theta[1]*theta[2]*theta[3]*(med.Rx.pos/lambda)^k)*theta[1]*theta[3]*(med.Rx.pos/lambda)^k,
                 -exp(-theta[1]*theta[2]*theta[3]*(med.Rx.pos/lambda)^k)*theta[1]*theta[2]*(med.Rx.pos/lambda)^k)
  
  H.mat[5,] <- c(p.neg*exp(-(med.C/lambda)^k)*k*(med.C/lambda)^(k-1)*(med.C/lambda^2)
                 +(1-p.neg)*exp(-theta[2]*(med.C/lambda)^k)*theta[2]*k*(med.C/lambda)^(k-1)*(med.C/lambda^2), 
                 -p.neg*exp(-(med.C/lambda)^k)*(med.C/lambda)^k*log(med.C/lambda)
                 -(1-p.neg)*exp(-theta[2]*(med.C/lambda)^k)*theta[2]*(med.C/lambda)^k*log(med.C/lambda),
                 0,
                 -(1-p.neg)*exp(-theta[2]*(med.C/lambda)^k)*(med.C/lambda)^k,
                 0)
  
  H.mat[6,] <- c(p.neg*exp(-theta[1]*(med.Rx/lambda)^k)*theta[1]*k*(med.Rx/lambda)^(k-1)*(med.Rx/lambda^2)
                 +(1-p.neg)*exp(-theta[1]*theta[2]*theta[3]*(med.Rx/lambda)^k)*theta[1]*theta[2]*theta[3]*k*(med.Rx/lambda)^(k-1)*(med.Rx/lambda^2),
                 -p.neg*exp(-theta[1]*(med.Rx/lambda)^k)*theta[1]*(med.Rx/lambda)^k*log(med.Rx/lambda) 
                 -(1-p.neg)*exp(-theta[1]*theta[2]*theta[3]*(med.Rx/lambda)^k)*theta[1]*theta[2]*theta[3]*(med.Rx/lambda)^k*log(med.Rx/lambda),   
                 -p.neg*exp(-theta[1]*(med.Rx/lambda)^k)*(med.Rx/lambda)^k-(1-p.neg)*exp(-theta[1]*theta[2]*theta[3]*(med.Rx/lambda)^k)*theta[2]*theta[3]*(med.Rx/lambda)^k,
                 -(1-p.neg)*exp(-theta[1]*theta[2]*theta[3]*(med.Rx/lambda)^k)*theta[1]*theta[3]*(med.Rx/lambda)^k,
                 -(1-p.neg)*exp(-theta[1]*theta[2]*theta[3]*(med.Rx/lambda)^k)*theta[1]*theta[2]*(med.Rx/lambda)^k)
  
  median.var <- J.mat.inv %*% H.mat %*% param.var %*% t(H.mat) %*% t(J.mat.inv)
  
  return(list(median.var = median.var))
}

#-------#
# ratio #
#-------#
# Use ratio (med.Rx.neg/med.C.neg, med.Rx.pos/med.C.pos, med.Rx/med.C) as the efficacy measure

log.ratio.median.surv.dt.dim3 <- function(median.surv,median.var){
  
  median.C.neg <- median.surv["median.C.neg"]
  median.C.pos <- median.surv["median.C.pos"]
  median.Rx.neg <- median.surv["median.Rx.neg"]
  median.Rx.pos <- median.surv["median.Rx.pos"]
  median.C <- median.surv["median.C"]
  median.Rx <- median.surv["median.Rx"]
  
  # ratio for M-, M+, combined
  gamma.neg <- median.Rx.neg/median.C.neg
  gamma.pos <- median.Rx.pos/median.C.pos
  gamma.mix <- median.Rx/median.C
  
  # log(ratio) for M-, M+, combined
  log.gamma.neg <- log(median.Rx.neg/median.C.neg)
  log.gamma.pos <- log(median.Rx.pos/median.C.pos)
  log.gamma.mix <- log(median.Rx/median.C)
  
  
  D.log.gamma <- t(matrix(c(-1/median.C.neg, 0, 1/median.Rx.neg, 0, 0, 0,
                            0, -1/median.C.pos, 0, 1/median.Rx.pos, 0, 0, 
                            0, 0, 0, 0,  -1/median.C, 1/median.Rx),3, 6, byrow =TRUE))
  # var.log.gamma is the variance/covariance matrix for log(ratio)
  var.log.gamma <- t(D.log.gamma) %*% median.var %*% (D.log.gamma)  
  
  return(list(log.ratio.neg=log.gamma.neg, log.ratio.pos=log.gamma.pos, log.ratio.mix=log.gamma.mix, 
              ratio.neg=gamma.neg, ratio.pos=gamma.pos, ratio.mix=gamma.mix,
              var.log.ratio=var.log.gamma))
}

# CI computation is based on log-scale ratio
log.ratio.simu.CI <- function(seed,alpha,log.ratio,var.log.ratio){
  # qmvnorm is simulation based, provide a seed
  set.seed(seed)
  log.CIs <- qmvnorm(1-alpha, corr=cov2cor(var.log.ratio), tail = "both")$quantile*sqrt(diag(var.log.ratio))
  
  ratio.CI.lbd <- exp(log.ratio - log.CIs)
  ratio.CI.ubd <- exp(log.ratio + log.CIs)
  
  return(list(ratio.CI.lbd=ratio.CI.lbd,ratio.CI.ubd=ratio.CI.ubd))
}



#------------#
# difference #
#------------#
# Use difference (med.Rx.neg-med.C.neg, med.Rx.pos-med.C.pos, med.Rx-med.C) as the efficacy measure

dif.median.surv.dt.dim3 <- function(median.surv, median.var){
	
	median.C.neg <- median.surv["median.C.neg"]
	median.C.pos <- median.surv["median.C.pos"]
	median.Rx.neg <- median.surv["median.Rx.neg"]
	median.Rx.pos <- median.surv["median.Rx.pos"]
	median.C <- median.surv["median.C"]
	median.Rx <- median.surv["median.Rx"]
	
	delta.neg <- median.Rx.neg-median.C.neg
	delta.pos <- median.Rx.pos-median.C.pos
	delta.mix <- median.Rx-median.C
	
	D.delta <- t(matrix(c(-1, 0, 1, 0, 0, 0,
							0, -1, 0, 1, 0, 0, 
							0, 0, 0, 0,	-1, 1),3, 6, byrow =TRUE))
	var.delta <- t(D.delta) %*% median.var %*% (D.delta)  
	
	return(list(dif.neg=delta.neg, dif.pos=delta.pos, dif.mix=delta.mix, var.dif=var.delta))
}


# CI computation is based on difference of median

dif.median.simu.CI <- function(seed, alpha, dif.median, var.dif.median){
  # qmvnorm is simulation based, provide a seed
  set.seed(seed)
  CIs <- qmvnorm(1-alpha, corr=cov2cor(var.dif.median), tail = "both")$quantile*sqrt(diag(var.dif.median))
  
  dif.median.CI.lbd <- dif.median - CIs
  dif.median.CI.ubd <- dif.median + CIs
  
  return(list(dif.median.CI.lbd=dif.median.CI.lbd, dif.median.CI.ubd=dif.median.CI.ubd))
}

