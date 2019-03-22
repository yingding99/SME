#-----------------------------------------#
# Weibull_Median_2grp.R includes 4 functions: 
#		1. median.surv		  
#		2. var.median.surv.dt.dim6
#		3. ratio.median.surv.dt.dim3
#		4. dif.median.surv.dt.dim3
#		5. log.ratio.median.surv.dt.dim3
#-----------------------------------------#


library(rootSolve) # need this package for the median.surv function
library(survival)
library(eha)


# calculate the median survival time for each subgroup and the combined groups
median.surv <- function(lambda,k,theta,p.neg){
  
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

var.median.surv.dt.dim6 <- function(med.C.neg, med.C.pos, med.Rx.neg, med.Rx.pos, med.C, med.Rx,p.neg,lambda,k,theta,param.var){
	
	# data have the following columns: Y, Delta, Trt, M, Trt_M
	
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


# Use ratio (med.Rx.neg/med.C.neg, med.Rx.pos/med.C.pos, med.Rx/med.C) as the efficacy measure

ratio.median.surv.dt.dim3 <- function(lambda,k,theta,p.neg,param.var){
	
	median.comp <- median.surv(lambda,k,theta,p.neg) # compute the median survival time for the combined groups
	
	median.C.neg <- median.comp$median.C.neg
	median.C.pos <- median.comp$median.C.pos
	median.Rx.neg <- median.comp$median.Rx.neg
	median.Rx.pos <- median.comp$median.Rx.pos
	median.C <- median.comp$median.C
	median.Rx <- median.comp$median.Rx
	
	gamma.neg <- median.Rx.neg/median.C.neg
	gamma.pos <- median.Rx.pos/median.C.pos
	gamma.mix <- median.Rx/median.C
	
	median.var.dim6 <- var.median.surv.dt.dim6(med.C.neg=median.C.neg, med.C.pos=median.C.pos,
			med.Rx.neg=median.Rx.neg, med.Rx.pos=median.Rx.pos, med.C=median.C, med.Rx=median.Rx,
			p.neg,lambda,k,theta,param.var)  # compute the variance estimates for the median survival time
	
	var.median <- median.var.dim6$median.var  # 6x6 matrix
	D.gamma <- t(matrix(c(-median.Rx.neg/(median.C.neg)^2, 0, 1/median.C.neg, 0, 0, 0,
			0, -median.Rx.pos/(median.C.pos)^2, 0, 1/median.C.pos, 0, 0, 
			0, 0, 0, 0,	-median.Rx/(median.C)^2, 1/median.C),3, 6, byrow =TRUE))
	var.gamma <- t(D.gamma) %*% var.median %*% (D.gamma)  
	
	return(list(ratio.neg=gamma.neg, ratio.pos=gamma.pos, ratio.mix=gamma.mix, var.ratio=var.gamma))
}

# Use difference (med.Rx.neg-med.C.neg, med.Rx.pos-med.C.pos, med.Rx-med.C) as the efficacy measure

dif.median.surv.dt.dim3 <- function(lambda,k,theta,p.neg,param.var){
	
	median.comp <- median.surv(lambda,k,theta,p.neg) # compute the median survival time for the combined groups
	
	median.C.neg <- median.comp$median.C.neg
	median.C.pos <- median.comp$median.C.pos
	median.Rx.neg <- median.comp$median.Rx.neg
	median.Rx.pos <- median.comp$median.Rx.pos
	median.C <- median.comp$median.C
	median.Rx <- median.comp$median.Rx
	
	delta.neg <- median.Rx.neg-median.C.neg
	delta.pos <- median.Rx.pos-median.C.pos
	delta.mix <- median.Rx-median.C
	
	median.var.dim6 <- var.median.surv.dt.dim6(med.C.neg=median.C.neg, med.C.pos=median.C.pos,
			med.Rx.neg=median.Rx.neg, med.Rx.pos=median.Rx.pos, med.C=median.C, med.Rx=median.Rx,
			p.neg,lambda,k,theta,param.var)  # compute the variance estimates for the median survival time
	
	var.median <- median.var.dim6$median.var  # 6x6 matrix
	D.delta <- t(matrix(c(-1, 0, 1, 0, 0, 0,
							0, -1, 0, 1, 0, 0, 
							0, 0, 0, 0,	-1, 1),3, 6, byrow =TRUE))
	var.delta <- t(D.delta) %*% var.median %*% (D.delta)  
	
	return(list(dif.neg=delta.neg, dif.pos=delta.pos, dif.mix=delta.mix, var.dif=var.delta))
}

# Use log ratio (med.Rx.neg/med.C.neg, med.Rx.pos/med.C.pos, med.Rx/med.C) as the efficacy measure

log.ratio.median.surv.dt.dim3 <- function(lambda,k,theta,p.neg,param.var){
	
	median.comp <- median.surv(lambda,k,theta,p.neg) # compute the median survival time for the combined groups
	
	median.C.neg <- median.comp$median.C.neg
	median.C.pos <- median.comp$median.C.pos
	median.Rx.neg <- median.comp$median.Rx.neg
	median.Rx.pos <- median.comp$median.Rx.pos
	median.C <- median.comp$median.C
	median.Rx <- median.comp$median.Rx
	
	log.gamma.neg <- log(median.Rx.neg/median.C.neg)
	log.gamma.pos <- log(median.Rx.pos/median.C.pos)
	log.gamma.mix <- log(median.Rx/median.C)
	
	median.var.dim6 <- var.median.surv.dt.dim6(med.C.neg=median.C.neg, med.C.pos=median.C.pos,
			med.Rx.neg=median.Rx.neg, med.Rx.pos=median.Rx.pos, med.C=median.C, med.Rx=median.Rx,
			p.neg,lambda,k,theta,param.var)  # compute the variance estimates for the median survival time
	
	var.median <- median.var.dim6$median.var  # 6x6 matrix
	D.gamma <- t(matrix(c(-1/median.C.neg, 0, 1/median.Rx.neg, 0, 0, 0,
							0, -1/median.C.pos, 0, 1/median.Rx.pos, 0, 0, 
							0, 0, 0, 0,	-1/median.C, 1/median.Rx),3, 6, byrow =TRUE))
	var.gamma <- t(D.gamma) %*% var.median %*% (D.gamma)  
	
	return(list(log.ratio.neg=log.gamma.neg, log.ratio.pos=log.gamma.pos, log.ratio.mix=log.gamma.mix, var.ratio=var.gamma))
}


