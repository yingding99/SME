#-----------------------------------------#
# Weibull_Fitting.R includes 4 functions: 
#		1. dataWeibull.unifCensor		  
#		2. dataWeibull.normalCensor
#		3. weib
#		4. surv.fun
#-----------------------------------------#


library(rootSolve) # need this package for the median.surv function
library(survival)
library(eha)

# "dataWeibull" generates the survival dataset

dataWeibull.unifCensor <- function(m,lambda,k,beta,n.obs,a,b){
  # assuming uniform random censoring
  # m: the total number of marker groups
  # lambda: Weibull distribution's scale parameter (for T)
  # k: Weibull distribution's shape parameter (for T)
  # a: lower bound of the uniform distribution (for C)
  # b: upper bound of the uniform distribution (for C)
  # beta: the slope parameter in the Cox PH model, a vector of length (2m-1)
  # n.obs: sample size vector (for each subgroup under each treatment arm), a vector of length 2m
  
  # the data set is sorted by Trt (first chunk is control, second chunk is treatment)
  # and then Marker (group 1 to m)
  Trt <- c(rep(0,sum(n.obs[1:m])),rep(1,sum(n.obs[(m+1):(2*m)])))
  M <- rep(0,sum(n.obs))  # M is the vector with values of 1 to m indicating the marker groups
  M.mat <- mat.or.vec(sum(n.obs),m) # M.mat is the matrix (jth column = 1 if M=j or 0)
  
  M[1:n.obs[1]] <- 1
  M[(sum(n.obs[1:m])+1) : sum(n.obs[1:(m+1)])] <- 1
  M.mat[1:n.obs[1],1] <- 1
  M.mat[(sum(n.obs[1:m])+1) : sum(n.obs[1:(m+1)]), 1] <- 1
  
  for(i in 2:m){
    M[(sum(n.obs[1:(i-1)])+1) : sum(n.obs[1:i])] <- i
    M[(sum(n.obs[1:(m+i-1)])+1) : sum(n.obs[1:(m+i)])] <- i
    
    M.mat[(sum(n.obs[1:(i-1)])+1) : sum(n.obs[1:i]), i] <- 1
    M.mat[(sum(n.obs[1:(m+i-1)])+1) : sum(n.obs[1:(m+i)]), i] <- 1
  }
  
  Trt_M <- Trt * M.mat  # this is the interaction matrix (dim=sum(n.obs)*m)
 
  T <- rep(0,sum(n.obs))
  C <- rep(0,sum(n.obs))
  Y <- rep(0,sum(n.obs))
  Delta <- rep(0,sum(n.obs))
  
  for(i in 1:sum(n.obs)){
    U <- runif(n=1)
    Cox.factor <- exp(as.vector(t(beta) %*% c(Trt[i],M.mat[i,-1],Trt_M[i,-1])))
    T[i] <- lambda*(-log(1-U)/Cox.factor)^(1/k)
    C[i] <- runif(n=1,min=a,max=b)
    Y[i] <- min(T[i],C[i])
    Delta[i] <- T[i]<=C[i]
  }
  survData <- data.frame(Subj=(1:sum(n.obs)),Y=Y,Delta=Delta,T=T,C=C,
                         Trt=Trt,M=M,M.mat=M.mat,Trt_M=Trt_M)
  return(survData) 
}


dataWeibull.normalCensor <- function(m,lambda,k,beta,n.obs,mu,sigma,a,b){
  # assuming normal random censoring
  # m: the total number of marker groups
  # lambda: Weibull distribution's scale parameter (for T)
  # k: Weibull distribution's shape parameter (for T)
  # beta: the slope parameter in the Cox PH model, a vector of length (2m-1)
  # n.obs: sample size vector (for each subgroup under each treatment arm), a vector of length 2m
  # mu: mean for the censoring distribution (assume Normal)
  # sigma: SD for the censoring distribution (assume Normal)
  
  # the data set is sorted by Trt (first chunk is control, second chunk is treatment)
  # and then Marker (group 1 to m)
  Trt <- c(rep(0,sum(n.obs[1:m])),rep(1,sum(n.obs[(m+1):(2*m)])))
  M <- rep(0,sum(n.obs))  # M is the vector with values of 1 to m indicating the marker groups
  M.mat <- mat.or.vec(sum(n.obs),m) # M.mat is the matrix (jth column = 1 if M=j or 0)
  
  M[1:n.obs[1]] <- 1
  M[(sum(n.obs[1:m])+1) : sum(n.obs[1:(m+1)])] <- 1
  M.mat[1:n.obs[1],1] <- 1
  M.mat[(sum(n.obs[1:m])+1) : sum(n.obs[1:(m+1)]), 1] <- 1
  
  for(i in 2:m){
    M[(sum(n.obs[1:(i-1)])+1) : sum(n.obs[1:i])] <- i
    M[(sum(n.obs[1:(m+i-1)])+1) : sum(n.obs[1:(m+i)])] <- i
    
    M.mat[(sum(n.obs[1:(i-1)])+1) : sum(n.obs[1:i]), i] <- 1
    M.mat[(sum(n.obs[1:(m+i-1)])+1) : sum(n.obs[1:(m+i)]), i] <- 1
  }
  
  Trt_M <- Trt * M.mat  # this is the interaction matrix (dim=sum(n.obs)*m)
  
  T <- rep(0,sum(n.obs))
  C <- rep(0,sum(n.obs))
  Y <- rep(0,sum(n.obs))
  Delta <- rep(0,sum(n.obs))
  
  for(i in 1:sum(n.obs)){
    U <- runif(n=1)
    Cox.factor <- exp(as.vector(t(beta) %*% c(Trt[i],M.mat[i,-1],Trt_M[i,-1])))
    T[i] <- lambda*(-log(1-U)/Cox.factor)^(1/k)
    C[i] <- rnorm(n=1,mean=mu,sd=sigma)
    if (C[i]<a) C[i]=a
    if (C[i]>b) C[i]=b
    Y[i] <- min(T[i],C[i])
    Delta[i] <- T[i]<=C[i]
  }
  survData <- data.frame(Subj=(1:sum(n.obs)),Y=Y,Delta=Delta,T=T,C=C,
                         Trt=Trt,M=M,M.mat=M.mat,Trt_M=Trt_M)
  return(survData) 
}

# fit a weibull regression and get the parameter estimates and their covariance matrix
weib <- function(data, m){
  # m: the total number of marker groups
  # data are generated from "dataWeibull" function
  
  fit.weib <- weibreg(Surv(Y,Delta) ~ Trt+as.factor(M)+Trt*as.factor(M), data=data)
  fit.LogCoeff <- fit.weib$coef
  fit.Coeff <- exp(fit.LogCoeff) # the coefficient in the "original" scale (exponentiated from the log scale)
  # change the order of the parameters (so that they correspond to "lambda","k","theta1",...,"theta_2m-1")
  fit.Coeff <- fit.Coeff[c(2*m, (2*m+1), 1:(2*m-1))]
  names(fit.Coeff)[1:2]<-c("Scale","Shape")
  # change the order of the var parameters accordingly 
  fit.LogVar <- fit.weib$var
  tmp1 <- fit.LogVar[1:(2*m-1), 1:(2*m-1)]
  tmp2 <- fit.LogVar[(2*m):(2*m+1), (2*m):(2*m+1)]
  tmp3 <- fit.LogVar[1:(2*m-1),(2*m):(2*m+1)]
  
  fit.LogVar.2 <- rbind(cbind(tmp2,t(tmp3)),cbind(tmp3,tmp1))
  # delta method to compute the Variance matrix for the parameter estimates in the "original" scale
  fit.Var <- diag(fit.Coeff) %*% fit.LogVar.2 %*% diag(fit.Coeff) 
  
  return(list(log.coef=fit.LogCoeff,log.var=fit.LogVar,coef=fit.Coeff,var=fit.Var))
}


# compute the survival function for each subgroup + the combined groups
surv.fun <- function(lambda,k,theta,p.neg,time){
  
  # lambda: Weibull distribution's scale parameter 
  # k: Weibull distribution's shape parameter 
  # theta=(theta1,theta2,theta3): =exp(beta) where beta is the slope parameter in the Cox PH model
  # p.neg: Prob(g-) (assume Prob(g-|C)=Prob(g-|Rx))
  # time: the time point to evaluate the survival functions, can be a scalar or vector
  
  S.C.neg <- exp(-(time/lambda)^k)
  S.C.pos <- S.C.neg^theta[2]
  S.Rx.neg <- S.C.neg^theta[1]
  S.Rx.pos <- S.C.neg^(theta[1]*theta[2]*theta[3])
  
  S.C <- S.C.neg*p.neg + S.C.pos*(1-p.neg)
  S.Rx <- S.Rx.neg*p.neg + S.Rx.pos*(1-p.neg)
  
  return(list(S.C.neg=S.C.neg,S.C.pos=S.C.pos,S.Rx.neg=S.Rx.neg,S.Rx.pos=S.Rx.pos,S.C=S.C,S.Rx=S.Rx))
  
}
