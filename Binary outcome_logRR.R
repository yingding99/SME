#-------------------------------------------------------------------------#
# Binary outcome_logRR.R includes 5 functions: 
#	1. logistic: fit a logistic regression model
# 2. logistic.log.RR.dt.dim3: delta method for log(RR) from logistic model
# 3. loglinear: fit a loglinear regression model
#	4. loglinear.log.RR.dt.dim3:: delta method for log(RR) from loglinear model
#	5. log.RR.simu.CI: simultaneous CIs computation for log(RR)
#-------------------------------------------------------------------------#

require(msm)
require(multcomp)

logistic <- function(data, m){
  # m: the total number of marker groups
  # data is a data frame with the following format:
  # Y: binary outcome (0/1)
  # Trt: treatment indicator (0=C, 1=Rx)
  # M: Marker (for example: 0=g-, 1=g+)
  
  fit.logit <- glm(Y ~ Trt+as.factor(M)+Trt*as.factor(M), family=binomial(logit), data=data)
  beta <- fit.logit$coef	# beta1, ... , beta4
  theta <- exp(beta)	# theta1, ... ,theta4
  
  fit.beta.Var <- vcov(fit.logit)

  # apply delta method : using the function "deltamethod" from {msm} package
  fit.theta.Var <- deltamethod ( list( ~exp(x1), ~exp(x2), ~exp(x3),~exp(x4)), beta, fit.beta.Var, ses=FALSE)		  
  
  return(list(beta=beta, fit.beta.Var=fit.beta.Var, theta=theta, fit.theta.Var=fit.theta.Var))
}


# Use Relative Response (RR) as the efficacy measure

logistic.log.RR.dt.dim3 <- function(theta, p.neg, param.var){
	# theta = exp(beta), beta are regression coefficient estimates from the logistic model
  # p.neg = population proportion of g-
  # param.var = estimated variance covariance for theta 
  
	# calculate response rate for each trt-by-marker subgroup
	response.C.neg = theta[1] / (1 + theta[1])
	response.C.pos = theta[1] * theta[3] / (1 + theta[1] * theta[3])
	response.Rx.neg = theta[1] * theta[2] / (1 + theta[1] * theta[2])
	response.Rx.pos = theta[1] * theta[2] * theta[3] * theta[4] / (1 + theta[1] * theta[2] * theta[3] * theta[4])
	
	# calculate response rate for each trt group (with marker subgroups combined)
	response.C = p.neg * (theta[1] / (1 + theta[1])) + (1-p.neg) * (theta[1] * theta[3] / (1 + theta[1] * theta[3]))
	response.Rx = p.neg * (theta[1] * theta[2] / (1 + theta[1] * theta[2])) + (1-p.neg) * (theta[1] * theta[2] * theta[3] * theta[4] / (1 + theta[1] * theta[2] * theta[3] * theta[4]))
	
	# calculate relative response (RR) for each marker group and the combined mixture group
	RR.neg <- response.Rx.neg/response.C.neg
	RR.pos <- response.Rx.pos/response.C.pos
	RR.mix <- response.Rx/response.C
	
	# log relative response
	log.RR.neg <- log(response.Rx.neg/response.C.neg)
	log.RR.pos <- log(response.Rx.pos/response.C.pos)
	log.RR.mix <- log(response.Rx/response.C)
	
	# compute the variance estimates for the relative response

	## extra variable p.neg build up the formula as a string, and convert to a formula
	form <- sprintf("~ log( ( %f*x2 / (1+x1*x2) + (1- %f)*x2*x3*x4 / (1+x1*x2*x3*x4) ) / ( %f / (1+x1) + (1- %f)*x3/ (1+x1*x3) ) )", p.neg, p.neg, p.neg, p.neg)
	var.log.RR <- deltamethod ( list( ~log( (x2*(1+x1)) / (1+x1*x2) ), ~log((x2*x4*(1+x1*x3)) / (1+x1*x2*x3*x4)), as.formula(form)), theta, param.var, ses=FALSE)
	
	return(list(response.C.neg=response.C.neg, 		response.C.pos=response.C.pos,
					response.Rx.neg=response.Rx.neg,	response.Rx.pos=response.Rx.pos,
					response.C=response.C, response.Rx=response.Rx,
					RR.neg=RR.neg, RR.pos=RR.pos, RR.mix=RR.mix,
					log.RR.neg=log.RR.neg, log.RR.pos=log.RR.pos, log.RR.mix=log.RR.mix, var.log.RR=var.log.RR))
}


loglinear <- function(data, m){
  # m: the total number of marker groups
  # data is a data frame with the following format:
  # Y: binary outcome (0/1)
  # Trt: treatment indicator (0=C, 1=Rx)
  # M: Marker (for example: 0=g-, 1=g+)
  
  fit.loglinear <- glm(Y ~ Trt+as.factor(M)+Trt*as.factor(M), family=poisson, data=data)
  beta <- fit.loglinear$coef	# beta1, ... , beta4
  theta <- exp(beta)	# theta1, ... ,theta4
  
  fit.beta.Var <- vcov(fit.loglinear)

  # apply delta method : using the function "deltamethod" from {msm} package
  fit.theta.Var <- deltamethod ( list( ~exp(x1), ~exp(x2), ~exp(x3),~exp(x4)), beta, fit.beta.Var, ses=FALSE)
  
  return(list(beta=beta, fit.beta.Var=fit.beta.Var, theta=theta, fit.theta.Var=fit.theta.Var))
}


loglinear.log.RR.dt.dim3 <- function(theta, p.neg, param.var){
  # theta = exp(beta), beta are regression coefficient estimates from the loglinear model
  # p.neg = population proportion of g-
  # param.var = estimated variance covariance for theta 
  
	# calculate response rate for each trt-by-marker subgroup
	response.C.neg = theta[1]
	response.C.pos = theta[1]*theta[3]
	response.Rx.neg = theta[1]*theta[2]
	response.Rx.pos = theta[1]*theta[2]*theta[3]*theta[4]
	
	# calculate response rate for each trt group (with marker subgroups combined)
	response.C = p.neg*theta[1] + (1-p.neg)*theta[1]*theta[3]
	response.Rx = p.neg*theta[1]*theta[2] + (1-p.neg)*theta[1]*theta[2]*theta[3]*theta[4]
	
	# calculate relative response (RR) for each marker group and the combined mixture group
	RR.neg <- response.Rx.neg/response.C.neg
	RR.pos <- response.Rx.pos/response.C.pos
	RR.mix <- response.Rx/response.C
	
	# log relative response
	log.RR.neg <- log(response.Rx.neg/response.C.neg)
	log.RR.pos <- log(response.Rx.pos/response.C.pos)
	log.RR.mix <- log(response.Rx/response.C)
	
	# compute the variance estimates for the relative response

	## extra variable p.neg build up the formula as a string, and convert to a formula.
	form <- sprintf("~ log(( %f * x2 + (1- %f) * x2 * x3 * x4)/( %f + (1- %f) * x3))", p.neg, p.neg, p.neg, p.neg)
	var.log.RR <- deltamethod ( list( ~log(x2), ~log(x2 * x4), as.formula(form)), theta, param.var, ses=FALSE)
	
	return(list(response.C.neg=response.C.neg, 		response.C.pos=response.C.pos,
				response.Rx.neg=response.Rx.neg,	response.Rx.pos=response.Rx.pos,
				response.C=response.C, response.Rx=response.Rx,
				RR.neg=RR.neg, RR.pos=RR.pos, RR.mix=RR.mix,
				log.RR.neg=log.RR.neg, log.RR.pos=log.RR.pos, log.RR.mix=log.RR.mix, var.log.RR=var.log.RR))
}


# Simulatenous CIs computation using "qmvnorm" from the {multcomp} package
log.RR.simu.CI <- function(seed=12345, alpha, log.RR, var.log.RR){
  # qmvnorm is simulation based, provide a seed (for reproducibility)
  # input includes: 
  # alpha: desired tail level
  # log.RR: log(Relative Response) for each marker subgroup and the overall population 
  
  set.seed(seed)
  log.CIs <- qmvnorm(1-alpha, corr=cov2cor(var.log.RR), tail = "both")$quantile*sqrt(diag(var.log.RR))

  RR.CI.lbd <- exp(log.RR - log.CIs)
  RR.CI.ubd <- exp(log.RR + log.CIs)
  
  return(list(RR.CI.lbd=RR.CI.lbd, RR.CI.ubd=RR.CI.ubd))
}

