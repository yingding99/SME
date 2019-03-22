setwd("C:\\Users\\Ying Ding\\Box Sync\\SubgroupBookChapter\\Data_Rcode\\example_1.3.4")


source("Weibull_Fitting.R") # need to install packages survival, eha, rootSolve
source("Weibull_Median_2grp.R")

library(mvtnorm) # mvtnorm needs a R version 3.5.0

# read in data
JCO.pos <- read.csv("JCO2013_IHC_METpositive.csv")
JCO.pos <- as.data.frame(apply(JCO.pos,2,as.numeric))


p.neg <- sum(JCO.pos$IHC==2)/nrow(JCO.pos) 
JCO.pos$M[JCO.pos$IHC==2] <- 0 # g- group
JCO.pos$M[JCO.pos$IHC==3] <- 1 # g+ group
JCO.pos$Trt_M <- JCO.pos$Trt*JCO.pos$M

#------------------------#
# Fit weibull regression #
#------------------------#
fit <- weib(JCO.pos,2)
# log.param are the log scaled parameter estimates from weibull fitting
log.param.est <- fit$log.coef
log.param.var <- diag(fit$log.var)
# param are the original scale (exponential of log.param) parameter estimates
param.est <- fit$coef
param.var <- fit$var
param.var.diag <- diag(param.var)

#-------------------------------------------------------------------------------#
# compute median survival times and their associated variance covariance matrix #
#-------------------------------------------------------------------------------#
comp.median <- median.surv(lambda=param.est[1],k=param.est[2],theta=param.est[3:5],p.neg=p.neg)
median.surv.est <- unlist(comp.median)
median.surv.est###this returns the median survival time for Trt and Ctrl group in 2+, 3+ and {2+,3+}

comp.median.var.dt<-var.median.surv.dt.dim6(med.C.neg=median.surv.est[1], med.C.pos=median.surv.est[2], 
                        med.Rx.neg=median.surv.est[3], med.Rx.pos=median.surv.est[4],
                        med.C=median.surv.est[5],med.Rx=median.surv.est[6],p.neg=p.neg,
                        lambda=param.est[1],k=param.est[2],theta=param.est[3:5],
                        param.var=fit$var)$median.var
median.surv.var.dt<-diag(comp.median.var.dt)

#------------------------------------------------------------------------------------#
# calculate the CIs for the log ratio of median survival times then take exponantial #
#------------------------------------------------------------------------------------#

median.ratio <- log.ratio.median.surv.dt.dim3(lambda=param.est[1],k=param.est[2],theta=param.est[3:5],p.neg,param.var)
exp(unlist(median.ratio)[1:3]) #this returns the ratio of median survival in 2+, 3+ and {2+,3+}


# the following returns the 3x3 variance covariance matrix for the ratios of medians (g-, g+, and combined)
sigma.mat.ratio <- log.ratio.median.surv.dt.dim3(lambda=param.est[1],k=param.est[2],theta=param.est[3:5],p.neg,param.var)$var.ratio


#the following line returns the lower bound of simultaneous 95% CI
round(exp(unlist(median.ratio)[1:3]-qmvnorm(0.95, corr=cov2cor(sigma.mat.ratio), tail = "both")$quantile*sqrt(diag(sigma.mat.ratio))),2)
#the following line returns the higher bound of simultaneous 95% CI
round(exp(unlist(median.ratio)[1:3]+qmvnorm(0.95, corr=cov2cor(sigma.mat.ratio), tail = "both")$quantile*sqrt(diag(sigma.mat.ratio))),2)  			



