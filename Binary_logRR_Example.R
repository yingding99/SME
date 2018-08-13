

data.1 <- read.csv("PDL1_NEJM_1pct.csv",head=T)
source("Binary outcome_logRR.R")

# fit logistic model
fit.logis <- logistic(data=data.1, m=2)

logRR.all3.logis <- logistic.log.RR.dt.dim3(theta=fit.logis$theta, 
                                            p.neg=1-sum(data.1$M)/dim(data.1)[1], 
                                            param.var=fit.logis$fit.theta.Var)

RR.CI3.logis <- log.RR.simu.CI(seed=304, alpha=0.05, 
                               log.RR=c(logRR.all3.logis$log.RR.neg,logRR.all3.logis$log.RR.pos,logRR.all3.logis$log.RR.mix), 
                               var.log.RR=logRR.all3.logis$var.log.RR)

# fit loglinear model
fit.loglin <- loglinear(data=data.1, m=2)

logRR.all3.loglin <- loglinear.log.RR.dt.dim3(theta=fit.loglin$theta, 
                                              p.neg=1-sum(data.10$M)/dim(data.10)[1], 
                                              param.var=fit.loglin$fit.theta.Var)

RR.CI3.loglin <- log.RR.simu.CI(seed=304, alpha=0.05, 
                                log.RR=c(logRR.all3.loglin$log.RR.neg,logRR.all3.loglin$log.RR.pos,logRR.all3.loglin$log.RR.mix), 
                                var.log.RR=logRR.all3.loglin$var.log.RR)

