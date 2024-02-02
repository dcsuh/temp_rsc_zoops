setwd("C:/Users/straussa/Documents/Research/Hall Lab/f spores/data")
library(bbmle) # optimizer wrapper
library(deSolve) # integrates ODEs

dataset <- read.csv("Traits1-3.csv", na.strings=c("NA", "NULL"))
dataset$index <- paste(dataset$round,dataset$clone) # create unique id for each clone in each round
dataset <- dataset[dataset$clone!="Down.282" , ] # round 3 has 0 infections; round 2 really low sample size
dataset <- dataset[dataset$clone!="Cerio" , ] # different species

dataset[dataset$death.age<11 & !is.na(dataset$death.age) , ]$infected <- NA # list the early deaths as NA
dataset[dataset$death.age<11 & !is.na(dataset$death.age) , ]$uninfected <- NA
length(dataset[dataset$death.age<11 & !is.na(dataset$death.age) , ]$infected)

clones <- unique(dataset$index) # how many unique?
vol <- 15 # in ml
f.min <- 3.6 # this is the bottom 2.5% 

options(warn=1) #now warnings print as they happen; set warn=0 if want at end


######################################################################################################
######################################################################################################
#################################                          ###########################################
#################################      GENOTYPE LOOPS      ###########################################
#################################                          ###########################################
######################################################################################################
######################################################################################################

# these functions loop through each genotype and fit 2-4 parameters for each:
# constant or exponential foraging crossed by constant or exponential u 


######################################################################################################
##################     CONSTANT FORAGING/EXPOSURE; CONSTANT SUSCEPTIBILITY     #######################
######################################################################################################

con.f.con.u.nll <- function(f.0.hat, u){ # the likelihood function - will fit f and u
  simulator <- function(times,y,params){  # simulate foraging and infection dynamics
    Z = y[1]; A = y[2]; S = y[3]; I = y[4]
    with(as.list(params),{
      dZ = -(S+I)*Z*(f.0.hat*(length^2)) # constant foraging
      dA = -(S+I)*A*(f.0.hat*(length^2)) # constant foraging
      dS = -u*S*Z*(f.0.hat*(length^2)) # constant foraging and u
      dI = u*S*Z*(f.0.hat*(length^2)) # constant foraging and u
      res = c(dZ, dA, dS, dI)
      list(res)})}
  transm <- function(Z, A, length, time, f.0.hat, u){ # function to vectorize the simulation
    sim <- lsoda(y=c(Z=Z, A=A, S=1, I=0), func=simulator, rtol=1e-4, atol=1e-4, # larger is faster but less accurate
              parms=c(f.0.hat=f.0.hat/vol, u=u/10000, length=length), times=seq(from=0, to=time, by=time))}
  transm.sum <- data.frame(t(mapply(transm, Z=g.dat$exposure*vol, A=g.dat$control, length=g.dat$mm, time=g.dat$time,
                                    f.0.hat=f.0.hat, u=u)))
  colnames(transm.sum) <- c("X1", "time", "Z.start", "Z.end", "A.start", "A.end", "S.start", "S.end", "I.start", "I.end")
  transm.sum$f.resid <- log(transm.sum$A.end) - log (g.dat$fluor.) # residuals of feeding rate
  transm.sum$u.LL <- dbinom(x=g.dat$infected, size=1, prob=transm.sum$I.end, log=TRUE) # log likelihood of infection data
  transm.sum$f.LL <- dnorm(x=transm.sum$f.resid, mean=0, sd=sd(transm.sum$f.resid, na.rm=TRUE), log=TRUE) # ll of feeding data
  # throw out nonsensical parameters like negative u or negative infection prevalence
  ifelse(u < 0 | any(transm.sum$I.end<0), 10000, # return a large number to prevent nonsense fit
         -sum(transm.sum$u.LL, na.rm=T) + -sum(transm.sum$f.LL, na.rm=T)) # otherwise sum log liklihoods and send to optimizer
}

summary.con.f.con.u <- data.frame(index = clones) # summary spreadsheet to fill out
for (k in 1:length(clones)){ # loop through genotypes
  g.dat <- dataset[dataset$index==clones[k] , ] # subset by genotype-round combination
  print(g.dat$index[1]) # to identify where warnings come from
  con.f.con.u.fit = mle2(minuslogl=con.f.con.u.nll, method="Nelder-Mead", skip.hessian=T, 
                         start=list(f.0.hat=8, u=4.5), 
                         control=list(parscale=c(f.0.hat=8, u=4.5), maxit=10000)) 
  summary.con.f.con.u$func[k]="con.f.con.u"
  summary.con.f.con.u$f.0.hat[k]=as.numeric(coef(con.f.con.u.fit)[1])
  summary.con.f.con.u$alpha[k]=NA
  summary.con.f.con.u$u[k]=as.numeric(coef(con.f.con.u.fit)[2])/10000 # rescaling u back
  summary.con.f.con.u$w[k]=NA
  summary.con.f.con.u$length[k]=mean(g.dat$mm)
  summary.con.f.con.u$time[k]=mean(g.dat$time)
  summary.con.f.con.u$loglik[k]=summary(con.f.con.u.fit)@m2logL/-2
  } # weird notation to extract log likelihood
 
write.csv(summary.con.f.con.u, file = "summary.con.f.con.u.csv")
# with rtol & atol = 1e-4, all 27 clones take only 4 minutes. 


######################################################################################################
################     EXPONENTIAL FORAGING/EXPOSURE; CONSTANT SUSCEPTIBILITY     ######################
######################################################################################################

exp.f.con.u.nll <- function(f.0.hat, alpha, u){ # will fit f, alpha, and u
  simulator <- function(times,y,params){
    Z = y[1]; A = y[2]; S = y[3]; I = y[4]
    with(as.list(params),{
      dZ = -(S+I)*Z*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol) # exponential f; con u
      dA = -(S+I)*A*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol) # exponential f; con u
      dS = -u*S*Z*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol) # exponential f; con u
      dI = u*S*Z*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol) # exponential f; con u
      res = c(dZ, dA, dS, dI)
      list(res)})}
  transm <- function(Z, A, length, time, f.0.hat, u, alpha){ # function to vectorize the simulation
    sim <- lsoda(y=c(Z=Z, A=A, S=1, I=0), func=simulator, rtol=1e-4, atol=1e-4, # larger is faster but less accurate
              parms=c(f.0.hat=f.0.hat/vol, u=u/10000, alpha=alpha/(vol*1000),
                      # Z in sim in sp/tube; want alpha estimated as ml per spore
                      length=length), times=seq(from=0, to=time, by=time))}
  transm.sum <- data.frame(t(mapply(transm, Z=g.dat$exposure*vol, A=g.dat$control, length=g.dat$mm, time=g.dat$time,
                                    f.0.hat=f.0.hat, u=u, alpha=alpha)))
  colnames(transm.sum) <- c("X1", "time", "Z.start", "Z.end", "A.start", "A.end", "S.start", "S.end", "I.start", "I.end")
  transm.sum$f.resid <- log(transm.sum$A.end) - log (g.dat$fluor.) # residuals of feeding rate
  transm.sum$u.LL <- dbinom(x=g.dat$infected, size=1, prob=transm.sum$I.end, log=TRUE) # log likelihood of infection data
  transm.sum$f.LL <- dnorm(x=transm.sum$f.resid, mean=0, sd=sd(transm.sum$f.resid, na.rm=TRUE), log=TRUE) # ll of feeding data
  # throw out nonsensical parameters like negative u or negative infection prevalence
  ifelse(u < 0 | any(transm.sum$I.end<0), 10000, # return a large number to prevent nonsense fit
         -sum(transm.sum$u.LL, na.rm=T) + -sum(transm.sum$f.LL, na.rm=T)) # otherwise sum log liklihoods and send to optimizer
  }

summary.exp.f.con.u=data.frame(index=clones)
  for (k in 1:length(clones)){
    g.dat <- dataset[dataset$index==clones[k] , ]
    print(g.dat$index[1])
    exp.f.con.u.fit = mle2(minuslogl=exp.f.con.u.nll, method="Nelder-Mead", skip.hessian=T,
                         start=list(f.0.hat=11, alpha=-2, u=5),
                         control=list(parscale=c(f.0.hat=11, alpha=2, u=5), maxit=10000))
    summary.exp.f.con.u$func[k]="exp.f.con.u"
    summary.exp.f.con.u$f.0.hat[k]=as.numeric(coef(exp.f.con.u.fit)[1])
    summary.exp.f.con.u$alpha[k]=as.numeric(coef(exp.f.con.u.fit)[2])/1000
    summary.exp.f.con.u$u[k]=as.numeric(coef(exp.f.con.u.fit)[3])/10000
    summary.exp.f.con.u$w[k]=NA
    summary.exp.f.con.u$length[k]=mean(g.dat$mm)
    summary.exp.f.con.u$time[k]=mean(g.dat$time)
    summary.exp.f.con.u$loglik[k]=summary(exp.f.con.u.fit)@m2logL/-2
  }

# WITH rtol & atol = 1e-4: all 27 clones take 13 minutes
summary.exp.f.con.u
write.csv(summary.exp.f.con.u, file = "summary.exp.f.con.u.csv")


######################################################################################################
###############     CONSTANT FORAGING/EXPOSURE; EXPONENTIAL SUSCEPTIBILITY     #######################
######################################################################################################

con.f.exp.u.nll <- function(f.0.hat, u, w){ # will fit f, u, and w
  simulator <- function(times,y,params){
    Z = y[1]; A = y[2]; S = y[3]; I = y[4]
    with(as.list(params),{
      dZ = -(S+I)*Z*f.0.hat*(length^2) # specify form: con f; exp u
      dA = -(S+I)*A*f.0.hat*(length^2) # specify form: con f; exp u
      dS = -u*exp(w*Z*f.0.hat*(length^2))*S*Z*f.0.hat*(length^2) # specify form: con f; exp u
      dI = u*exp(w*Z*f.0.hat*(length^2))*S*Z*f.0.hat*(length^2) # specify form: con f; exp u
      res = c(dZ, dA, dS, dI)
      list(res)})}
  transm <- function(Z, A, length, time, f.0.hat, u, w){ # function to vectorize the simulation
    sim <- lsoda(y=c(Z=Z, A=A, S=1, I=0), func=simulator, rtol=1e-4, atol=1e-4, # larger is faster but less accurate
              parms=c(f.0.hat=f.0.hat/vol, u=u/10000, w=w/10000,
                      length=length), times=seq(from=0, to=time, by=time))}
  transm.sum <- data.frame(t(mapply(transm, Z=g.dat$exposure*vol, A=g.dat$control, length=g.dat$mm, time=g.dat$time,
                                    f.0.hat=f.0.hat, u=u, w=w)))
  colnames(transm.sum) <- c("X1", "time", "Z.start", "Z.end", "A.start", "A.end", "S.start", "S.end", "I.start", "I.end")
  transm.sum$f.resid <- log(transm.sum$A.end) - log (g.dat$fluor.) # residuals of feeding rate
  transm.sum$u.LL <- dbinom(x=g.dat$infected, size=1, prob=transm.sum$I.end, log=TRUE) # log likelihood of infection data
  transm.sum$f.LL <- dnorm(x=transm.sum$f.resid, mean=0, sd=sd(transm.sum$f.resid, na.rm=TRUE), log=TRUE) # ll of feeding data
  # throw out nonsensical parameters like negative u or negative infection prevalence
  ifelse(u < 0 | any(transm.sum$I.end<0), 10000, # return a large number to prevent nonsense fit
         -sum(transm.sum$u.LL, na.rm=T) + -sum(transm.sum$f.LL, na.rm=T)) # otherwise sum log liklihoods and send to optimizer
}

summary.con.f.exp.u=data.frame(index=clones)
for (k in 1:length(clones)){
  g.dat=dataset[dataset$index==clones[k] , ]
  print(g.dat$index[1])
  con.f.exp.u.fit = mle2(minuslogl=con.f.exp.u.nll, skip.hessian=T, method="Nelder-Mead",
                         start=list(f.0.hat=8, u=15, w=-5), # start u and w small to avoid warnings
                         control=list(parscale=c(f.0.hat=8, u=15, w=5), maxit=10000))
  summary.con.f.exp.u$func[k]="con.f.exp.u"
  summary.con.f.exp.u$f.0.hat[k]=as.numeric(coef(con.f.exp.u.fit)[1])
  summary.con.f.exp.u$alpha[k]=NA
  summary.con.f.exp.u$u[k]=as.numeric(coef(con.f.exp.u.fit)[2])/10000 # rescaling back to u
  summary.con.f.exp.u$w[k]=as.numeric(coef(con.f.exp.u.fit)[3])/10000 # rescaling back to ml per spore
  summary.con.f.exp.u$length[k]=mean(g.dat$mm)
  summary.con.f.exp.u$time[k]=mean(g.dat$time)
  summary.con.f.exp.u$loglik[k]=summary(con.f.exp.u.fit)@m2logL/-2
}

write.csv(summary.con.f.exp.u, file = "summary.con.f.exp.u.csv")
# all 27 clones take 10 minutes 


######################################################################################################
################     EXPONENTIAL FORAGING/EXPOSURE; EXPONENTIAL SUSCEPTIBILITY     ###################
######################################################################################################

exp.f.exp.u.nll <- function(f.0.hat, alpha, u, w){ 
  simulator <- function(times,y,params){  
    Z = y[1]; A = y[2]; S = y[3]; I = y[4]
    with(as.list(params),{
      dZ = -(S+I)*Z*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol) 
      dA = -(S+I)*A*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol) 
      dS = -u*exp(w*Z*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol))*
        S*Z*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol) 
      dI = u*exp(w*Z*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol))*
        S*Z*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol) 
      res = c(dZ, dA, dS, dI)
      list(res)})}
  transm <- function(Z, A, length, time, f.0.hat, alpha, u, w){ # function to vectorize the simulation
    sim <- lsoda(y=c(Z=Z, A=A, S=1, I=0), func=simulator, rtol=1e-4, atol=1e-4, # larger is faster but less accurate
              parms=c(f.0.hat=f.0.hat/vol, alpha=alpha/(vol*1000), u=u/10000, w=w/10000, 
                      length=length), times=seq(from=0, to=time, by=time))} 
  transm.sum <- data.frame(t(mapply(transm, Z=g.dat$exposure*vol, A=g.dat$control, length=g.dat$mm, time=g.dat$time,
                                    f.0.hat=f.0.hat, alpha=alpha, u=u, w=w)))
  colnames(transm.sum) <- c("X1", "time", "Z.start", "Z.end", "A.start", "A.end", "S.start", "S.end", "I.start", "I.end")
  transm.sum$f.resid <- log(transm.sum$A.end) - log (g.dat$fluor.) # residuals of feeding rate
  transm.sum$u.LL <- dbinom(x=g.dat$infected, size=1, prob=transm.sum$I.end, log=TRUE) # log likelihood of infection data
  transm.sum$f.LL <- dnorm(x=transm.sum$f.resid, mean=0, sd=sd(transm.sum$f.resid, na.rm=TRUE), log=TRUE) # ll of feeding data
  # throw out nonsensical parameters like negative u or negative infection prevalence
  ifelse(u < 0 | any(transm.sum$I.end<0), 10000, # return a large number to prevent nonsense fit
         -sum(transm.sum$u.LL, na.rm=T) + -sum(transm.sum$f.LL, na.rm=T)) # otherwise sum log liklihoods and send to optimizer
  }

summary.exp.f.exp.u <- data.frame(index=clones) 
  for (k in 1:length(clones)){ 
    g.dat <- dataset[dataset$index==clones[k] , ]
    print(g.dat$index[1])
    exp.f.exp.u.fit = mle2(minuslogl=exp.f.exp.u.nll, method = "Nelder-Mead", skip.hessian=T,
                           start=list(f.0.hat=12, alpha=-1.2, u=11, w=-3),
                           control=list(parscale=c(f.0.hat=12, alpha=1.2, u=11, w=3), maxit=10000))
    summary.exp.f.exp.u$func[k]="exp.f.exp.u"
    summary.exp.f.exp.u$f.0.hat[k]=as.numeric(coef(exp.f.exp.u.fit)[1])
    summary.exp.f.exp.u$alpha[k]=as.numeric(coef(exp.f.exp.u.fit)[2])/1000 # rescaling back to ml per spore
    summary.exp.f.exp.u$u[k]=as.numeric(coef(exp.f.exp.u.fit)[3])/10000 # rescaling back to u
    summary.exp.f.exp.u$w[k]=as.numeric(coef(exp.f.exp.u.fit)[4])/10000 # rescaling back to ml per spore
    summary.exp.f.exp.u$length[k]=mean(g.dat$mm)
    summary.exp.f.exp.u$time[k]=mean(g.dat$time)
    summary.exp.f.exp.u$loglik[k]=summary(exp.f.exp.u.fit)@m2logL/-2
    } 

summary.exp.f.exp.u
write.csv(summary.exp.f.exp.u, file = "summary.exp.f.exp.u.2.csv")


######################################################################################################
###################     EXPONENTIAL FORAGING/EXPOSURE; LINEAR SUSCEPTIBILITY     #####################
######################################################################################################

exp.f.lin.u.nll <- function(f.0.hat, alpha, u, w){ 
  simulator <- function(times,y,params){  
    Z = y[1]; A = y[2]; S = y[3]; I = y[4]
    with(as.list(params),{
      dZ = -(S+I)*Z*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol) # specify form: exponential alpha
      dA = -(S+I)*A*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol)
      dS = -u*(1+(w*Z*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol)))*
        S*Z*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol) # form: linear u
      dI = u*(1+(w*Z*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol)))*
        S*Z*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol) 
      res = c(dZ, dA, dS, dI)
      list(res)})}
  transm <- function(Z, A, length, time, f.0.hat, alpha, u, w){ # function to vectorize the simulation
    sim <- lsoda(y=c(Z=Z, A=A, S=1, I=0), func=simulator, rtol=1e-4, atol=1e-4, # larger is faster but less accurate
              parms=c(f.0.hat=f.0.hat/vol, alpha=alpha/(vol*1000), u=u/10000, w=w/10000, 
                      length=length), times=seq(from=0, to=time, by=time))} 
  transm.sum <- data.frame(t(mapply(transm, Z=g.dat$exposure*vol, A=g.dat$control, length=g.dat$mm, time=g.dat$time,
                                    f.0.hat=f.0.hat, alpha=alpha, u=u, w=w)))
  colnames(transm.sum) <- c("X1", "time", "Z.start", "Z.end", "A.start", "A.end", "S.start", "S.end", "I.start", "I.end")
  transm.sum$f.resid <- log(transm.sum$A.end) - log (g.dat$fluor.) # residuals of feeding rate
  transm.sum$u.LL <- dbinom(x=g.dat$infected, size=1, prob=transm.sum$I.end, log=TRUE) # log likelihood of infection data
  transm.sum$f.LL <- dnorm(x=transm.sum$f.resid, mean=0, sd=sd(transm.sum$f.resid, na.rm=TRUE), log=TRUE) # ll of feeding data
  # throw out nonsensical parameters like negative u or negative infection prevalence
  ifelse(u < 0 | any(transm.sum$I.end<0), 10000, # return a large number to prevent nonsense fit
         -sum(transm.sum$u.LL, na.rm=T) + -sum(transm.sum$f.LL, na.rm=T)) # otherwise sum log liklihoods and send to optimizer
}

summary.exp.f.lin.u <- data.frame(index=clones) 
for (k in 1:length(clones)){ 
  g.dat <- dataset[dataset$index==clones[k] , ]
  print(g.dat$index[1])
  exp.f.lin.u.fit = mle2(minuslogl=exp.f.lin.u.nll, skip.hessian=T, method="Nelder-Mead", 
                         start=list(f.0.hat=10, alpha=-2, u=9, w=-2),
                         control=list(parscale=c(f.0.hat=10, alpha=2, u=9, w=2), maxit=10000)) 
  summary.exp.f.lin.u$func[k]="exp.f.lin.u"
  summary.exp.f.lin.u$f.0.hat[k]=as.numeric(coef(exp.f.lin.u.fit)[1])
  summary.exp.f.lin.u$alpha[k]=as.numeric(coef(exp.f.lin.u.fit)[2])/1000 # rescaling back to ml per spore
  summary.exp.f.lin.u$u[k]=as.numeric(coef(exp.f.lin.u.fit)[3])/10000 # rescaling back to u
  summary.exp.f.lin.u$w[k]=as.numeric(coef(exp.f.lin.u.fit)[4])/10000 # rescaling back to ml per spore
  summary.exp.f.lin.u$length[k]=mean(g.dat$mm)
  summary.exp.f.lin.u$time[k]=mean(g.dat$time)
  summary.exp.f.lin.u$loglik[k]=summary(exp.f.lin.u.fit)@m2logL/-2
} 

summary.exp.f.lin.u
# all 27 clones take 55 minutes. 
write.csv(summary.exp.f.lin.u, file = "summary.exp.f.lin.u.csv")


######################################################################################################
##################     LINEAR FORAGING/EXPOSURE; CONSTANT SUSCEPTIBILITY     #########################
######################################################################################################

lin.f.con.u.nll <- function(f.0.hat, alpha, u){ # will fit f, alpha, and u 
  simulator <- function(times,y,params){  
    Z = y[1]; A = y[2]; S = y[3]; I = y[4]
    with(as.list(params),{
      dZ = -(S+I)*Z*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol) # specify form: linear alpha 
      dA = -(S+I)*A*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol) 
      dS = -u*S*Z*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol) 
      dI = u*S*Z*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol)
      res = c(dZ, dA, dS, dI)
      list(res)})}
  transm <- function(Z, A, length, time, f.0.hat, u, alpha){ # function to vectorize the simulation
    sim <- lsoda(y=c(Z=Z, A=A, S=1, I=0), func=simulator, rtol=1e-4, atol=1e-4, # larger is faster but less accurate
              parms=c(f.0.hat=f.0.hat/vol, u=u/10000, alpha=alpha/(vol*1000),   
                      # Z in sim in sp/tube; want alpha estimated as ml per spore 
                      length=length), times=seq(from=0, to=time, by=time))} 
  transm.sum <- data.frame(t(mapply(transm, Z=g.dat$exposure*vol, A=g.dat$control, length=g.dat$mm, time=g.dat$time,
                                    f.0.hat=f.0.hat, u=u, alpha=alpha)))
  colnames(transm.sum) <- c("X1", "time", "Z.start", "Z.end", "A.start", "A.end", "S.start", "S.end", "I.start", "I.end")
  transm.sum$f.resid <- log(transm.sum$A.end) - log (g.dat$fluor.) # residuals of feeding rate
  transm.sum$u.LL <- dbinom(x=g.dat$infected, size=1, prob=transm.sum$I.end, log=TRUE) # log likelihood of infection data
  transm.sum$f.LL <- dnorm(x=transm.sum$f.resid, mean=0, sd=sd(transm.sum$f.resid, na.rm=TRUE), log=TRUE) # ll of feeding data
  # throw out nonsensical parameters like negative u or negative infection prevalence
  ifelse(u < 0 | any(transm.sum$I.end<0), 10000, # return a large number to prevent nonsense fit
         -sum(transm.sum$u.LL, na.rm=T) + -sum(transm.sum$f.LL, na.rm=T)) # otherwise sum log liklihoods and send to optimizer
}

summary.lin.f.con.u=data.frame(index=clones) 
for (k in 1:length(clones)){ 
  g.dat <- dataset[dataset$index==clones[k] , ] 
  print(g.dat$index[1])
  lin.f.con.u.fit = mle2(minuslogl=lin.f.con.u.nll, skip.hessian=T, method="Nelder-Mead", 
                         start=list(f.0.hat=11, alpha=-0.1, u=5),
                         control=list(parscale=c(f.0.hat=11, alpha=0.1, u=5), maxit=10000))
  summary.lin.f.con.u$func[k]="lin.f.con.u"
  summary.lin.f.con.u$f.0.hat[k]=as.numeric(coef(lin.f.con.u.fit)[1])
  summary.lin.f.con.u$alpha[k]=as.numeric(coef(lin.f.con.u.fit)[2])/1000 
  summary.lin.f.con.u$u[k]=as.numeric(coef(lin.f.con.u.fit)[3])/10000 
  summary.lin.f.con.u$w[k]=NA
  summary.lin.f.con.u$length[k]=mean(g.dat$mm)
  summary.lin.f.con.u$time[k]=mean(g.dat$time)
  summary.lin.f.con.u$loglik[k]=summary(lin.f.con.u.fit)@m2logL/-2
}

# WITH rtol & atol = 1e-4, all 27 clones take < 40 min
write.csv(summary.lin.f.con.u, file = "summary.lin.f.con.u.csv")


######################################################################################################
###################     LINEAR FORAGING/EXPOSURE; EXPONENTIAL SUSCEPTIBILITY     #####################
######################################################################################################

lin.f.exp.u.nll <- function(f.0.hat, alpha, u, w){ 
  simulator <- function(times,y,params){  
    Z = y[1]; A = y[2]; S = y[3]; I = y[4]
    with(as.list(params),{
      dZ = -(S+I)*Z*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol) # specify form: linear foraging
      dA = -(S+I)*A*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol)
      dS = -u*exp(w*Z*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol))*
        S*Z*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol)
      dI = u*exp(w*Z*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol))*
        S*Z*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol) 
      res = c(dZ, dA, dS, dI)
      list(res)})}
  transm <- function(Z, A, length, time, f.0.hat, alpha, u, w){ # function to vectorize the simulation
    sim <- lsoda(y=c(Z=Z, A=A, S=1, I=0), func=simulator, rtol=1e-4, atol=1e-4, # larger is faster but less accurate
              parms=c(f.0.hat=f.0.hat/vol, alpha=alpha/(vol*1000), u=u/10000, w=w/10000, 
                      length=length), times=seq(from=0, to=time, by=time))} 
  transm.sum <- data.frame(t(mapply(transm, Z=g.dat$exposure*vol, A=g.dat$control, length=g.dat$mm, time=g.dat$time,
                                    f.0.hat=f.0.hat, alpha=alpha, u=u, w=w)))
  colnames(transm.sum) <- c("X1", "time", "Z.start", "Z.end", "A.start", "A.end", "S.start", "S.end", "I.start", "I.end")
  transm.sum$f.resid <- log(transm.sum$A.end) - log (g.dat$fluor.) # residuals of feeding rate
  transm.sum$u.LL <- dbinom(x=g.dat$infected, size=1, prob=transm.sum$I.end, log=TRUE) # log likelihood of infection data
  transm.sum$f.LL <- dnorm(x=transm.sum$f.resid, mean=0, sd=sd(transm.sum$f.resid, na.rm=TRUE), log=TRUE) # ll of feeding data
  # throw out nonsensical parameters like negative u or negative infection prevalence
  ifelse(u < 0 | any(transm.sum$I.end<0), 10000, # return a large number to prevent nonsense fit
         -sum(transm.sum$u.LL, na.rm=T) + -sum(transm.sum$f.LL, na.rm=T)) # otherwise sum log liklihoods and send to optimizer
}

summary.lin.f.exp.u <- data.frame(index=clones) 
for (k in 1:length(clones)){ 
  g.dat <- dataset[dataset$index==clones[k] , ]
  print(g.dat$index[1])
  lin.f.exp.u.fit = mle2(minuslogl=lin.f.exp.u.nll, skip.hessian=T, method="Nelder-Mead", 
                         start=list(f.0.hat=11, alpha=-0.8, u=18, w=-4),
                         control=list(parscale=c(f.0.hat=11, alpha=0.8, u=18, w=4), maxit=10000)) 
  summary.lin.f.exp.u$func[k]="lin.f.exp.u"
  summary.lin.f.exp.u$f.0.hat[k]=as.numeric(coef(lin.f.exp.u.fit)[1])
  summary.lin.f.exp.u$alpha[k]=as.numeric(coef(lin.f.exp.u.fit)[2])/1000 # rescaling back to ml per spore
  summary.lin.f.exp.u$u[k]=as.numeric(coef(lin.f.exp.u.fit)[3])/10000 # rescaling back to u
  summary.lin.f.exp.u$w[k]=as.numeric(coef(lin.f.exp.u.fit)[4])/10000 # rescaling back to ml per spore
  summary.lin.f.exp.u$length[k]=mean(g.dat$mm)
  summary.lin.f.exp.u$time[k]=mean(g.dat$time)
  summary.lin.f.exp.u$loglik[k]=summary(lin.f.exp.u.fit)@m2logL/-2
} 


# all 27 clones take: 65 min
write.csv(summary.lin.f.exp.u, file = "summary.lin.f.exp.u.csv")


######################################################################################################
###################     CONSTANT FORAGING/EXPOSURE; LINEAR SUSCEPTIBILITY     ########################
######################################################################################################

con.f.lin.u.nll <- function(f.0.hat, u, w){ # will fit f, u, and w
  simulator <- function(times,y,params){  
    Z = y[1]; A = y[2]; S = y[3]; I = y[4]
    with(as.list(params),{
      dZ = -(S+I)*Z*f.0.hat*(length^2) # specify form: con f
      dA = -(S+I)*A*f.0.hat*(length^2)
      dS = -u*S*Z*f.0.hat*(length^2)*(1+(w*Z*f.0.hat*(length^2))) # specify form: con f; lin u
      dI = u*S*Z*f.0.hat*(length^2)*(1+(w*Z*f.0.hat*(length^2)))   
      res = c(dZ, dA, dS, dI)
      list(res)})}
  transm <- function(Z, A, length, time, f.0.hat, u, w){ # function to vectorize the simulation
    sim <- lsoda(y=c(Z=Z, A=A, S=1, I=0), func=simulator, rtol=1e-4, atol=1e-4, # larger is faster but less accurate
              parms=c(f.0.hat=f.0.hat/vol, u=u/10000, w=w/10000,          
                      length=length), times=seq(from=0, to=time, by=time))} 
  transm.sum <- data.frame(t(mapply(transm, Z=g.dat$exposure*vol, A=g.dat$control, length=g.dat$mm, time=g.dat$time,
                                    f.0.hat=f.0.hat, u=u, w=w)))
  colnames(transm.sum) <- c("X1", "time", "Z.start", "Z.end", "A.start", "A.end", "S.start", "S.end", "I.start", "I.end")
  transm.sum$f.resid <- log(transm.sum$A.end) - log (g.dat$fluor.) # residuals of feeding rate
  transm.sum$u.LL <- dbinom(x=g.dat$infected, size=1, prob=transm.sum$I.end, log=TRUE)
  transm.sum$f.LL <- dnorm(x=transm.sum$f.resid, mean=0, sd=sd(transm.sum$f.resid, na.rm=TRUE), log=TRUE) # ll of feeding data
  # throw out nonsensical parameters like negative u or negative infection prevalence
  ifelse(u < 0 | any(transm.sum$I.end<0), 10000, # return a large number to prevent nonsense fit
         -sum(transm.sum$u.LL, na.rm=T) + -sum(transm.sum$f.LL, na.rm=T)) # otherwise sum log liklihoods and send to optimizer
}  

summary.con.f.lin.u=data.frame(index=clones)
for (k in 1:length(clones)){
  g.dat=dataset[dataset$index==clones[k] , ]
  print(g.dat$index[1])
  con.f.lin.u.fit = mle2(minuslogl=con.f.lin.u.nll, skip.hessian=T, method="Nelder-Mead",
                         start=list(f.0.hat=5, u=9, w=-2),
                         control=list(parscale=c(f.0.hat=5, u=9, w=2), maxit=10000))
  summary.con.f.lin.u$func[k]="con.f.lin.u"
  summary.con.f.lin.u$f.0.hat[k]=as.numeric(coef(con.f.lin.u.fit)[1])
  summary.con.f.lin.u$alpha[k]=NA
  summary.con.f.lin.u$u[k]=as.numeric(coef(con.f.lin.u.fit)[2])/10000 # rescaling back to u
  summary.con.f.lin.u$w[k]=as.numeric(coef(con.f.lin.u.fit)[3])/10000 # rescaling back to ml per spore
  summary.con.f.lin.u$length[k]=mean(g.dat$mm)
  summary.con.f.lin.u$time[k]=mean(g.dat$time)
  summary.con.f.lin.u$loglik[k]=summary(con.f.lin.u.fit)@m2logL/-2
}

write.csv(summary.con.f.lin.u, file = "summary.con.f.lin.u.csv") 


######################################################################################################
###################     LINEAR FORAGING/EXPOSURE; LINEAR SUSCEPTIBILITY     ##########################
######################################################################################################

lin.f.lin.u.nll <- function(f.0.hat, alpha, u, w){ 
  simulator <- function(times,y,params){  
    Z = y[1]; A = y[2]; S = y[3]; I = y[4]
    with(as.list(params),{
      dZ = -(S+I)*Z*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol) # specify form: linear foraging
      dA = -(S+I)*A*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol)
      dS = -u*(1+(w*Z*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol)))*
        S*Z*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol)
      dI = u*(1+(w*Z*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol)))*
        S*Z*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol)
      res = c(dZ, dA, dS, dI)
      list(res)})}
  transm <- function(Z, A, length, time, f.0.hat, alpha, u, w){ # function to vectorize the simulation
    sim <- lsoda(y=c(Z=Z, A=A, S=1, I=0), func=simulator, rtol=1e-4, atol=1e-4, # larger is faster but less accurate
              parms=c(f.0.hat=f.0.hat/vol, alpha=alpha/(vol*1000), u=u/10000, w=w/10000, 
                      length=length), times=seq(from=0, to=time, by=time))} 
  transm.sum <- data.frame(t(mapply(transm, Z=g.dat$exposure*vol, A=g.dat$control, length=g.dat$mm, time=g.dat$time,
                                    f.0.hat=f.0.hat, alpha=alpha, u=u, w=w)))
  colnames(transm.sum) <- c("X1", "time", "Z.start", "Z.end", "A.start", "A.end", "S.start", "S.end", "I.start", "I.end")
  transm.sum$f.resid <- log(transm.sum$A.end) - log (g.dat$fluor.) # residuals of feeding rate
  transm.sum$u.LL <- dbinom(x=g.dat$infected, size=1, prob=transm.sum$I.end, log=TRUE) # log likelihood of infection data
  transm.sum$f.LL <- dnorm(x=transm.sum$f.resid, mean=0, sd=sd(transm.sum$f.resid, na.rm=TRUE), log=TRUE) # ll of feeding data
  # throw out nonsensical parameters like negative u or negative infection prevalence
  ifelse(u < 0 | any(transm.sum$I.end<0), 100000, # return a large number to prevent nonsense fit
         -sum(transm.sum$u.LL, na.rm=T) + -sum(transm.sum$f.LL, na.rm=T)) # otherwise sum log liklihoods and send to optimizer
}

summary.lin.f.lin.u <- data.frame(index=clones) 
for (k in 1:length(clones)){ 
  g.dat <- dataset[dataset$index==clones[k] , ]
  print(g.dat$index[1])
  lin.f.lin.u.fit = mle2(minuslogl=lin.f.lin.u.nll, skip.hessian=T, method="Nelder-Mead", 
                         start=list(f.0.hat=12, u=5, w=0.5, alpha=ifelse(g.dat$round[1]==1, -0.5, -3.5)),
                         # trying if/else here, because really hard to get all genotype to converge on sensible params
                         control=list(parscale=c(f.0.hat=12, u=5, w=0.5, alpha=ifelse(g.dat$round[1]==1, 0.5, 3.5)), maxit=10000)) 
  summary.lin.f.lin.u$func[k]="lin.f.lin.u"
  summary.lin.f.lin.u$f.0.hat[k]=as.numeric(coef(lin.f.lin.u.fit)[1])
  summary.lin.f.lin.u$alpha[k]=as.numeric(coef(lin.f.lin.u.fit)[2])/1000 # rescaling back to ml per spore
  summary.lin.f.lin.u$u[k]=as.numeric(coef(lin.f.lin.u.fit)[3])/10000 # rescaling back to u
  summary.lin.f.lin.u$w[k]=as.numeric(coef(lin.f.lin.u.fit)[4])/10000 # rescaling back to ml per spore
  summary.lin.f.lin.u$length[k]=mean(g.dat$mm)
  summary.lin.f.lin.u$time[k]=mean(g.dat$time)
  summary.lin.f.lin.u$loglik[k]=summary(lin.f.lin.u.fit)@m2logL/-2
} 

# all 27 clones take 90 minutes
summary.lin.f.lin.u
write.csv(summary.lin.f.lin.u, file = "summary.lin.f.lin.u.csv")



######################################################################################################
######################                                                ################################
###################     GROUPING BY BLOCK, GENOTYPE, OR ALL TOGETHER              ######################################
######################                                                ######################################################################################################################################
######################################################################################################

# only for the best F and U functions.

# three variations:
# 1) look at each genotype separately (but across all rounds if repeated)
# 2) look at each round spearately (ignoring differences in genotypes within round)
# 3) look at all data together (ignoring differences in rounds and genotypes within rounds)

# together, these results will tell us how much the round-to-round variation matters,
# and how repeatable genotypes are within rounds. 


######################################################################################################
##################################   GENOTYPE FITS (AVERAGE THE REPEATS)    ##########################
######################################################################################################

clones.avg <- unique(dataset$clone) # how many unique?

summary.clones.avg <- data.frame(index=clones.avg) 
for (k in 1:length(clones.avg)){ 
  g.dat <- dataset[dataset$clone==clones.avg[k] , ]
  print(g.dat$clone[1])
  exp.f.exp.u.fit = mle2(minuslogl=exp.f.exp.u.nll, skip.hessian=T, method="Nelder-Mead", 
                         start=list(f.0.hat=10, alpha=-2, u=10, w=-3),
                         control=list(parscale=c(f.0.hat=10, alpha=2, u=10, w=3), maxit=10000)) 
  summary.clones.avg$func[k]="exp.f.exp.u"
  summary.clones.avg$f.0.hat[k]=as.numeric(coef(exp.f.exp.u.fit)[1])
  summary.clones.avg$alpha[k]=as.numeric(coef(exp.f.exp.u.fit)[2])/1000 # rescaling back to ml per spore
  summary.clones.avg$u[k]=as.numeric(coef(exp.f.exp.u.fit)[3])/10000 # rescaling back to u
  summary.clones.avg$w[k]=as.numeric(coef(exp.f.exp.u.fit)[4])/10000 # rescaling back to ml per spore
  summary.clones.avg$length[k]=mean(g.dat$mm)
  summary.clones.avg$time[k]=mean(g.dat$time)
  summary.clones.avg$loglik[k]=summary(exp.f.exp.u.fit)@m2logL/-2
} 

summary.clones.avg

write.csv(summary.clones.avg, file = "summary.clones.avg.csv")


######################################################################################################
############################   ROUND FITS (IGNORE GENOTYPES WITHIN ROUND)    #########################
######################################################################################################

rounds <- unique(dataset$round) # how many unique?

summary.round <- data.frame(index=rounds) 
for (k in 1:length(rounds)){ 
  g.dat <- dataset[dataset$round==rounds[k] , ]
  print(g.dat$round[1])
  exp.f.exp.u.fit = mle2(minuslogl=exp.f.exp.u.nll, skip.hessian=T, method="Nelder-Mead", 
                         start=list(f.0.hat=10, alpha=-2, u=10, w=-3),
                         control=list(parscale=c(f.0.hat=10, alpha=2, u=10, w=3), maxit=10000)) 
  summary.round$func[k]="exp.f.exp.u"
  summary.round$f.0.hat[k]=as.numeric(coef(exp.f.exp.u.fit)[1])
  summary.round$alpha[k]=as.numeric(coef(exp.f.exp.u.fit)[2])/1000 # rescaling back to ml per spore
  summary.round$u[k]=as.numeric(coef(exp.f.exp.u.fit)[3])/10000 # rescaling back to u
  summary.round$w[k]=as.numeric(coef(exp.f.exp.u.fit)[4])/10000 # rescaling back to ml per spore
  summary.round$length[k]=mean(g.dat$mm)
  summary.round$time[k]=mean(g.dat$time)
  summary.round$loglik[k]=summary(exp.f.exp.u.fit)@m2logL/-2
} 

summary.round

write.csv(summary.round, file = "summary.round.csv")


######################################################################################################
#######################################   SPECIES LEVEL FIT    #######################################
######################################################################################################

summary.sp=data.frame(index="species") # just picking the best model to compare 1 vs. all genotypes.
g.dat=dataset # since not looping through anymore, "g.dat" is the entire dataset

exp.f.exp.u.fit = mle2(minuslogl=exp.f.exp.u.nll, skip.hessian=T, method="Nelder-Mead", 
                       start=list(f.0.hat=10, alpha=-2, u=10, w=-3),
                       control=list(parscale=c(f.0.hat=10, alpha=2, u=10, w=3), maxit=10000)) 

summary.sp$func[1]="exp.f.exp.u"
summary.sp$f.0.hat[1]=as.numeric(coef(exp.f.exp.u.fit)[1])
summary.sp$alpha[1]=as.numeric(coef(exp.f.exp.u.fit)[2])/1000 # rescaling back to ml per spore
summary.sp$u[1]=as.numeric(coef(exp.f.exp.u.fit)[3])/10000 # rescaling back to u
summary.sp$w[1]=as.numeric(coef(exp.f.exp.u.fit)[4])/10000 # rescaling back to ml per spore
summary.sp$length[1]=mean(g.dat$mm)
summary.sp$time[1]=mean(g.dat$time)
summary.sp$loglik[1]=summary(exp.f.exp.u.fit)@m2logL/-2

summary.sp

write.csv(summary.sp, file = "summary.sp.csv")


#######################################################################################################
#######################################################################################################
#######################################################################################################

# read summary fits for each: 
summary.con.f.con.u <- read.csv("summary.con.f.con.u.csv")
summary.exp.f.con.u <- read.csv("summary.exp.f.con.u.csv")
summary.con.f.exp.u <- read.csv("summary.con.f.exp.u.csv")
summary.exp.f.exp.u <- read.csv("summary.exp.f.exp.u.csv")

summary.exp.f.lin.u <- read.csv("summary.exp.f.lin.u.csv")
summary.lin.f.exp.u <- read.csv("summary.lin.f.exp.u.csv")
summary.lin.f.lin.u <- read.csv("summary.lin.f.lin.u.csv")
summary.con.f.lin.u <- read.csv("summary.con.f.lin.u.csv")
summary.lin.f.con.u <- read.csv("summary.lin.f.con.u.csv")

summary.clones.avg <- read.csv("summary.clones.avg.csv")
summary.round <- read.csv("summary.round.csv")
summary.sp = read.csv("summary.sp.csv")

# for AIC table:
sum(summary.exp.f.exp.u$loglik)
sum(summary.exp.f.con.u$loglik)
sum(summary.con.f.exp.u$loglik)
sum(summary.con.f.con.u$loglik)

sum(summary.exp.f.lin.u$loglik)
sum(summary.lin.f.exp.u$loglik)
sum(summary.lin.f.lin.u$loglik)
sum(summary.con.f.lin.u$loglik)
sum(summary.lin.f.con.u$loglik)

sum(summary.clones.avg$loglik)
sum(summary.round$loglik)
sum(summary.sp$loglik)


# save as master summary:
fits.summary=rbind(summary.con.f.con.u, summary.exp.f.con.u, 
                   summary.con.f.exp.u, summary.exp.f.exp.u,
                   summary.lin.f.exp.u, summary.exp.f.lin.u,
                   summary.lin.f.con.u, summary.con.f.lin.u,
                   summary.lin.f.lin.u)

write.csv(fits.summary, file = "fits.summary.csv")

