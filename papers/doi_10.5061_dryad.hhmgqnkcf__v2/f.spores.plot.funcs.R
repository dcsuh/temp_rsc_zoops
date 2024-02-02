setwd("C:/Users/straussa/Documents/Research/Hall Lab/f spores/data")
library(deSolve)
library(Hmisc)
dataset <- read.csv("Traits1-3.csv")
dataset <- dataset[dataset$clone != "Cerio" , ] # different species
dataset <- dataset[dataset$clone != "Down.282" , ] # really low sample size both times; don't trust
dataset$round <- as.factor(dataset$round)
dataset$clone <- as.factor(dataset$clone)
dataset$index <- paste(dataset$round,dataset$clone) # unique round-clone combo

dataset[dataset$death.age<11 & !is.na(dataset$death.age) , ]$infected <- NA # list the early deaths as NA
dataset[dataset$death.age<11 & !is.na(dataset$death.age) , ]$uninfected <- NA

clones <- unique(dataset$index)
zs <- unique(dataset$exposure)
inf.data <- dataset[!(is.na(dataset$infected)) , ]
f.data <- dataset[!(is.na(dataset$ss.f)) , ]
d.data <- dataset[dataset$infected==0 , ]
colnames(dataset)


################################################################################################


# this function saves prevalence and mean foraging data for each clone at each spore level
data.summary <- data.frame(round=numeric(length(clones)*length(zs)))
for(j in 1:length(clones)){
  f.clonedata <- f.data[f.data$index==clones[j] , ] # subset by spore level (because different foraging predictions)
  inf.clonedata <- inf.data[inf.data$index==clones[j] , ] # subset by spore level (because different foraging predictions)
  d.clonedata <- d.data[d.data$index==clones[j] , ] # including animals still alive
  for (i in 1:length(zs)){
    data.summary$round[i+((j-1)*4)] <- f.clonedata$round[1]
    data.summary$clone[i+((j-1)*4)] <- as.character(f.clonedata$clone[1])
    z.f.data <- f.clonedata[f.clonedata$exposure==zs[i] , ]
    z.inf.data <- inf.clonedata[inf.clonedata$exposure==zs[i] , ]
    z.d.data <- d.clonedata[d.clonedata$exposure==zs[i] , ]
    data.summary$exposure[i+((j-1)*4)] <- z.inf.data$exposure[1]
    prev <- sum(z.inf.data$infected)/length(z.inf.data$infected)
    data.summary$prev[i+((j-1)*4)] <- prev
    data.summary$prev.stderr[i+((j-1)*4)] <- sqrt(prev*(1-prev)/nrow(z.inf.data))
    data.summary$f.mean[i+((j-1)*4)] <- mean(z.f.data$f)
    data.summary$f.stderr[i+((j-1)*4)] <- sd(z.f.data$f)/sqrt(length(z.f.data$f))
    died.early <- z.d.data[z.d.data$death.age<11 , ]
    data.summary$prop.early.dead[i+((j-1)*4)] <- nrow(died.early)/nrow(z.d.data)
    data.summary$int[i+((j-1)*4)] <- coef(lm(mult.ed~exposure,data=d.clonedata))[1]
    data.summary$slope[i+((j-1)*4)] <- coef(lm(mult.ed~exposure,data=d.clonedata))[2]
    }}
data.summary$index <- paste(data.summary$round,data.summary$clone) # unique round-clone combo
data.summary$round <- as.factor(data.summary$round)


############################################################################################

# need to impose a minimum feeding rate:
# genotype averages:
hist(data.summary$f.mean)
min(data.summary$f.mean, na.rm=T) # so minimum is 2.8 ml/day, across all clones. 
quantile(data.summary$f.mean, na.rm=T, probs=seq(0, 0.1, 0.025))
f.min <- 3.6 # this is the bottom 2.5% 


############################################################################################

# what percent drop in feeding rate for each clone?

min(1-(data.summary[data.summary$exposure==393,]$f.mean / data.summary[data.summary$exposure==0,]$f.mean))
median(1-(data.summary[data.summary$exposure==393,]$f.mean / data.summary[data.summary$exposure==0,]$f.mean))
max(1-(data.summary[data.summary$exposure==393,]$f.mean / data.summary[data.summary$exposure==0,]$f.mean))

############################################################################################

library(car) # for type III sum squares aov
library(lme4)
library(lmerTest) # for p values in lme4 (but not perfect)

# analysis of deviance tables:
fmod <- glm(ss.f ~ clone * exposure + clone * round, data=f.data)
Anova(fmod, type="II")  
shapiro.test(resid(fmod))

imod <- glm(infected ~ clone * exposure + clone * round, data=inf.data, 
            family=binomial(link="logit"))
Anova(imod, type="II") 

################################################################################################

# now, need functions to plot PREDICTED foraging and infection prevalence

fits <- read.csv("fits.summary.csv")
fits <- cbind(fits, t(matrix(unlist(strsplit(as.character(fits$index), " ")), nrow=2)))
colnames(fits)[12] <- "round"
colnames(fits)[13] <- "clone"

vol <- 15
z.seq <- seq(from=0,to=393, length.out=20) # sequence to plot for x axis


###################
# CONSTANT F AND U
###################

cfcu <-fits[fits$func=="con.f.con.u" , ]
cfcu.f <- matrix(ncol=length(z.seq), nrow=length(cfcu$func))
rownames(cfcu.f) <- cfcu$index
cfcu.i <- cfcu.f; cfcu.u <- cfcu.f; cfcu.z <- cfcu.f
for (i in 1:length(cfcu$func)){cfcu.f[i,]=cfcu$f.0.hat[i]*(cfcu$length[i]^2)}
for (i in 1:length(cfcu$func)){cfcu.u[i,]=cfcu$u[i]}
for (i in 1:length(cfcu$func)){
  cfcu.i[i,]=1 - exp((cfcu$u[i])*z.seq*vol*(exp(-(cfcu$f.0.hat[i]/vol)*(cfcu$length[i]^2)*cfcu$time[i])-1))}
for (i in 1:length(cfcu$func)){cfcu.z[i,]=cfcu$f.0.hat[i]*(cfcu$length[i]^2)*z.seq}

##############
# EXP F CON U:
##############

efcu <- fits[fits$func=="exp.f.con.u" , ]
efcu.sim <- function(times,y,params){
  Z = y[1]; A = y[2]; S = y[3]; I = y[4]
  with(as.list(params),{
    dZ = -(S+I)*Z*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol) # exponential f; con u
    dA = -(S+I)*A*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol) # exponential f; con u
    dS = -u*S*Z*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol) # exponential f; con u
    dI = u*S*Z*max(f.0.hat*(length^2)*exp(alpha*Z*(length^2)), f.min/vol) # exponential f; con u
    res = c(dZ, dA, dS, dI)
    list(res)})}
efcu.f=matrix(ncol=length(z.seq), nrow=length(efcu$func))
rownames(efcu.f) <- efcu$index
efcu.i <- efcu.f; efcu.u <- efcu.f; efcu.z <- efcu.f
transm <- function(Z, length, time, f.0.hat, alpha, u){ # function to vectorize the simulation
  sim <- rk(y=c(Z=Z, A=200, S=1, I=0), func=efcu.sim, method="ode45", rtol=1e-4, atol=1e-4,
            parms=c(f.0.hat=f.0.hat/vol, alpha=alpha/vol, u=u, length=length),
            times=seq(from=0, to=time, by=time))}
for(i in 1:nrow(efcu)){
  out <- mapply(transm, Z=z.seq*vol, length=efcu$length[i], time=efcu$time[i],
                f.0.hat=efcu$f.0.hat[i], alpha=efcu$alpha[i], u=efcu$u[i])
  efcu.f[i, ] <- pmax(efcu$f.0.hat[i]*(efcu$length[i]^2)*exp(efcu$alpha[i]*z.seq*(efcu$length[i]^2)), f.min)
  efcu.i[i, ] <- out[10, ] # final infection prevalence
  efcu.u[i, ] <- efcu$u[i]
  efcu.z[i, ] <- pmax(efcu$f.0.hat[i]*(efcu$length[i]^2)*exp(efcu$alpha[i]*z.seq*(efcu$length[i]^2)), f.min)*z.seq
  }


##############
# CON F EXP U:
##############

cfeu <- fits[fits$func=="con.f.exp.u" , ]
cfeu.sim <- function(times,y,params){
  Z = y[1]; A = y[2]; S = y[3]; I = y[4]
  with(as.list(params),{
    dZ = -(S+I)*Z*f.0.hat*(length^2) # specify form: con f; exp u
    dA = -(S+I)*A*f.0.hat*(length^2) # specify form: con f; exp u
    dS = -u*exp(w*Z*f.0.hat*(length^2))*S*Z*f.0.hat*(length^2) # specify form: con f; exp u
    dI = u*exp(w*Z*f.0.hat*(length^2))*S*Z*f.0.hat*(length^2) # specify form: con f; exp u
    res = c(dZ, dA, dS, dI)
    list(res)})}
cfeu.f=matrix(ncol=length(z.seq), nrow=length(cfeu$func))
rownames(cfeu.f) <- cfeu$index
cfeu.i <- cfeu.f; cfeu.u <- cfeu.f; cfeu.z <- cfeu.f
transm <- function(Z, length, time, f.0.hat, u, w){ # function to vectorize the simulation
  sim <- rk(y=c(Z=Z, A=200, S=1, I=0), func=cfeu.sim, method="ode45", rtol=1e-4, atol=1e-4,
            parms=c(f.0.hat=f.0.hat/vol, u=u, w=w, length=length),
            times=seq(from=0, to=time, by=time))}
for(i in 1:nrow(cfeu)){
  out <- mapply(transm, Z=z.seq*vol, length=cfeu$length[i], time=cfeu$time[i],
                f.0.hat=cfeu$f.0.hat[i], u=cfeu$u[i], w=cfeu$w[i])
  cfeu.f[i, ] <- cfeu$f.0.hat[i]*(cfeu$length[i]^2)
  cfeu.i[i, ] <- out[10, ] # final infection prevalence
  cfeu.u[i, ] <- cfeu$u[i]*exp(cfeu$w[i]*z.seq*cfeu$f.0.hat[i]*(cfeu$length[i]^2))
  cfeu.z[i, ] <- cfeu$f.0.hat[i]*(cfeu$length[i]^2)*z.seq
  }


#########################
# EXPONENTIAL ALPHA AND U:
#########################

efeu <- fits[fits$func=="exp.f.exp.u" , ]
efeu.sim <- function(times,y,params){  
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
efeu.f=matrix(ncol=length(z.seq), nrow=length(efeu$func))
rownames(efeu.f) <- efeu$index
efeu.i <- efeu.f; efeu.u <- efeu.f; efeu.z <- efeu.f
transm <- function(Z, length, time, f.0.hat, alpha, u, w){ # function to vectorize the simulation
  sim <- rk(y=c(Z=Z, A=200, S=1, I=0), func=efeu.sim, method="ode45", rtol=1e-4, atol=1e-4, 
            parms=c(f.0.hat=f.0.hat/vol, alpha=alpha/vol, u=u, w=w, length=length), 
            times=seq(from=0, to=time, by=time))} 
for(i in 1:nrow(efeu)){
  out <- mapply(transm, Z=z.seq*vol, length=efeu$length[i], time=efeu$time[i], 
                f.0.hat=efeu$f.0.hat[i], alpha=efeu$alpha[i], u=efeu$u[i], w=efeu$w[i])
  efeu.f[i, ] <- pmax(efeu$f.0.hat[i]*(efeu$length[i]^2)*
                       exp(efeu$alpha[i]*(efeu$length[i]^2)*z.seq), f.min)
  efeu.i[i, ] <- out[10, ] # final infection prevalence
  efeu.u[i, ] <- efeu$u[i]*exp(efeu$w[i]*z.seq*pmax(efeu$f.0.hat[i]*(efeu$length[i]^2)*
                                 exp(efeu$alpha[i]*(efeu$length[i]^2)*z.seq), f.min))
  efeu.z[i, ] <- pmax(efeu$f.0.hat[i]*(efeu$length[i]^2)*
                       exp(efeu$alpha[i]*(efeu$length[i]^2)*z.seq), f.min)*z.seq
}


#########################
# EXPONENTIAL F; LINEAR U:
#########################

eflu <- fits[fits$func=="exp.f.lin.u" , ]
eflu.sim <- function(times,y,params){  
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
eflu.f=matrix(ncol=length(z.seq), nrow=length(eflu$func))
rownames(eflu.f) <- eflu$index
eflu.i <- eflu.f; eflu.u <- eflu.f; eflu.z <- eflu.f
transm <- function(Z, length, time, f.0.hat, alpha, u, w){ # function to vectorize the simulation
  sim <- rk(y=c(Z=Z, A=200, S=1, I=0), func=eflu.sim, method="ode45", rtol=1e-4, atol=1e-4, 
            parms=c(f.0.hat=f.0.hat/vol, alpha=alpha/vol, u=u, w=w, length=length), 
            times=seq(from=0, to=time, by=time))} 
for(i in 1:nrow(eflu)){
  out <- mapply(transm, Z=z.seq*vol, length=eflu$length[i], time=eflu$time[i], 
                f.0.hat=eflu$f.0.hat[i], alpha=eflu$alpha[i], u=eflu$u[i], w=eflu$w[i])
  eflu.f[i, ] <- pmax(eflu$f.0.hat[i]*(eflu$length[i]^2)*exp(eflu$alpha[i]*(eflu$length[i]^2)*z.seq), f.min)
  eflu.i[i, ] <- out[10, ] # final infection prevalence
  eflu.u[i, ] <- eflu$u[i]*(1+(eflu$w[i]*z.seq*pmax(eflu$f.0.hat[i]*(eflu$length[i]^2)*exp(eflu$alpha[i]*(eflu$length[i]^2)*z.seq), f.min)))
  eflu.z[i, ] <- pmax(eflu$f.0.hat[i]*(eflu$length[i]^2)*exp(eflu$alpha[i]*(eflu$length[i]^2)*z.seq), f.min)*z.seq
}


##############
# LIN F CON U: 
##############

lfcu <- fits[fits$func=="lin.f.con.u" , ]
lfcu.sim <- function(times,y,params){  
  Z = y[1]; A = y[2]; S = y[3]; I = y[4]
  with(as.list(params),{
    dZ = -(S+I)*Z*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol) # specify form: linear alpha 
    dA = -(S+I)*A*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol) 
    dS = -u*S*Z*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol) 
    dI = u*S*Z*max(f.0.hat*(length^2)*(1+(alpha*Z*(length^2))), f.min/vol)
    res = c(dZ, dA, dS, dI)
    list(res)})}
lfcu.f=matrix(ncol=length(z.seq), nrow=length(lfcu$func)); rownames(lfcu.f) <- lfcu$index
lfcu.i <- lfcu.f; lfcu.u <- lfcu.f; lfcu.z <- lfcu.f
transm <- function(Z, length, time, f.0.hat, alpha, u){ # function to vectorize the simulation
  sim <- rk(y=c(Z=Z, A=200, S=1, I=0), func=lfcu.sim, method="ode45", rtol=1e-4, atol=1e-4, 
            parms=c(f.0.hat=f.0.hat/vol, alpha=alpha/vol, u=u, length=length), 
            times=seq(from=0, to=time, by=time))} 
# 200 is arbitrary fluoresence. just need some positive number for foraging to deplete.
for(i in 1:nrow(lfcu)){
  out <- mapply(transm, Z=z.seq*vol, length=lfcu$length[i], time=lfcu$time[i], 
                f.0.hat=lfcu$f.0.hat[i], alpha=lfcu$alpha[i], u=lfcu$u[i])
  lfcu.f[i, ] <- pmax(lfcu$f.0.hat[i]*(lfcu$length[i]^2)*(1+(lfcu$alpha[i]*z.seq*lfcu$length[i]^2)), f.min)
  lfcu.i[i, ] <- out[10, ] # final infection prevalence
  lfcu.u[i, ] <- lfcu$u[i]
  lfcu.z[i, ] <- pmax(lfcu$f.0.hat[i]*(lfcu$length[i]^2)*(1+(lfcu$alpha[i]*z.seq*lfcu$length[i]^2)), f.min)*z.seq
}


##############
# CON F LIN U: 
##############

cflu <- fits[fits$func=="con.f.lin.u" , ]
cflu.sim <- function(times,y,params){  
  Z = y[1]; A = y[2]; S = y[3]; I = y[4]
  with(as.list(params),{
    dZ = -(S+I)*Z*f.0.hat*(length^2) # specify form: con f
    dA = -(S+I)*A*f.0.hat*(length^2)
    dS = -u*S*Z*f.0.hat*(length^2)*(1+(w*Z*f.0.hat*(length^2))) 
    dI = u*S*Z*f.0.hat*(length^2)*(1+(w*Z*f.0.hat*(length^2)))   
    res = c(dZ, dA, dS, dI)
    list(res)})}
cflu.f=matrix(ncol=length(z.seq), nrow=length(cflu$func))
rownames(cflu.f) <- cflu$index
cflu.i <- cflu.f; cflu.u <- cflu.f; cflu.z <- cflu.f
transm <- function(Z, length, time, f.0.hat, u, w){ # function to vectorize the simulation
  sim <- rk(y=c(Z=Z, A=200, S=1, I=0), func=cflu.sim, method="ode45", rtol=1e-4, atol=1e-4, 
            parms=c(f.0.hat=f.0.hat/vol, u=u, w=w, length=length), 
            times=seq(from=0, to=time, by=time))} 
for(i in 1:nrow(cflu)){
  out <- mapply(transm, Z=z.seq*vol, length=cflu$length[i], time=cflu$time[i], 
                f.0.hat=cflu$f.0.hat[i], u=cflu$u[i], w=cflu$w[i])
  cflu.f[i, ] <- cflu$f.0.hat[i]*(cflu$length[i]^2)
  cflu.i[i, ] <- out[10, ] # final infection prevalence
  cflu.u[i, ] <- cflu$u[i]*(1+(cflu$w[i]*z.seq*cflu$f.0.hat[i]*(cflu$length[i]^2)))
  cflu.z[i, ] <- cflu$f.0.hat[i]*(cflu$length[i]^2)*z.seq
}


#########################
# LINEAR F; EXPONENTIAL U:
#########################

lfeu <- fits[fits$func=="lin.f.exp.u" , ]
lfeu.sim <- function(times,y,params){  
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
lfeu.f=matrix(ncol=length(z.seq), nrow=length(lfeu$func))
rownames(lfeu.f) <- lfeu$index
lfeu.i <- lfeu.f; lfeu.u <- lfeu.f; lfeu.z <- lfeu.f
transm <- function(Z, length, time, f.0.hat, alpha, u, w){ # function to vectorize the simulation
  sim <- rk(y=c(Z=Z, A=200, S=1, I=0), func=lfeu.sim, method="ode45", rtol=1e-4, atol=1e-4, 
            parms=c(f.0.hat=f.0.hat/vol, alpha=alpha/vol, u=u, w=w, length=length), 
            times=seq(from=0, to=time, by=time))} 
for(i in 1:nrow(lfeu)){
  out <- mapply(transm, Z=z.seq*vol, length=lfeu$length[i], time=lfeu$time[i], 
                f.0.hat=lfeu$f.0.hat[i], alpha=lfeu$alpha[i], u=lfeu$u[i], w=lfeu$w[i])
  lfeu.f[i, ] <- pmax(lfeu$f.0.hat[i]*(lfeu$length[i]^2)*(1+(lfeu$alpha[i]*(lfeu$length[i]^2)*z.seq)), f.min)
  lfeu.i[i, ] <- out[10, ] # final infection prevalence
  lfeu.u[i, ] <- lfeu$u[i]*exp(lfeu$w[i]*z.seq*pmax(lfeu$f.0.hat[i]*(lfeu$length[i]^2)*(1+(lfeu$alpha[i]*(lfeu$length[i]^2)*z.seq)), f.min))
  lfeu.z[i, ] <- pmax(lfeu$f.0.hat[i]*(lfeu$length[i]^2)*(1+(lfeu$alpha[i]*(lfeu$length[i]^2)*z.seq)), f.min)*z.seq
}


#########################
# LINEAR F; LINEAR U:
#########################

lflu <- fits[fits$func=="lin.f.lin.u" , ]
lflu.sim <- function(times,y,params){  
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
lflu.f=matrix(ncol=length(z.seq), nrow=length(lflu$func))
rownames(lflu.f) <- lflu$index
lflu.i <- lflu.f; lflu.u <- lflu.f; lflu.z <- lflu.f
transm <- function(Z, length, time, f.0.hat, alpha, u, w){ # function to vectorize the simulation
  sim <- rk(y=c(Z=Z, A=200, S=1, I=0), func=lflu.sim, method="ode45", rtol=1e-4, atol=1e-4, 
            parms=c(f.0.hat=f.0.hat/vol, alpha=alpha/vol, u=u, w=w, length=length), 
            times=seq(from=0, to=time, by=time))} 
for(i in 1:nrow(lflu)){
  out <- mapply(transm, Z=z.seq*vol, length=lflu$length[i], time=lflu$time[i], 
                f.0.hat=lflu$f.0.hat[i], alpha=lflu$alpha[i], u=lflu$u[i], w=lflu$w[i])
  lflu.f[i, ] <- pmax(lflu$f.0.hat[i]*(lflu$length[i]^2)*(1+(lflu$alpha[i]*(lflu$length[i]^2)*z.seq)), f.min)
  lflu.i[i, ] <- out[10, ] # final infection prevalence
  lflu.u[i, ] <- lflu$u[i]*(1+(lflu$w[i]*z.seq*pmax(lflu$f.0.hat[i]*(lflu$length[i]^2)*(1+(lflu$alpha[i]*(lflu$length[i]^2)*z.seq)), f.min)))
  lflu.z[i, ] <- pmax(lflu$f.0.hat[i]*(lflu$length[i]^2)*(1+(lflu$alpha[i]*(lflu$length[i]^2)*z.seq)), f.min)*z.seq
}


######################################################################
# plot:

plot.fit = function(g.index, yax){

 par(mar=c(0,0,5,0.5))
 data <- data.summary[data.summary$index==g.index , ]
 plot(data$exposure,data$f.mean, pch=20, cex=4, xlim=c(-20,420), ylim=c(0,22), # ymax=28 for all; 22 for fig 1
      xlab="", ylab="", xaxt="n", yaxt="n", col="white")
 points(data$exposure,data$f.mean,pch=15,cex=2.5) 
 with(data=data,expr=errbar(data$exposure,data$f.mean, data$f.mean+data$f.stderr, 
                            data$f.mean-data$f.stderr, add=T,pch=1,cap=0))
 if(yax=="yes"){axis(side=2, at=c(0,10,20), cex.axis=1.8)} # cex 1.8 for fig 1; 1.2 for r1-r3
 lines(z.seq, cfcu.f[which(rownames(cfcu.f)==g.index) , ], col="orange", lwd=2, lty=1)  
 lines(z.seq, efeu.f[which(rownames(efeu.f)==g.index) , ], col="deepskyblue", lwd=5, lty=1)  
 lines(z.seq, efcu.f[which(rownames(efcu.f)==g.index) , ], col="dark green", lwd=3, lty=3)  
 lines(z.seq, cfeu.f[which(rownames(cfeu.f)==g.index) , ], col="purple", lwd=3, lty=2)
 #lines(z.seq, eflu.f[which(rownames(eflu.f)==g.index) , ], col="blue", lwd=5, lty=1)  
 #lines(z.seq, cflu.f[which(rownames(cflu.f)==g.index) , ], col="purple", lwd=3, lty=2)  
 #lines(z.seq, lflu.f[which(rownames(lflu.f)==g.index) , ], col="black", lwd=3, lty=2)  
 #lines(z.seq, lfcu.f[which(rownames(lfcu.f)==g.index) , ], col="pink", lwd=3, lty=1)  
 #lines(z.seq, lfeu.f[which(rownames(lfeu.f)==g.index) , ], col="red", lwd=3, lty=1)  
 
 # SHOW SPORES CONSUMPTION RATE?
 # plot(z.seq, z.seq, xlim=c(-20,420), ylim=c(0,5000), xlab="", ylab="", # ymax=0.0015 for fig2; ymax=0.003 for all
 #      xaxt="n", yaxt="n", col="white")
 # #axis(side=2, at=c(0,1000,2000,3000), cex.axis=1.6)
 # #axis(side=1, at=c(0,200,400), cex.axis=1.6)
 # lines(z.seq, efeu.z[which(rownames(efeu.z)==g.index) , ], col="deepskyblue", lwd=5, lty=1)  
 
 par(mar=c(5,0,0.5,0.5))
 plot(data$exposure,data$prev, pch=20, cex=4, xlim=c(-20,420), ylim=c(0,0.7), # ymax=1 for all; 0.7 for fig 1
      xlab="", ylab="", xaxt="n", yaxt="n", col="white")
 points(data$exposure,data$prev,pch=15,cex=2.5) 
 with(data=data,expr=errbar(data$exposure,data$prev, data$prev+data$prev.stderr, 
                            data$prev-data$prev.stderr, add=T,pch=1,cap=0))
 if(yax=="yes"){axis(side=2, at=c(0,0.5,1), cex.axis=1.8)} # cex 1.8 for fig 1; 1.2 for r1-r3
 axis(side=1, at=c(0,200,400), cex.axis=1.8) # cex 1.8 for fig 1; 1.2 for r1-r3
 lines(z.seq, cfcu.i[which(rownames(cfcu.i)==g.index) , ], col="orange", lwd=2, lty=1)  
 lines(z.seq, efeu.i[which(rownames(efeu.i)==g.index) , ], col="deepskyblue", lwd=5, lty=1)  
 lines(z.seq, efcu.i[which(rownames(efcu.i)==g.index) , ], col="dark green", lwd=3, lty=3)
 lines(z.seq, cfeu.i[which(rownames(cfeu.i)==g.index) , ], col="purple", lwd=3, lty=2)  
 #lines(z.seq, eflu.i[which(rownames(eflu.i)==g.index) , ], col="blue", lwd=5, lty=1)  
 #lines(z.seq, cflu.i[which(rownames(cflu.i)==g.index) , ], col="purple", lwd=3, lty=2)  
 #lines(z.seq, lflu.i[which(rownames(lflu.i)==g.index) , ], col="black", lwd=3, lty=2)  
 #lines(z.seq, lfcu.i[which(rownames(lfcu.i)==g.index) , ], col="pink", lwd=3, lty=1)  
 #lines(z.seq, lfeu.i[which(rownames(lfeu.i)==g.index) , ], col="red", lwd=3, lty=1)  
 
 plot(z.seq, z.seq, xlim=c(-20,8000), ylim=c(0,0.003), xlab="", ylab="", 
      xaxt="n", yaxt="n", col="white")
 #abline(h=0)
 axis(side=1, at=c(0,3000,6000), cex.axis=1.8) # cex 1.8 for fig 1; 1.2 for r1-r3
 if(yax=="yes"){axis(side=2, at=c(0,0.002), cex.axis=1.8)} # cex 1.8 for fig 1; 1.2 for r1-r3
 lines(cfcu.z[which(rownames(cfcu.z)==g.index) , ], cfcu.u[which(rownames(cfcu.u)==g.index) , ], col="orange", lwd=2, lty=1)
 lines(efeu.z[which(rownames(efeu.z)==g.index) , ], efeu.u[which(rownames(efeu.u)==g.index) , ], col="deepskyblue", lwd=5, lty=1)
 lines(efcu.z[which(rownames(efcu.z)==g.index) , ], efcu.u[which(rownames(efcu.u)==g.index) , ], col="dark green", lwd=3, lty=3)
 lines(cfeu.z[which(rownames(cfeu.z)==g.index) , ], cfeu.u[which(rownames(cfeu.u)==g.index) , ], col="purple", lwd=3, lty=2)
 #lines(eflu.z[which(rownames(eflu.z)==g.index) , ], eflu.u[which(rownames(eflu.u)==g.index) , ], col="blue", lwd=5, lty=1)
 #lines(cflu.z[which(rownames(cflu.z)==g.index) , ], cflu.u[which(rownames(cflu.u)==g.index) , ], col="purple", lwd=3, lty=2)
 #lines(lflu.z[which(rownames(lflu.z)==g.index) , ], lflu.u[which(rownames(lflu.u)==g.index) , ], col="black", lwd=3, lty=2)
 #lines(lfcu.z[which(rownames(lfcu.z)==g.index) , ], lfcu.u[which(rownames(lfcu.u)==g.index) , ], col="pink", lwd=3, lty=1)
 #lines(lfeu.z[which(rownames(lfeu.z)==g.index) , ], lfeu.u[which(rownames(lfeu.u)==g.index) , ], col="red", lwd=3, lty=1)
 } 

# 3 example genotypes: 

layout(matrix(c(1:9), 3, 3, byrow = FALSE))
par(mar=c(0,0,0.5,0.5),oma=c(0.5,3,0,0)) 
plot.fit("3 Mid.252", yax="yes") # good example of flat foraging and u
plot.fit("3 Mid.244", yax="no") # good example of curvy u
plot.fit("2 Bris.10", yax="no") # good example of exponential depression

# genotypes by round: 

layout(matrix(c(1:27), 3, 9, byrow = FALSE))
par(mar=c(0,0,0.5,0.5),oma=c(0.5,2,0,0)) 
plot.fit("1 A45", yax="yes")
plot.fit("1 Bris.10", yax="no")
plot.fit("1 Bris.112", yax="no")
plot.fit("1 Bris.6", yax="no")
plot.fit("1 Cback.256", yax="no")
plot.fit("1 Dog.4", yax="no") 
plot.fit("1 Mid.276", yax="no")
plot.fit("1 Mid.277", yax="no")
plot.fit("1 War.5", yax="no") 

layout(matrix(c(1:27), 3, 9, byrow = FALSE))
par(mar=c(0,0,0.5,0.5),oma=c(0.5,2,0,0)) 
plot.fit("2 A45", yax="yes")
plot.fit("2 A48", yax="no")
plot.fit("2 Bris.10", yax="no")
plot.fit("2 Dog.4", yax="no")
plot.fit("2 Is.278", yax="no") 
plot.fit("2 Mid.263", yax="no")
plot.fit("2 Mid.273", yax="no") 
plot.fit("3 A43", yax="no")
plot.fit("3 A45", yax="no")

layout(matrix(c(1:27), 3, 9, byrow = FALSE))
par(mar=c(0,0,0.5,0.5),oma=c(0.5,2,0,0)) 
plot.fit("3 Bris.10", yax="yes")
plot.fit("3 Bris.111", yax="no")
plot.fit("3 Bris.112", yax="no")
plot.fit("3 Bris.6", yax="no")
plot.fit("3 Cback.276", yax="no")
plot.fit("3 Mid.244", yax="no") 
plot.fit("3 Mid.252", yax="no") 
plot.fit("3 Mid.281", yax="no")
plot.fit("3 War.5", yax="no")



