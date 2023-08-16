#Daniel Suh
#August 16, 2023

#estimate per-spore susceptibility using infection and foraging rate data
#beta can be estimaed using infection data from a life table infection assay
#foraging rate is measured separately
#beta is per-spore susceptibility * foraging rate
#so if we can estimate beta, then we can tease out per-spore susceptibility because we have foraging rate as well

#there seems to be two ways to do this. One way is to estimate beta first and then break it down into its components using the mean of those estimated values (i.e., a constant). The other way is to bootstrap both at the same time so that we are pulling beta from a bootstrapped distribution (?) rather than treating it as a constant. I think this would make sense if we had foraging rates associated with the same individuals that we are measuring infection for but we are not. In this case, I think it may be safe to just treat beta as a constant and then generate our bootstrapped distribution for per-spore susceptibility using our sampled foraging rate data.

library(here)
library(tidyverse)
library(magrittr)
library(bbmle)
dataset <- read.csv(here("workshop/post_death_meas.csv")) 
foraging <- read_csv(here("raw_data","foraging.csv"))
set.seed(8878896)


#calculate foraging rates for each individual
foraging %<>% mutate(length = as.numeric(length),
                 start_min = time_in_hr*60+time_in_min,
                 end_min = time_measured_hr*60+time_measured_min,
                 time = end_min-start_min,
                 vol = 15)

control <- foraging %>% filter(trt=="C") %>% 
  group_by(temp, resource, block) %>% 
  summarize(control_read = mean(read, na.rm=T))

foraging %<>% left_join(., control)

foraging %<>% filter(trt=="trt") %>% 
  mutate(rate = log(control_read/read)*(vol/time),
         amt_time = (control_read - read)/time)


dataset$uninf <- 1-dataset$inf # could do it either way
dataset$time <- 1 # best guess at duration of exposure?
dataset <- dataset[!is.na(dataset$uninf),] # remove NAs
dataset <- dataset[dataset$spore_exposure > 0,] # can't estimate beta if no exposure!
dataset$spore_exposure <- as.numeric(dataset$spore_exposure) 

# this is the function that estimates beta. i'll put it in a couple loops below. 
# right here it is estimating a single beta for the entire dataset (not very useful,
# but it gives me a sense of what order of magnitude to start the optimizer on)
beta.est <- mle2(uninf ~ dbinom(size=1, prob=exp(-beta*spore_exposure*1000*time)), 
                 start=list(beta = 0.000001), data=dataset, control=list(parscale = c(beta = 0.000001)),
                 skip.hessian=F,method="L-BFGS-B", lower=.00000001, upper=0.0001)
print(coef(beta.est))

# will loop through all treatments:
dataset$trt <- paste(dataset$temp, dataset$resource, dataset$species, sep = "_")
trts <- unique(dataset$trt)
trts
beta.summary <- data.frame(trt = trts)

for(i in 1:length(trts)){
  dsub <- dataset[dataset$trt==trts[i],]
  beta.summary$temp[i] <- dsub$temp[1]
  beta.summary$resource[i] <- dsub$resource[1]
  beta.summary$species[i] <- dsub$species[1]
  beta.summary$trt[i] <- dsub$trt[1]
  beta.summary$prev[i] <- sum(dsub$inf)/nrow(dsub)
  
  beta.est <- mle2(uninf ~ dbinom(size=1, prob=exp(-beta*spore_exposure*1000*time)),
                   start=list(beta = 0.000001), data=dsub, control=list(parscale = c(beta = 0.000001)),
                   skip.hessian=F,method="L-BFGS-B", lower=.00000001, upper=0.0001)
  beta.summary$beta.est[i] <- coef(beta.est)
}

warnings()  # these are fine. 
beta.summary  
plot(beta.summary$prev, beta.summary$beta.est) # make sure that beta and prevalence are 
# positively correlated. yup, checks out. 
# note about units: units for time in the spreadsheet were already in day; exposure was per ml,
# but after multiplying by 1000, now per L. Therefore, units on beta end up being L/day. 


# okay, next, I'll bootstrap confidence intervals around each of those. 
iterations <- 1000   # 100 iterations will take ~2 minutes to run through all treatments.
# i would recommend 1,000 or even 10,000 to be really confident on the confidence intervals

for(i in 1:length(trts)){
  print(i)
  dsub <- dataset[dataset$trt==trts[i],]
  beta.list <- data.frame(beta=numeric(iterations)) #create empty list to fill in w bootstrap
  for(j in 1:iterations){
    boot.sample <- dsub[sample.int(nrow(dsub),size=nrow(dsub), replace=TRUE),]
    beta.boot <- mle2(uninf ~ dbinom(size=1, prob=exp(-beta*spore_exposure*1000*time)),
                      start=list(beta = 0.000001), data=boot.sample, control=list(parscale = c(beta = 0.000001)),
                      skip.hessian=F,method="L-BFGS-B", lower=.00000001, upper=0.0001)
    beta.list$beta[j] <- coef(beta.boot)[1]
  }
  beta.summary$beta.025[i] <- quantile(beta.list$beta, probs=seq(0.025, 0.975, 0.95))[1] # lower 95% confidence interval
  beta.summary$beta.975[i] <- quantile(beta.list$beta, probs=seq(0.025, 0.975, 0.95))[2] # upper 95% confidence interval
}

foraging %<>% filter(trt == "trt") %>% mutate(ID = paste(temp, resource, "D", sep = "_"))
beta.summary %<>% filter(temp %in% c(15, 20, 25)) %>% filter(species == "D")
beta.summary$ID <- beta.summary$trt

foraging %<>% left_join(., beta.summary, by = "ID")

foraging %<>% mutate(susc = beta.est/rate)

IDs <- unique(foraging$ID)

fora.summary <- data.frame(ID = IDs)

for(i in 1:length(trts)){
  print(i)
  dsub <- foraging[foraging$ID==IDs[i],]
  fora.list <- data.frame(susc=numeric(iterations),
                          fora=numeric(iterations))
  for(j in 1:iterations){
    boot.sample <- dsub[sample.int(nrow(dsub), size=nrow(dsub), replace = TRUE),]
    fora.list$susc[j] <- mean(boot.sample$susc)
    fora.list$fora[j] <- mean(boot.sample$rate)
  }
  fora.summary$susc.025[i] <- quantile(fora.list$susc, probs=seq(0.025, 0.975, 0.95))[1] 
  fora.summary$susc.975[i] <- quantile(fora.list$susc, probs=seq(0.025, 0.975, 0.95))[2] 
  fora.summary$fora.025[i] <- quantile(fora.list$fora, probs=seq(0.025, 0.975, 0.95))[1] 
  fora.summary$fora.975[i] <- quantile(fora.list$fora, probs=seq(0.025, 0.975, 0.95))[2] 
}

beta.summary
