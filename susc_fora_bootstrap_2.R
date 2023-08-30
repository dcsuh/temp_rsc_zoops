#Daniel Suh
#August 16, 2023

#estimate per-spore susceptibility using infection and foraging rate data
#beta can be estimaed using infection data from a life table infection assay
#foraging rate is measured separately
#beta is per-spore susceptibility * foraging rate
#so if we can estimate beta, then we can tease out per-spore susceptibility because we have foraging rate as well

#there seems to be two ways to do this. One way is to estimate beta first and then break it down into its components using the mean of those estimated values (i.e., a constant). The other way is to bootstrap both at the same time so that we are pulling beta from a bootstrapped distribution (?) rather than treating it as a constant. I think this would make sense if we had foraging rates associated with the same individuals that we are measuring infection for but we are not. In this case, I think it may be safe to just treat beta as a constant and then generate our bootstrapped distribution for per-spore susceptibility using our sampled foraging rate data.

#this is round 2. I will try to estimate per-spore susceptibility when also estimating beta. The results from the other method don't look right. Not sure yet why but let me try this first

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

foraging %<>% mutate(rate = rate*1440) #convert rate from ml/min to ml/day


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
dataset %<>% filter(temp %in% c(15, 20, 25)) %>% filter(species == "D")


dataset$trt <- paste(dataset$temp, dataset$resource, dataset$species, sep = "_")
foraging %<>% filter(trt == "trt") %>% mutate(ID = paste(temp, resource, "D", sep = "_"))
IDs <- unique(foraging$ID)
trts <- unique(dataset$trt)
trts
beta.summary <- data.frame(trt = trts)
beta.summary$ID <- beta.summary$trt


for(i in 1:length(trts)){
  dsub <- dataset[dataset$trt==trts[i],]
  fora_dsub <- foraging[foraging$ID==IDs[i],]
  beta.summary$temp[i] <- dsub$temp[1]
  beta.summary$resource[i] <- dsub$resource[1]
  beta.summary$species[i] <- dsub$species[1]
  beta.summary$trt[i] <- dsub$trt[1]
  beta.summary$prev[i] <- sum(dsub$inf)/nrow(dsub)
  
  beta.est <- mle2(uninf ~ dbinom(size=1, prob=exp(-beta*spore_exposure*1000*time)),
                   start=list(beta = 0.000001), data=dsub, control=list(parscale = c(beta = 0.000001)),
                   skip.hessian=F,method="L-BFGS-B", lower=.00000001, upper=0.0001)
  beta.summary$beta.est[i] <- coef(beta.est)
  beta.summary$fora.est[i] <- mean(fora_dsub$rate)
  beta.summary$susc.est[i] <- beta.summary$beta.est[i]/beta.summary$fora.est[i]
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
  dsub <- dataset[dataset$ID==IDs[i],]
  fora_dsub <- foraging[foraging$ID==IDs[i],]
  beta.list <- data.frame(beta=numeric(iterations),
                          susc=numeric(iterations),
                          fora=numeric(iterations)) #create empty list to fill in w bootstrap
  for(j in 1:iterations){
    boot.sample <- dsub[sample.int(nrow(dsub),size=nrow(dsub), replace=TRUE),]
    fora.boot.sample <- fora_dsub[sample.int(nrow(fora_dsub),size=nrow(fora_dsub), replace=TRUE),]
    beta.boot <- mle2(uninf ~ dbinom(size=1, prob=exp(-beta*spore_exposure*1000*time)),
                      start=list(beta = 0.000001), data=boot.sample, control=list(parscale = c(beta = 0.000001)),
                      skip.hessian=F,method="L-BFGS-B", lower=.00000001, upper=0.0001)
    beta.list$beta[j] <- coef(beta.boot)[1]
    beta.list$fora[j] <- mean(fora.boot.sample$rate)
    beta.list$susc[j] <- beta.list$beta[j]/beta.list$fora[j]
  }
  beta.summary$beta.025[i] <- quantile(beta.list$beta, probs=seq(0.025, 0.975, 0.95))[1] # lower 95% confidence interval
  beta.summary$beta.975[i] <- quantile(beta.list$beta, probs=seq(0.025, 0.975, 0.95))[2] # upper 95% confidence interval
  beta.summary$fora.025[i] <- quantile(beta.list$fora, probs=seq(0.025, 0.975, 0.95))[1]
  beta.summary$fora.975[i] <- quantile(beta.list$fora, probs=seq(0.025, 0.975, 0.95))[2]
  beta.summary$susc.025[i] <- quantile(beta.list$susc, probs=seq(0.025, 0.975, 0.95))[1]
  beta.summary$susc.975[i] <- quantile(beta.list$susc, probs=seq(0.025, 0.975, 0.95))[2]
}

beta.summary %<>% mutate(beta.prod = susc.est*fora.est,
                         beta.diff = beta.est-beta.prod)

beta.summary %>% ggplot(.,aes(x=ID, y=beta.est)) + geom_point() + geom_pointrange(aes(ymin = beta.025, ymax = beta.975))

beta.summary %>% ggplot(.,aes(x=ID, y=fora.est)) + geom_point() + geom_pointrange(aes(ymin = fora.025, ymax = fora.975))

beta.summary %>% ggplot(.,aes(x=ID, y=susc.est)) + geom_point() + geom_pointrange(aes(ymin = susc.025, ymax = susc.975))

beta.summary %>% ggplot(.,aes(x=ID, y=beta.prod)) + geom_point()

beta.summary %>% ggplot(.,aes(x=ID, y=beta.diff)) + geom_point()
