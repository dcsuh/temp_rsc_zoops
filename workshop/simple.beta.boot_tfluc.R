
############################################
### Estimating Beta from Infection Assay ###
############################################

library(here)
library(tidyverse)
library(magrittr)
library(bbmle)
dataset <- read.csv(here("workshop/post_death_meas.csv")) 
set.seed(8878896)

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

beta.summary

saveRDS(beta.summary, file = here("processed_data","beta_summary.rds"))


beta.summary %>% filter(species=="D" & trt %in% c("2V 1 D", "6V 1 D", "14V 1 D", "20 1 D")) %>% ggplot(.,aes(x=trt, y = beta.est)) +
  geom_point(aes(), position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = beta.025, ymax = beta.975), position = position_dodge(width = 0.5)) +
  theme_minimal()

beta.summary %>% ggplot(.,aes(x=trt, y = beta.est, color = temp)) +
  geom_point(aes(), position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = beta.025, ymax = beta.975), position = position_dodge(width = 0.5)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5))



