#estimating density of susceptible hosts


library(here)
library(tidyverse)
library(magrittr)
library(bbmle)


dataset <- read.csv(here("raw_data/main_fitness_edit.csv"))


set.seed(8878896)

#DCS Note: if using read_csv so that dataset is a tibble then dates are properly encoded
#however, this somehow changes how some of the NaN's are handled downstream (in the little.r.calculator) which breaks the function
#if NaN's are manually changed to 0's then this then alters results

dataset %<>% mutate(birthdate = ifelse(species == "daphnia", "4/5/22", "4/6/22"), 
                    lifespan = as.numeric(mdy(final_date) - mdy(birthdate))) %>% filter(is.na(male))

#dataset <- read.csv("/Users/dcsuh/Desktop/Traits1-3.csv", na.strings=c("NA", "NULL"))

dataset %<>% mutate(inf_status = replace_na(inf, 0), dead = ifelse(is.na(REMOVED) & is.na(KBP), 1, 0)) %>% mutate(ID = paste(species, temp_id, resource, inf_status, sep = "_"))

treatment_factors <- dataset %>% select(ID, species, temp_id, mean_temp, temp_var, resource, inf_status) %>% distinct()

# this dataset includes results of a life table experiment that I did in the Hall 
# lab. We manipulated spore exposure dose and host genotype, and conducted the 
# experiment over 3 temporal blocks. here, I'll just walk through how to estimate
# little r, little b, and little d for each genotype
treatments <- unique(dataset$ID) 
length(treatments) # how many unique? # this tells you how big to make the summary 
# dataframe where we will store all of the results. we will be looping through. 

lt.summary_factors <- readRDS(here("processed_data", "lt_summary.rds"))
#this data includes bootstrapped values for host intrinsic growth rate (r)
#alternatively, we can bootstrap this when we bootstrap for host density but for now we will just use the already-bootstrapped values


# dataset$uninf <- 1-dataset$inf # could do it either way
# dataset$time <- 1 # best guess at duration of exposure?
# dataset <- dataset[!is.na(dataset$uninf),] # remove NAs
# dataset <- dataset[dataset$spore_exposure > 0,] # can't estimate beta if no exposure!
# dataset$spore_exposure <- as.numeric(dataset$spore_exposure) 
# 
# # this is the function that estimates beta. i'll put it in a couple loops below. 
# # right here it is estimating a single beta for the entire dataset (not very useful,
# # but it gives me a sense of what order of magnitude to start the optimizer on)
# beta.est <- mle2(uninf ~ dbinom(size=1, prob=exp(-beta*spore_exposure*1000*time)), 
#                  start=list(beta = 0.000001), data=dataset, control=list(parscale = c(beta = 0.000001)),
#                  skip.hessian=F,method="L-BFGS-B", lower=.00000001, upper=0.0001)
# print(coef(beta.est))
# 
# # will loop through all treatments:
# dataset$trt <- paste(dataset$temp, dataset$resource, dataset$species, sep = "_")
# trts <- unique(dataset$trt)
# trts
# beta.summary <- data.frame(trt = trts)
# 
# for(i in 1:length(trts)){
#   dsub <- dataset[dataset$trt==trts[i],]
#   beta.summary$temp[i] <- dsub$temp[1]
#   beta.summary$resource[i] <- dsub$resource[1]
#   beta.summary$species[i] <- dsub$species[1]
#   beta.summary$trt[i] <- dsub$trt[1]
#   beta.summary$prev[i] <- sum(dsub$inf)/nrow(dsub)
#   
#   beta.est <- mle2(uninf ~ dbinom(size=1, prob=exp(-beta*spore_exposure*1000*time)),
#                    start=list(beta = 0.000001), data=dsub, control=list(parscale = c(beta = 0.000001)),
#                    skip.hessian=F,method="L-BFGS-B", lower=.00000001, upper=0.0001)
#   beta.summary$beta.est[i] <- coef(beta.est)
# }
