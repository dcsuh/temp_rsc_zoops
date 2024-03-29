################################################
################################################
# Code to estimate parameters from life table 
################################################
################################################


library(bbmle) # install this package before loading it the first times
# this package uses optimizers to find the best parameter(s) that maximizes the 
# likelihood of a function that you define (techncially, minimizes the negaive
# log likelihood)


##########################################################
# 1) read in data and acquaint yourself with its structure

library(tidyverse)
library(magrittr)
library(lubridate)
library(here)

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


#########################################################
# 2) some meta-information to now about the experiment:

age.at.start <- 5  #this tells the function what the starting age of animals in the spread sheet is
# day.columns <- which(names(dataset)%in%c("4/10/22",       "4/11/22",       "4/12/22",       "4/13/22",
#                                          "4/14/22",       "4/15/22",       "4/16/22",       "4/17/22",
#                                          "4/18/22",       "4/19/22",       "4/20/22",       "4/21/22",
#                                          "4/22/22",       "4/23/22",       "4/24/22",       "4/25/22",
#                                          "4/26/22",       "4/27/22",       "4/28/22",       "4/29/22",
#                                          "4/30/22",       "5/1/22",        "5/2/22",        "5/3/22",
#                                          "5/4/22",        "5/5/22",        "5/6/22",        "5/7/22",
#                                          "5/8/22",        "5/9/22",        "5/10/22",       "5/11/22",
#                                          "5/12/22",       "5/13/22",       "5/14/22",       "5/15/22",
#                                          "5/16/22",       "5/17/22",       "5/18/22",       "5/19/22",
#                                          "5/20/22",       "5/21/22",       "5/22/22",       "5/23/22"))

day.columns <- which(names(dataset)%in%c("X4.10.22","X4.11.22","X4.12.22","X4.13.22","X4.14.22","X4.15.22","X4.16.22",
                                         "X4.17.22","X4.18.22","X4.19.22","X4.20.22","X4.21.22","X4.22.22","X4.23.22",
                                         "X4.24.22","X4.25.22","X4.26.22","X4.27.22","X4.28.22","X4.29.22","X4.30.22",
                                         "X5.1.22","X5.2.22","X5.3.22","X5.4.22","X5.5.22","X5.6.22","X5.7.22",
                                         "X5.8.22","X5.9.22","X5.10.22","X5.11.22","X5.12.22","X5.13.22","X5.14.22",
                                         "X5.15.22","X5.16.22","X5.17.22","X5.18.22","X5.19.22","X5.20.22","X5.21.22",
                                         "X5.22.22","X5.23.22"))
# which columns in the spreadsheet have entries that correspond to # babies on each day?
max.age <- length(day.columns) + age.at.start - 1 # we use this in the calculation of death
# rate, in case some hosts never died


#########################################################
# 3) define functions that will estimate r and d

Euler = function(r, column, day){
  sum(column*exp(-r*day)) - 1}
# define the euler lotka equation. this would be a complex polynomial and really
# difficult/impossible to solve by hand... so we just ask the computer to do it for us

d.NLL = function(d){
  d <- exp(d) # apply exp transform to ensure positive death rate
  observed <- ddata$dead
  day <- ddata$lifespan
  NLL <- 0
  for (i in 1:length(observed)){
    NLL <- ifelse(observed[i]==1, NLL - log(d*exp(-d*day[i])),
            NLL - (-d*max.age))}
    NLL}
# the if/else statement says lets us leverage data from hosts that died (and which day 
# they died) but also data for hosts that were still alive at the end

little.r.calculator = function(rdat){
  # fx is mean fecundity produced by surviving moms on each day x
  fx <- colMeans(rdat[day.columns], na.rm=T) # note that moms that 
  # had previously died should be blank (or NA) so they are removed from this mean
  
  #DCS note: daphnia and ceriodaphnia started on different dates so they must use different age.at.start
  if(unique(rdat$species=="cerio")){
    age.at.start <- 4
  }
  
    # z is total number of moms at start of the experiment
  z <- nrow(rdat)
  # sx is proportion of original moms surviving to day x
  Sx <- numeric(length(day.columns))
  for (i in day.columns){
    Sx[(i-day.columns[1]+1)] <- length(na.omit(rdat[,i]))/z}
  
  # now combine survivorship and fecundity into an expression to solve for r
  day <- age.at.start:(length(Sx)+age.at.start-1)
  both <- Sx*fx
    for( k in 1:length(both)){
    if(Sx[k] == 0)
      {both[k] = 0}
    else{both[k]=Sx[k]*fx[k]}}
  
  # finally use the function from earlier to find out what value of 'r' solves the euler-lotka
  little.r <<- uniroot(Euler, interval= c(-2, 2), column=both, day=day)$root
  # need to use <<- to define outside of function (global environment)
  little.r 
  }


#########################################################
# 4) loop through clones and calculate for each one:

lt.summary=data.frame(ID=treatments) # will fill out this summary spreadsheet as
# we loop through each clone
iterations = 1000 # how many samples to use in bootstrap

for (j in 1:length(treatments)){
  print(j) # ask R to tell you where it is in the loop 
  #j=1   # if something is broken in your loop in can be helpful to set the value 
  # of j and go through the inside of the loop line by line to figure out what is wrong
  
  clonedata <- dataset[dataset$ID==treatments[j],]
  # subset the data so you are just looking at clone j. this will get copied over each time
  S.data <- clonedata #[clonedata$infected==0,] # just looking at uninfected (for now)... 
  
  little.r.calculator(rdat=S.data) # this function saves r to global environment
    lt.summary$S.r[j] <- little.r # add it to the summary spreadsheet
  
    ddata <- S.data[!is.na(S.data$dead),] # for death data, we need to remove observations 
    # where we don't know whether host died or not (in practice, user error during experiment)
    d <- ifelse(mean(ddata$dead) == 0, 0,   #this in case nothing died..
             exp(coef(mle2(minuslogl=d.NLL, start=list(d=log(0.05)), 
                            data=list(ddata=ddata), skip.hessian=T))))
  lt.summary$S.d[j] = as.numeric(d) # save to summary
  lt.summary$S.b[j] = lt.summary$S.r[j] + lt.summary$S.d[j] # and calculate b as b = r+d
  
  # now create another summary sheet where we'll save all of the bootstrapped
  # parameter estimates.. make it as long as # iterations you decided earlier
  boot.summary <- data.frame(iterations=1:iterations)
  for(k in 1:iterations){
    boot.rows <- sample(1:nrow(S.data), size=nrow(S.data), replace=TRUE)
    S.boot.data <- S.data[boot.rows,]
    little.r.calculator(rdat=S.boot.data)
    boot.summary$r.boot[k] <- little.r
    ddata <- S.boot.data #[!is.na(S.boot.data$dead),]
    d <- ifelse(mean(ddata$dead) == 0, 0,   #this in case nothing died..
                exp(coef(mle2(minuslogl=d.NLL, start=list(d=log(0.05)),
                              data=list(ddata=ddata), skip.hessian=T))))
    boot.summary$d.boot[k] <- as.numeric(d)
    boot.summary$b.boot[k] <- boot.summary$d.boot[k] + boot.summary$r.boot[k]
  }
  
  # save upper and lower 95% CI's from bootstrapped params:
  lt.summary$S.r.975[j] <- quantile(boot.summary$r.boot, probs=seq(0.025, 0.975, 0.95))[2]  
  lt.summary$S.r.025[j] <- quantile(boot.summary$r.boot, probs=seq(0.025, 0.975, 0.95))[1]  
  lt.summary$S.d.975[j] <- quantile(boot.summary$d.boot, probs=seq(0.025, 0.975, 0.95))[2]  
  lt.summary$S.d.025[j] <- quantile(boot.summary$d.boot, probs=seq(0.025, 0.975, 0.95))[1]  
  lt.summary$S.b.975[j] <- quantile(boot.summary$b.boot, probs=seq(0.025, 0.975, 0.95))[2]  
  lt.summary$S.b.025[j] <- quantile(boot.summary$b.boot, probs=seq(0.025, 0.975, 0.95))[1]  
  
}

# take a look at the results!  (you can save as csv if you want)
lt.summary

lt.summary_factors <- left_join(lt.summary, treatment_factors, by = "ID")

lt.summary_factors$ID <- factor(lt.summary_factors$ID, levels = c(   "daphnia_15_0.1_0", "daphnia_15_0.1_1", 
                                                                     "daphnia_15_0.5_0", "daphnia_15_0.5_1", 
                                                                     "daphnia_15_1_0",   "daphnia_15_1_1",
                                                                     "cerio_15_1_0",     "cerio_15_1_1", 
                                                                     "daphnia_20_0.1_0", "daphnia_20_0.1_1", 
                                                                     "daphnia_20_0.5_0", "daphnia_20_0.5_1",
                                                                     "daphnia_20_1_0",   "daphnia_20_1_1",
                                                                     "cerio_20_1_0",     "cerio_20_1_1",
                                                                     "daphnia_2V_1_0",   "daphnia_2V_1_1",  
                                                                     "daphnia_6V_1_0",   "daphnia_6V_1_1",   
                                                                     "daphnia_14V_1_0",  "daphnia_14V_1_1",     
                                                                     "daphnia_25_0.1_0", "daphnia_25_0.1_1", 
                                                                     "daphnia_25_0.5_0", "daphnia_25_0.5_1", 
                                                                     "daphnia_25_1_0",   "daphnia_25_1_1",     
                                                                     "cerio_25_1_0",     "cerio_25_1_1"))
lt.summary_factors$inf_status <- ifelse(lt.summary_factors$inf_status==1,"I","U")

saveRDS(lt.summary_factors, file = here("processed_data","lt_summary.rds"))

lt.summary_factors %>% ggplot(.,aes(x=ID,y=S.r,color=inf_status, shape=species)) + 
  geom_point() + 
  geom_linerange(aes(ymax=S.r.975, ymin=S.r.025)) +
  scale_color_manual(values = c("forestgreen", "red")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=0.5))


lt.summary_factors %>% filter(species=="daphnia", resource == 1, mean_temp==20) %>% ggplot(.,aes(x=ID,y=S.r,color=inf_status)) + 
  geom_point() + 
  geom_linerange(aes(ymax=S.r.975, ymin=S.r.025)) +
  scale_color_manual(values = c("forestgreen", "red")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle=70, vjust = 0.5, hjust=0.5))


lt.summary_factors %>% filter(resource == 1, temp_var == 0) %>% ggplot(.,aes(x=ID,y=S.r,color=inf_status)) + 
  geom_point() + 
  geom_linerange(aes(ymax=S.r.975, ymin=S.r.025)) +
  scale_color_manual(values = c("forestgreen", "red")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle=70, vjust = 0.5, hjust=0.5))


lt.summary_factors %>% filter(species=="daphnia", temp_var == 0) %>% ggplot(.,aes(x=ID,y=S.r,color=inf_status)) + 
  geom_point() + 
  geom_linerange(aes(ymax=S.r.975, ymin=S.r.025)) +
  scale_color_manual(values = c("forestgreen", "red")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle=70, vjust = 0.5, hjust=0.5))
