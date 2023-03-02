#estimating density of susceptible hosts


library(here)
library(tidyverse)
library(magrittr)
library(bbmle)



lt.summary_factors <- readRDS(here("processed_data", "lt_summary.rds"))
#this data includes bootstrapped values for host intrinsic growth rate (r)
#alternatively, we can bootstrap this when we bootstrap for host density but for now we will just use the already-bootstrapped values


#parameters used
#d = death rate (taken from bootstrap estimate)
#b = birth rate (taken from bootstrap estimate)
#K = carrying capacity (assumed constant)
#r = resource growth rate (assumed constant)
#f = foraging rate (assumed constant)

r <- 1
K <- 1
f <- 10

lt.summary_factors %<>% mutate(R = S.d/S.b,
                               S = (r*(1-((S.d/S.b)/K)))/f)

