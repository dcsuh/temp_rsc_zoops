#Processing data for analyses
#Daniel Suh

library(here)

source(here("base","src.R"))

#read data
fitness <- read_csv(here("raw_data", "main_fitness_data.csv"))
mort <- readRDS(here("raw_data", "main_mort_edit.rds"))
data <- read_csv(here("raw_data","foraging.csv"))


# Cleaning ----------------------------------------------------------------



#change classes as necessary
mort$resource <- as.factor(mort$resource)
fitness$resource <- as.factor(fitness$resource)

#remove males and missing
mort %<>% filter(!male %in% c(1) & !missing %in% c(1)) %>% select(-c(male,missing))


##make treatmentID and reorder
#add factor for infected ("I") or uninfected ("U") to make it easier to make treatmentID
mort %<>% mutate(inf_status = ifelse(inf==1, "I", "U")) %>% filter(inf_status %in% c("I","U"))
#make treatment ID (technically this lumps exposed but uninfected in with unexposed and uninfected)
mort %<>% mutate(ID = paste(temp, resource, species, inf_status,sep = "_"))
#I don't know any better way to set this order
newOrder <- c("15_0.1_D_U","15_0.1_D_I","15_0.5_D_U","15_0.5_D_I","15_1_D_U","15_1_D_I",
              "20_0.1_D_U","20_0.1_D_I","20_0.5_D_U","20_0.5_D_I","20_1_D_U","20_1_D_I",
              "25_0.1_D_U","25_0.1_D_I","25_0.5_D_U","25_0.5_D_I","25_1_D_U","25_1_D_I",
              "15_1_C_U","15_1_C_I","20_1_C_U","20_1_C_I","25_1_C_U","25_1_C_I",
              "2V_1_D_U","2V_1_D_I","6V_1_D_U","6V_1_D_I","14V_1_D_U","14V_1_D_I")
#reorder according to newOrder
mort$ID <- factor(mort$ID, levels = newOrder)

#make df for treatment factors
treatment_factors <- mort %>% dplyr::select(ID,temp,resource, species, exposed, inf, inf_status)
treatment_factors %<>% distinct()
#reorder this again
treatment_factors$ID <- factor(treatment_factors$ID, levels = newOrder)


# Prevalence --------------------------------------------------------------



# Prevalence (probability of infection)
prevalence <- mort %>% drop_na(inf) %>% filter(exposed==1) #remove observations with NA for infection status
prevalence %<>% mutate(ID = paste(temp, resource, species,sep = "_"))

#make order for prevalence data
prevOrder <- c("15_0.1_D","15_0.5_D","15_1_D",
               "20_0.1_D","20_0.5_D","20_1_D",
               "25_0.1_D","25_0.5_D","25_1_D",
               "15_1_C","20_1_C","25_1_C",
               "2V_1_D","6V_1_D","14V_1_D")
prevalence$ID <- factor(prevalence$ID, levels = prevOrder)


prev_treatment_factors <- prevalence %>% dplyr::select(ID,temp,resource,species)
prev_treatment_factors %<>% distinct()
prev_treatment_factors$ID <- factor(prev_treatment_factors$ID, levels = prevOrder)

prevalence %<>% 
  group_by(ID) %>% 
  summarize(n = n(), 
            infected = sum(inf), 
            prev = sum(inf)/n(),
            conf = binom.exact(infected, n, conf.level = 0.95))
prevalence %<>% left_join(.,prev_treatment_factors)





# Body Size ---------------------------------------------------------------

length <- mort %>% filter(!removed %in% c(1) & !KBP %in% c(1)) %>% group_by(ID) %>% 
  summarize(n=n(), 
            mean_length = mean(length, na.rm = T), 
            var = var(length, na.rm = T), 
            se = sqrt(var(length,na.rm = T)/n()))
length %<>% left_join(.,treatment_factors)




# Foraging ----------------------------------------------------------------


data %<>% mutate(length = as.numeric(length),
                 start_min = time_in_hr*60+time_in_min,
                 end_min = time_measured_hr*60+time_measured_min,
                 time = end_min-start_min,
                 time_day = end_min/60/24 - start_min/60/24, #convert to day
                 vol = 15) #volume in mL

control <- 
  data %>% 
  filter(trt=="C") %>% 
  group_by(temp, resource, block) %>% 
  summarize(control_read = mean(read, na.rm=T))

data %<>% 
  left_join(., control) %>% 
  filter(trt=="trt")

#this calculation comes from Hite et al. 2020 and assumes a strong linear relationship between fluorescence and algae conc
#perform linear regression on experimental blocks separately
data_A <- data %>% filter(block=="A")
data_B <- data %>% filter(block=="B")
control_coef_A <- summary(lm(control_read~resource, data=data_A))$coef[2]
control_intercept_A <- summary(lm(control_read~resource, data=data_A))$coef[1]
control_coef_B <- summary(lm(control_read~resource, data=data_B))$coef[2]
control_intercept_B <- summary(lm(control_read~resource, data=data_B))$coef[1]

data_A$read_coef <- control_coef_A
data_A$read_intercept <- control_intercept_A
data_B$read_coef <- control_coef_B
data_B$read_intercept <- control_intercept_B

data <- rbind(data_A, data_B)

data %<>% mutate(mm = as.numeric(length)*17.86/1000, 
                 #at default magnification (5.6x), 1 unit in micrometer is equal to 17.86 micron
                 rate = log(control_read/read)*(vol/time), #vol is mL and time is min for a rate of mL/min
                 rate_len = rate/mm^2,
                 amt_rem = read/read_coef, #final resource concentration
                 amt_rem_tot = amt_rem*vol,
                 amt_init = control_read/read_coef,
                 amt_consumed = amt_init - amt_rem, #final resource amount
                 resource_tot = resource*vol,
                 treatment_ID = paste(temp, resource, sep = "_")) 



data_summ <- data %>% 
  group_by(temp, resource) %>% 
  summarize(rate_mean = mean(rate),
            rate_var = var(rate),
            rate_sd = sd(rate),
            rate_se = sd(rate)/sqrt(n()),
            rate_len_mean = mean(rate_len, na.rm=T),
            rate_len_mean_var = var(rate_len, na.rm=T),
            rate_len_mean_sd = sd(rate_len, na.rm=T),
            rate_len_mean_se = sd(rate_len, na.rm=T)/sqrt(sum(!is.na(rate_len))),
            conc_mean = mean(amt_rem),
            conc_var = var(amt_rem),
            conc_sd = sd(amt_rem),
            conc_se = sd(amt_rem)/sqrt(n()),
            amt_mean = mean(amt_rem_tot),
            amt_var = var(amt_rem_tot),
            amt_sd = sd(amt_rem_tot),
            amt_se = sd(amt_rem_tot)/sqrt(n()),
            amt_init_mean = mean(amt_init),
            amt_init_var = var(amt_init),
            amt_init_sd = sd(amt_init),
            amt_init_se = sd(amt_init)/sqrt(n()),
            length_mean = mean(length, na.rm = T),
            length_var = var(length, na.rm = T),
            length_sd = sd(length, na.rm = T),
            length_se = sd(length, na.rm = T)/sqrt(n()),
            mm_mean = mean(mm, na.rm = T),
            time_mean = mean(time_day, na.rm = T)) %>% ungroup()

data_summ %<>% mutate(species = "D",
                      ID = paste(temp,resource, sep = "_"))

#Some of the lengths are missing so we impute these missing values by using the average length from that treatment
mean_length_summ <- data_summ %>% dplyr::select(temp, resource, mm_mean)
data %<>% left_join(., mean_length_summ)
data %<>% mutate(mm = if_else(!is.na(mm), mm, mm_mean))






# please process data before this line 
if(dir.exists(here("processed_data")) == FALSE) {
  message("Welcome! Let's make some room for the processed data.")
  dir.create(here("processed_data")) 
} else {
  message("/processed_data exists! Proceeeding to save.")
}

saveRDS(prevalence, file = here("processed_data","prevalence.rds"))
saveRDS(length, file = here("processed_data","length.rds"))
saveRDS(data, file = here("processed_data", "foraging_raw.rds"))
saveRDS(data_summ, file = here("processed_data", "foraging.rds"))

