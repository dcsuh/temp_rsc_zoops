#Processing data for analyses
#Daniel Suh


library(here)

source(here("base","src.R"))

#read data
fitness <- read_csv(here("raw_data/main_fitness_data.csv"))
mort <- readRDS(here("raw_data/main_mort_edit.rds"))


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

##Body Size
length <- mort %>% filter(!removed %in% c(1) & !KBP %in% c(1)) %>% group_by(ID) %>% 
  summarize(n=n(), 
            mean_length = mean(length, na.rm = T), 
            var = var(length, na.rm = T), 
            se = sqrt(var(length,na.rm = T)/n()))
length %<>% left_join(.,treatment_factors)






# please process data before this line 
if(dir.exists(here("processed_data")) == FALSE) {
  message("Welcome! Let's make some room for the processed data.")
  dir.create(here("processed_data")) 
} else {
  message("/processed_data exists! Proceeeding to save.")
}

saveRDS(prevalence, file = here("processed_data","prevalence.rds"))
saveRDS(length, file = here("processed_data","length.rds"))

