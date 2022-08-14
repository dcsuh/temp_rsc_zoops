library(here)

source(here("base","src.R"))

fitness <- read_csv(here("raw_data/main_fitness_data.csv"))
mort <- readRDS(here("raw_data/main_mort_edit.rds"))

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


##Prevalence
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


##Spore Yield
spores <- mort %>% drop_na(inf) %>% filter(exposed==1) %>%
  mutate(spore_yield = ((spore_RAW*(spore_water_added+1))/8)*10000) 
#this doubles the spore count if the water was added because that would make it twice as dilute
#spore yield is calculated by dividing total counted spores by 8 (number of cells used in hemocytomer) and then multiplying by 10000 to get spores/mL
spores %<>% mutate(log_yield = log(spore_yield),ID = paste(temp, resource, species,sep = "_"))

spores %<>% filter(spore_yield>0) %>% group_by(ID) %>% 
  summarize(mean_yield = mean(spore_yield, na.rm=T), 
            # NAs removed but there should be no NAs after the filter
            var = var(spore_yield, na.rm = T), 
            se = sqrt(var(spore_yield, na.rm = T)/n()),
            log_mean_yield = mean(log_yield, na.rm=T), 
            # NAs removed but there should be no NAs after the filter
            log_var = var(log_yield, na.rm = T), 
            log_se = sqrt(var(log_yield, na.rm = T)/n()))
spores %<>% left_join(.,prev_treatment_factors)
spores$ID <- factor(spores$ID, levels = prevOrder)
spores$temp <- factor(spores$temp, levels = c("15", "20", "25", "2V", "6V", "14V"))


##Body Size
length <- mort %>% filter(!removed %in% c(1) & !KBP %in% c(1)) %>% group_by(ID) %>% 
  summarize(n=n(), 
            mean_length = mean(length, na.rm = T), 
            var = var(length, na.rm = T), 
            se = sqrt(var(length,na.rm = T)/n()))
length %<>% left_join(.,treatment_factors)


##Lifespan
lifespan <- mort %>% filter(!removed %in% c(1) & !KBP %in% c(1))
lifespan %<>% mutate(birth_day = ifelse(species == "D","4/5/22","4/6/22"))
lifespan %<>% transmute(tube=tube, ID=ID, birthday=mdy(birth_day), deathday=mdy(mortality_day))
lifespan %<>% mutate(span = deathday - birthday)

lifespan %<>% group_by(ID) %>% summarize(mean_span = mean(span, na.rm=T),
                                                     var = var(span, na.rm = T),
                                                     se = sqrt(var(span, na.rm = T)/n()))

lifespan %<>% left_join(.,treatment_factors)

##Average fecundity
day_columns <- c("4/10/22", "4/11/22", "4/12/22", "4/13/22",
                 "4/14/22", "4/15/22", "4/16/22", "4/17/22",
                 "4/18/22", "4/19/22", "4/20/22", "4/21/22",
                 "4/22/22", "4/23/22", "4/24/22", "4/25/22",
                 "4/26/22", "4/27/22", "4/28/22", "4/29/22",
                 "4/30/22", "5/1/22", "5/2/22", "5/3/22",
                 "5/4/22", "5/5/22", "5/6/22", "5/7/22",
                 "5/8/22", "5/9/22", "5/10/22", "5/11/22",
                 "5/12/22", "5/13/22", "5/14/22", "5/15/22",
                 "5/16/22", "5/17/22", "5/18/22", "5/19/22",
                 "5/20/22", "5/21/22", "5/22/22", "5/23/22")
fecundity <- inner_join(fitness,mort)
day_columns <- c(which(names(fecundity) %in% day_columns))
fitness_mat <- fecundity %>% select(day_columns)
fecundity %<>% select(tube, ID, temp, resource, species, male, REMOVED, KBP, inf_status)

birth_sums <- rowSums(fitness_mat, na.rm = T)
fecundity$tot_births <- birth_sums

fecundity %<>% filter(is.na(REMOVED) & is.na(KBP)) %>% 
  group_by(temp, resource, species, inf_status) %>% 
  summarize(mean_births = mean(tot_births, na.rm=T),
            var = var(tot_births, na.rm = T),
            se = sqrt(var(tot_births, na.rm = T)/n()))
fecundity$species <- factor(fecundity$species,  levels = c("D", "C"))


# please process data before this line 
if(dir.exists(here("processed_data")) == FALSE) {
  message("Welcome! Let's make some room for the processed data.")
  dir.create(here("processed_data")) 
} else {
  message("/processed_data exists! Proceeeding to save.")
}

saveRDS(prevalence, file = here("processed_data","prevalence.rds"))
saveRDS(spores, file = here("processed_data","spore_yield.rds"))
saveRDS(length, file = here("processed_data","length.rds"))
saveRDS(lifespan, file = here("processed_data","lifespan.rds"))
saveRDS(fecundity, file = here("processed_data","fecundity.rds"))

