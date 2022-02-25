library(here)

source(here("base","src.R"))

treatments <- read_csv(here("raw_data", "treatments.csv"))
fitness <- read_csv(here("raw_data", "fitness.csv"))
mort <- read_csv(here("raw_data", "mortality.csv"))

#reorder data
newOrder <- c("1C","1B","1A","2C","2B","2A","3C","3B","3A")
treatments$ID <- factor(treatments$ID, levels = newOrder)
fitness$ID <- factor(fitness$ID, levels = newOrder)
mort$ID <- factor(mort$ID, levels = newOrder)

treatment_factors <- treatments %>% dplyr::select(ID,resource,temperature)

#lifespan
lifespan <- left_join(mort, treatments, by = "ID")
lifespan %<>% select(ID, replicate, birth_date, date)
lifespan %<>% transmute(ID=ID, replicate=replicate, birthday=mdy(birth_date), deathday=mdy(date))
lifespan %<>% mutate(span = deathday - birthday)
lifespan %<>% left_join(.,treatment_factors)


#prevalence
prevalence <- mort %>% drop_na(infection_status) #remove observations with NA for infection status
prevalence %<>% group_by(ID) %>% summarize(n = n(), prev = sum(infection_status)/n())
prevalence %<>% left_join(.,treatment_factors)




# please process data before this line 
if(dir.exists(here("processed_data")) == FALSE) {
  message("Welcome! Let's make some room for the processed data.")
  dir.create(here("processed_data")) 
} else {
  message("/processed_data exists! Proceeeding to save.")
}

# please save data after this line 
saveRDS(lifespan, file = here("processed_data","lifespan.rds"))
saveRDS(fitness, file = here("processed_data","fitness.rds"))
saveRDS(mort, file = here("processed_data","mortality.rds"))
saveRDS(treatments, file = here("processed_data","treatment.rds"))
saveRDS(treatment_factors, file = here("processed_data","treatment_factors.rds"))
saveRDS(prevalence, file = here("processed_data", "prevalence.rds"))


