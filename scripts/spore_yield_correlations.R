##Correlation of spore yield with body size and lifespan

library(here)

source(here("base","src.R"))

fitness <- read_csv(here("raw_data/main_fitness_data.csv"))
mort <- readRDS(here("raw_data/main_mort_edit.rds"))

#change classes as necessary
mort$resource <- as.factor(mort$resource)
fitness$resource <- as.factor(fitness$resource)


#remove males, missing, removed, and killed by pipette
mort %<>% filter(!male %in% c(1) & !missing %in% c(1) & species == "D")

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

##Spore Yield
spores <- mort %>% drop_na(inf) %>% filter(exposed==1) %>%
  mutate(spore_conc = ((spore_RAW/8)*10000))
#this doubles the spore count if the water was added because that would make it twice as dilute
#spore yield is calculated by dividing total counted spores by 8 (number of cells used in hemocytomer) and then multiplying by 10000 to get spores/mL
spores %<>% mutate(spore_yield = ifelse(spore_water_added==1,spore_conc*0.5,spore_conc*0.25))


##Lifespan
lifespan <- mort %>% filter(!removed %in% c(1) & !KBP %in% c(1))
lifespan %<>% mutate(birth_day = ifelse(species == "D","4/5/22","4/6/22"))
lifespan %<>% transmute(tube=tube, ID=ID, birthday=mdy(birth_day), deathday=mdy(mortality_day))
lifespan %<>% mutate(span = deathday - birthday)

spores %<>% left_join(.,lifespan)

spores %<>% filter(inf==1 & !is.na(span) & 
                     !is.na(spore_yield) & 
                     !is.na(length) &
                     temp %in% const_temp &
                     species =="D")


spores %>% ggplot(.,aes(x = length, y = log(spore_yield))) + 
  geom_point(aes(color = temp, shape = resource), size = 3) +
  scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
  geom_smooth(method = "lm") +
  labs(x = "Length (mm)", y = "ln(Spore Yield)") +
  proj_theme

spores %>% ggplot(.,aes(x = span, y = log(spore_yield))) + 
  geom_point(aes(color = temp, shape = resource), size = 3) +
  scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
  geom_smooth(method = "lm") +
  labs(x = "Lifespan (days)", y = "ln(Spore Yield)") +
  proj_theme

spores %>% ggplot(.,aes(x = span, y = log(spore_yield))) + 
  geom_point(aes(color = temp, shape = resource), size = 3) +
  scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
  geom_smooth(method = "lm") +
  labs(x = "Lifespan (days)", y = "ln(Spore Yield)") +
  facet_wrap(vars(temp)) +
  proj_theme


