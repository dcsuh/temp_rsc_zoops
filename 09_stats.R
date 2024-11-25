#traditional statistics

library(here)

source(here("base","src.R"))

foraging <- readRDS(here("processed_data", "foraging_raw.rds"))
mort <- readRDS(here("raw_data/main_mort_edit.rds"))


# prep data ---------------------------------------------------------------

mort %<>% filter(exposed==1 & temp %in% const_temp)
foraging %<>% mutate(temp = as.factor(temp),
                     resource = as.factor(resource))



# length ~ temp. + resource -----------------------------------------------

size_lm <- lm(mm ~ as.numeric(resource) + as.numeric(temp), data = foraging)
size_lm_int <- lm(mm ~ as.numeric(resource) * as.numeric(temp), data = foraging) #interaction not significant

summary(size_lm)
summary(size_lm_int)

size_lm_15 <- lm(mm ~ as.numeric(resource), data = foraging %>% filter(temp == 15))
size_lm_20 <- lm(mm ~ as.numeric(resource), data = foraging %>% filter(temp == 20))
size_lm_25 <- lm(mm ~ as.numeric(resource), data = foraging %>% filter(temp == 25))

summary(size_lm_15)
summary(size_lm_20)
summary(size_lm_25)



# resources consumed ~ temp. + resource -----------------------------------

amt_lm <- lm(amt_consumed ~ as.numeric(resource) + as.numeric(temp), data = foraging)
amt_lm_int <- lm(amt_consumed ~ as.numeric(resource) * as.numeric(temp), data = foraging)

summary(amt_lm)
summary(amt_lm_int)



# prevalence ~ temp. + resource -------------------------------------------

inf_glm <- glm(inf ~ as.numeric(temp) + as.numeric(resource), family = binomial, data = mort)
inf_glm_int <- glm(inf ~ as.numeric(temp) * as.numeric(resource), family = binomial, data = mort)

summary(inf_glm)
summary(inf_glm_int)


inf_glm_15 <- glm(inf ~ as.numeric(resource), 
                  family = binomial, data = mort %>% filter(temp == 15))
inf_glm_20 <- glm(inf ~ as.numeric(resource), 
                  family = binomial, data = mort %>% filter(temp == 20))
inf_glm_25 <- glm(inf ~ as.numeric(resource), 
                  family = binomial, data = mort %>% filter(temp == 25))

summary(inf_glm_15)
summary(inf_glm_20)
summary(inf_glm_25)

