#estimating density of susceptible hosts


library(here)

source(here("base","src.R"))




lt.summary_factors <- readRDS(here("processed_data", "lt_summary.rds"))
#this data includes bootstrapped values for host intrinsic growth rate (r)
#alternatively, we can bootstrap this when we bootstrap for host density but for now we will just use the already-bootstrapped values


#parameters used
#d = death rate (taken from bootstrap estimate)
#b = birth rate (taken from bootstrap estimate)
#K = carrying capacity (assumed constant)
#r = resource growth rate (assumed constant)
#f = foraging rate (assumed constant)

r <- 0.9
K <- 250
f <- 2

density.summary <- lt.summary_factors %>% mutate(R = S.d/S.b,
                                                 S = (r*(1-((S.d/S.b)/K)))/f)

beta.summary <- readRDS(here("processed_data", "beta_summary.rds"))

density.summary %<>% mutate(temp=temp_id,
                            species_ID = species,
                            species = gsub("daphnia", "D",
                                      gsub("cerio", "C", species_ID)))

r0.summary <- left_join(density.summary, beta.summary)

spores <- readRDS(here("processed_data", "spore_yield.rds"))
spores$resource <- as.numeric(as.character(spores$resource))
spores$temp <- as.character(spores$temp)

r0.summary %<>% left_join(.,spores, by = join_by(species,resource,temp))

#parameters used for R0 include
#S = susceptible host density
#B = transmission parameter
#A = spore yield
#M = degradation rate

M <- 0.2

r0.summary %<>% mutate(R.naught = (S*beta.est*mean_yield)/M)

r0.summary %>% filter(temp %in% const_temp, inf_status=="I") %>% ggplot(.,aes(x=ID.x, y=R.naught, fill = temp)) + 
  geom_col() + theme(axis.text.x = element_text(angle = 75, hjust = 1))
