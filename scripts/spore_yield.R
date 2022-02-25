library(here)

source(here("base","src.R"))

mort <- readRDS(here("processed_data", "mortality.rds"))
lifespan <- readRDS(here("processed_data","lifespan.rds"))


#Effects of lifespan and length on spore yield

#set threshold for spore yield to be included
threshold = 1000

body_size_spores <- mort %>% filter(spore_yield>threshold) %>% drop_na() %>% 
  ggplot(.,aes(x=length, y=log(spore_yield))) + geom_point() + proj_theme

lifespan_spores <- lifespan %>% filter(spore_yield>threshold) %>% drop_na() %>%
  ggplot(.,aes(x=as.numeric(span), y=log(spore_yield), 
               color=as.factor(temperature), 
               shape = as.factor(resource))) + geom_point(size=5) + proj_theme

save(body_size_spores, file = here("figures","body_size_spores.RData"))
ggsave("body_size_spores.png", body_size_spores, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))

save(lifespan_spores, file = here("figures","lifespan_spores.RData"))
ggsave("lifespan_spores.png", lifespan_spores, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))