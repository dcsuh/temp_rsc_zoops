library(here)

source(here("base","src.R"))

mort <- readRDS(here("processed_data", "mortality.rds"))
treatment_factors <- readRDS(here("processed_data", "treatment_factors.rds"))

#effects of resource and temperature on length
#still need to figure out transformation to get raw numbers from scope to actual length

length <- mort %>% group_by(ID) %>% summarize(n=n(), mean_length = mean(length, na.rm = T), 
                                              var = var(length, na.rm = T), 
                                              se = sqrt(var(length,na.rm = T)/n()))
length %<>% left_join(.,treatment_factors)

body_size <- length %>% ggplot(., aes(x=ID, y=mean_length, fill=as.factor(resource))) +
  geom_col(position="dodge") + 
  geom_linerange(aes(ymin = mean_length-se, ymax = mean_length+se)) +
  labs(title = "Average Total Length w/ SE", fill = "Resource mgC/L", x = "ID", y="Average Length") + proj_theme

save(body_size, file = here("figures","body_size.RData"))
ggsave("body_size.png", body_size, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))
