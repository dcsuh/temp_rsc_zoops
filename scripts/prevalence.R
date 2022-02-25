library(here)

source(here("base","src.R"))

mort <- readRDS(here("processed_data", "mortality.rds"))
prevalence <- readRDS(here("processed_data", "prevalence.rds"))
treatment_factors <- readRDS(here("processed_data", "treatment_factors.rds"))

#effects of resource and temperature on prevalence

prev <- prevalence %>% ggplot(.,aes(x=as.factor(temperature),y=prev, fill=as.factor(resource))) +
  geom_col(position="dodge") + 
  labs(title = "Prevalence", fill = "Resource mgC/L", x = "Temperature", y = "Prevalence") + theme_minimal()


#optional set prev threshold

threshold = 3125

prev_adj <- get_threshold(threshold=threshold)

prev_adj_title <- paste("Threshold-adjusted Prevalence: Threshold = ",threshold, " spores", sep="")

adjusted_prev <- prev_adj %>% ggplot(.,aes(x=as.factor(temperature),y=prev, fill=as.factor(resource))) +
  geom_col(position="dodge") + 
  labs(title = prev_adj_title, fill = "Resource mgC/L", x = "Temperature", y = "Prevalence") + theme_minimal()

save(prev, file = here("figures","prev.RData"))
ggsave("prev.png", prev, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))

save(adjusted_prev, file = here("figures","adjusted_prev.RData"))
ggsave("adjusted_prev.png", adjusted_prev, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))