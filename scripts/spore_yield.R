##spore yield figure

library(here)

source(here("base","src.R"))

spores <- readRDS(here("processed_data", "spore_yield.rds"))

spores %>% filter(temp %in% const_temp & species =="D") %>%  
  ggplot(., aes(x=resource, y = log_mean_yield, color = T)) + 
  geom_point(aes(), position = position_dodge(width = 0.5)) + 
  geom_line(aes(group = temp), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = log_mean_yield-log_se, ymax = log_mean_yield+log_se), position = position_dodge(width = 0.5)) +
  facet_wrap(vars(temp)) + 
  scale_color_manual(values = "red") + 
  labs(y = "Average Spore Yield", x = "Resource Concentration mgC/L") +
  theme_bw() 
