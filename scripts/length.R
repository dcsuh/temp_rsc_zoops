##Body Size figure

library(here)

source(here("base","src.R"))

length <- readRDS(here("processed_data", "length.rds"))

length %>% filter(temp %in% const_temp & species =="D") %>% 
  ggplot(., aes(x=resource, y = mean_length, color = inf_status)) + 
  geom_point(aes(), position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = mean_length-se, ymax = mean_length+se), position = position_dodge(width = 0.5)) + 
  facet_wrap(vars(temp)) +
  scale_color_manual(values = c("red", "blue")) + 
  labs(y = "Average Length (mm)", x = "Resource Concentration mgC/L", color = "Infection\nStatus" ) +
  theme_bw()

