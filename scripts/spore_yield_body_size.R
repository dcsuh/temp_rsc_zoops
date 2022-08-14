##Spore yield ~ body size figure

library(here)

source(here("base","src.R"))

length <- readRDS(here("processed_data", "length.rds"))
spores <- readRDS(here("processed_data", "spore_yield.rds"))

length %<>% filter(exposed == 1 & inf == 1)

spores_length <- left_join(length, spores, by = c("temp", "resource", "species"))

spores_length %>% filter(temp %in% const_temp & species =="D") %>% 
  ggplot(.,aes(x=mean_length, y=log_mean_yield, color = temp)) + 
  geom_point() +
  geom_errorbar(aes(ymin = log_mean_yield-log_se, ymax = log_mean_yield+log_se)) +
  geom_errorbarh(aes(xmin = mean_length-se.x, xmax = mean_length+se.x)) +
  scale_color_manual() +
  labs(x = "Average Length (mm)", y = "ln(Average Spore Yield)")
  

length %>% filter(temp %in% const_temp & species =="D") %>% 
  ggplot(., aes(x=resource, y = mean_length, color = inf_status)) + 
  geom_point(aes(), position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = mean_length-se, ymax = mean_length+se), position = position_dodge(width = 0.5)) + 
  facet_wrap(vars(temp)) +
  scale_color_manual(values = c("red", "blue")) + 
  labs(y = "Average Length (mm)", x = "Resource Concentration mgC/L", color = "Infection\nStatus" ) +
  theme_bw()