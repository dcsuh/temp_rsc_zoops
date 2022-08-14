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
  geom_errorbar(aes(ymin = log_mean_yield-log_se, ymax = log_mean_yield+log_se, width = 0.01)) +
  geom_errorbarh(aes(xmin = mean_length-se.x, xmax = mean_length+se.x)) +
  scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
  labs(x = "Average Length (mm)", y = "ln(Average Spore Yield)") +
  theme_bw()
  
