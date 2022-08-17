##Spore yield ~ body size figure

library(here)

source(here("base","src.R"))

length <- readRDS(here("processed_data", "length.rds"))
spores <- readRDS(here("processed_data", "spore_yield.rds"))

length %<>% filter(exposed == 1 & inf == 1)

spores_length <- left_join(length, spores, by = c("temp", "resource", "species"))

spores_length_fig <- spores_length %>% filter(temp %in% const_temp & species =="D") %>% 
  ggplot(.,aes(x=mean_length, y=log_mean_yield, color = temp, shape = resource)) + 
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = log_mean_yield-log_se, ymax = log_mean_yield+log_se, width = 0.01)) +
  geom_errorbarh(aes(xmin = mean_length-se.x, xmax = mean_length+se.x)) +
  scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
  labs(x = "Average Length (mm)", y = "ln(Average Spore Yield)", color = "Temp", shape = "Resource") +
  proj_theme
  
ggsave("spores_length.png", spores_length_fig, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))
