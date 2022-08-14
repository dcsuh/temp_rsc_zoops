##Spore yield ~ lifespan figure

library(here)

source(here("base","src.R"))

lifespan <- readRDS(here("processed_data", "lifespan.rds"))
spores <- readRDS(here("processed_data", "spore_yield.rds"))

lifespan %<>% filter(exposed == 1 & inf == 1)

spores_lifespan <- left_join(lifespan, spores, by = c("temp", "resource", "species"))

spores_lifespan %>% filter(temp %in% const_temp & species =="D") %>% 
  ggplot(.,aes(x=mean_span, y=log_mean_yield, color = temp)) + 
  geom_point() +
  geom_errorbar(aes(ymin = log_mean_yield-log_se, ymax = log_mean_yield+log_se, width = 0.01)) +
  geom_errorbarh(aes(xmin = mean_span-se.x, xmax = mean_span+se.x)) +
  scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
  labs(x = "Average Lifespan (days)", y = "ln(Average Spore Yield)") +
  theme_bw()
