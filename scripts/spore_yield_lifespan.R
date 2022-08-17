##Spore yield ~ lifespan figure

library(here)

source(here("base","src.R"))

lifespan <- readRDS(here("processed_data", "lifespan.rds"))
spores <- readRDS(here("processed_data", "spore_yield.rds"))

lifespan %<>% filter(exposed == 1 & inf == 1)

spores_lifespan <- left_join(lifespan, spores, by = c("temp", "resource", "species"))

spores_lifespan_fig <- spores_lifespan %>% filter(temp %in% const_temp & species =="D") %>% 
  ggplot(.,aes(x=mean_span, y=log_mean_yield, color = temp, shape = resource)) + 
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = log_mean_yield-log_se, ymax = log_mean_yield+log_se, width = 0.3)) +
  geom_errorbarh(aes(xmin = mean_span-se.x, xmax = mean_span+se.x)) +
  scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
  labs(x = "Average Lifespan (days)", y = "ln(Average Spore Yield)", color = "Temp", shape = "Resource") +
  proj_theme

ggsave("spores_lifespan.png", spores_lifespan_fig, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))
