##Lifespan ~ body size figure

library(here)

source(here("base","src.R"))

length <- readRDS(here("processed_data", "length.rds"))
lifespan <- readRDS(here("processed_data", "lifespan.rds"))

length_span <- left_join(length, lifespan, by = c("ID", "temp", "resource", "species", "inf", "inf_status"))

length_span_fig <- length_span %>% filter(temp %in% const_temp & species =="D") %>% 
  ggplot(.,aes(x=mean_length, y=mean_span, shape = resource, color = temp)) + 
  geom_point(aes(),size = 3) +
  geom_errorbar(aes(ymin = mean_span-se.y, ymax = mean_span+se.y, width = 0.01)) +
  geom_errorbarh(aes(xmin = mean_length-se.x, xmax = mean_length+se.x)) +
  facet_wrap(. ~ inf_status, nrow=1) +
  scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
  labs(x = "Average Length (mm)", y = "Average lifespan (days)") +
  proj_theme

