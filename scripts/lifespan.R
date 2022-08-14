##lifespan figure

library(here)

source(here("base","src.R"))

lifespan <- readRDS(here("processed_data", "lifespan.rds"))


lifespan %>% filter(temp %in% const_temp & species == "D") %>%
  ggplot(.,aes(x=resource, y = as.numeric(mean_span),color = inf_status)) + 
  geom_point(aes(), position = position_dodge(width = 0.5)) + 
  geom_line(aes(group = interaction(temp, inf_status)), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = mean_span-se, ymax = mean_span+se), position = position_dodge(width = 0.5)) +
  facet_wrap(. ~ temp, nrow=1) +
  scale_color_manual(values = c("red", "blue")) + 
  labs(y = "Average Lifespan (Days)", x = "Resource Concentration mgC/L", color = "Infection\nStatus" ) +
  theme_bw()
