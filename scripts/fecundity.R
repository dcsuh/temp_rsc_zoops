##Fecundity figure

library(here)

source(here("base","src.R"))

fecundity <- readRDS(here("processed_data", "fecundity.rds"))

fecundity_fig <- fecundity %>% filter(temp %in% const_temp & species == "D") %>%
  ggplot(.,aes(x=resource, y = as.numeric(mean_births),color = as.factor(inf))) + 
  geom_point(aes(), position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = mean_births-se, ymax = mean_births+se), position = position_dodge(width = 0.5)) +
  facet_wrap(. ~ temp, nrow=1) +
  scale_color_manual(values = c("blue", "red")) + 
  labs(y = "Average Fecundity (total births)", x = "Resource Concentration mgC/L", color = "Infection\nStatus" ) +
  proj_theme

ggsave("fecundity.png", fecundity_fig, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))
