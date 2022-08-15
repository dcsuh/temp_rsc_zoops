##prevalence figure

library(here)

source(here("base","src.R"))

prevalence <- readRDS(here("processed_data", "prevalence.rds"))

prob_inf_fig <- prevalence %>% filter(temp %in% const_temp & species =="D") %>%
  ggplot(.,aes(x=resource, y = prev, group = temp, color = temp)) +
  geom_bar(stat = "identity", aes(fill = temp), width = 0.5, position = position_dodge(width = 0.5), color = "black") +
  geom_linerange(aes(ymax=conf$upper, ymin=conf$lower), position = position_dodge(width = 0.5), color = "black") +
  scale_fill_manual(values = c("#619CFF", "#00BA38", "#F8766D")) + 
  labs(y = "Probability of Infection", x = "Resource Concentration (mgC/L)", fill = "Temperature") + 
  proj_theme

ggsave("prob_inf.png", prob_inf_fig, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))

