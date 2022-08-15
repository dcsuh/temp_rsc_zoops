##spore yield figure

library(here)

source(here("base","src.R"))

spores <- readRDS(here("processed_data", "spore_yield.rds"))

spores_fig <- spores %>% filter(temp %in% const_temp & species =="D") %>%  
  ggplot(., aes(x=resource, y = log_mean_yield)) + 
  geom_point(aes(), position = position_dodge(width = 0.5)) + 
  geom_line(aes(group = temp), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = log_mean_yield-log_se, ymax = log_mean_yield+log_se), position = position_dodge(width = 0.5)) +
  facet_wrap(vars(temp)) + 
  labs(y = "ln(Average Spore Yield)", x = "Resource Concentration (mgC/L)") +
  proj_theme

spores_fig_20_25 <- spores %>% filter(temp %in% c(20,25) & species =="D") %>%  
  ggplot(., aes(x=resource, y = log_mean_yield)) + 
  geom_point(aes(), position = position_dodge(width = 0.5)) + 
  geom_line(aes(group = temp), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = log_mean_yield-log_se, ymax = log_mean_yield+log_se), position = position_dodge(width = 0.5)) +
  facet_wrap(vars(temp)) + 
  labs(y = "ln(Average Spore Yield)", x = "Resource Concentration (mgC/L)") +
  proj_theme


ggsave("spore_yield.png", spores_fig, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))
ggsave("spore_yield_20_25.png", spores_fig_20_25, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))
