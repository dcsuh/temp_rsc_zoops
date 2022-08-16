##Figure for net host reproductive rate

library(here)

source(here("base","src.R"))

lt.summary_factors <- readRDS(here("processed_data", "lt_summary.rds"))

little_r_fig <- lt.summary_factors %>% filter(species=="daphnia", temp_var == 0) %>% 
  ggplot(.,aes(x=as.factor(resource),y=S.r,color=inf_status)) + 
  geom_point(aes(), position = position_dodge(width = 0.5)) + 
  geom_linerange(aes(ymax=S.r.975, ymin=S.r.025), position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("blue", "red")) + 
  facet_wrap(. ~ temp_id, nrow = 1) + 
  labs(y = "Host net reproductive rate", x = "Resource Concentration mgC/L", color = "Infection\nStatus" ) +
  theme_bw()

ggsave("little_r.png", little_r_fig, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))
