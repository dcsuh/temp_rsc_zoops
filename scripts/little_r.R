##Figure for net host reproductive rate

library(here)

source(here("base","src.R"))

lt.summary_factors <- readRDS(here("processed_data", "lt_summary.rds"))

little_r_fig <- lt.summary_factors %>% filter(species=="daphnia", temp_var == 0) %>% 
  ggplot(.,aes(x=as.factor(resource),y=S.r,color=inf_status)) + 
  geom_point(aes(), position = position_dodge(width = 0.5)) + 
  geom_linerange(aes(ymax=S.r.975, ymin=S.r.025), position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("red", "blue")) + 
  facet_wrap(. ~ temp_id, nrow = 1) + 
  labs(y = "Host net reproductive rate", x = "Resource Concentration mgC/L", color = "Infection\nStatus" )

ggsave("little_r.png", little_r_fig, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))

little_r_I <- lt.summary_factors %>% filter(species=="daphnia", temp_var == 0) %>% 
  ggplot(.,aes(x=as.factor(resource),y=S.r,color=inf_status)) + 
  geom_point(aes(), position = position_dodge(width = 0.5)) + 
  geom_linerange(aes(ymax=S.r.975, ymin=S.r.025), position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("red", "white")) + 
  facet_wrap(. ~ temp_id, nrow = 1) + 
  labs(y = "Host net reproductive rate", x = "Resource Concentration mgC/L", color = "Infection\nStatus" )

ggsave("little_r_I.png", little_r_I, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))

little_r_U <- lt.summary_factors %>% filter(species=="daphnia", temp_var == 0) %>% 
  ggplot(.,aes(x=as.factor(resource),y=S.r,color=inf_status)) + 
  geom_point(aes(), position = position_dodge(width = 0.5)) + 
  geom_linerange(aes(ymax=S.r.975, ymin=S.r.025), position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("white", "blue")) + 
  facet_wrap(. ~ temp_id, nrow = 1) + 
  labs(y = "Host net reproductive rate", x = "Resource Concentration mgC/L", color = "Infection\nStatus" )

ggsave("little_r_U.png", little_r_U, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))

