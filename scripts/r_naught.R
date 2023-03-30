
library(here)

source(here("base","src.R"))

lt_full <- readRDS(here("processed_data", "lt_full_summary.rds"))



r_naught_fig <- lt_full %>% filter(species=="D", temp_var == 0) %>% 
  ggplot(.,aes(x=as.factor(resource),y=S.R_naught)) + 
  geom_point(aes(), position = position_dodge(width = 0.5)) + 
  geom_linerange(aes(ymin=S.R_naught.025, ymax=S.R_naught.975)) + 
  facet_wrap(. ~ temp_id, nrow = 1) + 
  labs(y = "R0", x = "Resource Concentration mgC/L" )

ggsave("R0.png", r_naught_fig, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))


