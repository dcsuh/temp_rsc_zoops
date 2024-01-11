
library(here)

library(plotly)

source(here("base","src.R"))

lt_full <- readRDS(here("processed_data", "lt_full_summary.rds"))



r_naught_fig <- lt_full %>% filter(species=="D", temp_var == 0) %>% 
  ggplot(.,aes(x=as.factor(resource),y=S.R_naught, color=temp_id)) + 
  geom_point(aes(), position = position_dodge(width = 0.5)) + 
  geom_linerange(aes(ymin=S.R_naught.025, ymax=S.R_naught.975)) + 
  scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
  guides(color="none")+
  facet_wrap(. ~ temp_id, nrow = 1) + 
  labs(y = "R0", x = "Resource Concentration mgC/L" )

ggsave("R0.png", r_naught_fig, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))


r_naught_fig_temp <- lt_full %>% filter(species=="D", temp_var == 0) %>% 
  ggplot(.,aes(x=as.factor(temp_id),y=S.R_naught, color=temp_id)) + 
  geom_point(aes(), position = position_dodge(width = 0.5)) + 
  geom_linerange(aes(ymin=S.R_naught.025, ymax=S.R_naught.975)) + 
  scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) + 
  facet_wrap(. ~ as.factor(resource), nrow = 1) + 
  labs(y = "R0", x = "Temp (C)" )

ggsave("R0_temp.png", r_naught_fig_temp, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))


r_naught_fill <- lt_full %>% filter(species=="D", temp_var == 0) %>% 
  ggplot(.,aes(x=as.factor(temp_id),y=as.factor(resource),fill=S.R_naught)) +
  geom_tile() +
  scale_fill_gradient(low = "orange", high= "red") + 
  labs(x = "Temperature", y = "Resource", fill = "R0")

plot_ly(data=lt_full, x=as.factor(lt_full$temp_id), y=as.factor(lt_full$resource), z=lt_full$S.R_naught,
        type="scatter3d", color=lt_full$species)
