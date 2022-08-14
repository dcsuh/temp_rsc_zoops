library(here)

source(here("base","src.R"))

lifespan <- readRDS(here("processed_data","lifespan.rds"))

mean_span <- lifespan %>% group_by(ID) %>% summarize(mean_span = mean(span, na.rm=T),
                                                     var = var(span, na.rm = T),
                                                     se = sqrt(var(span, na.rm = T)/n()),
                                                     resource = unique(resource),
                                                     temperature = unique(temperature))

avg_span <- lifespan %>% drop_na() %>% ggplot(.,aes(x=ID, y = as.numeric(span), fill=as.factor(resource))) + 
  geom_boxplot(position = "dodge") + stat_summary(geom="text", fun=quantile, 
                                                  aes(label=sprintf("%1.1f", ..y..)), 
                                                  position=position_nudge(y=0.35, x=-0.2), size=3) + 
  labs(title = "Host lifespan", x = "Treatment", y = "Time to death in days", fill = "Resource mg C/L") +
  guides(color = "none") + proj_theme

save(avg_span, file = here("figures","avg_span.RData"))
ggsave("lifespan.png", avg_span, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))
