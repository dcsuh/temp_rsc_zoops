library(here)

source(here("base","src.R"))

lifespan <- readRDS(here("processed_data","lifespan.rds"))

mean_span <- lifespan %>% group_by(ID) %>% summarize(mean_span = mean(span, na.rm=T),
                                                     var = var(span, na.rm = T),
                                                     se = sqrt(var(span, na.rm = T)/n()),
                                                     resource = unique(resource),
                                                     temperature = unique(temperature))

mean_span %>% ggplot(.,aes(x=ID, y=as.numeric(mean_span), fill=as.factor(resource))) +
  geom_col(position="dodge") + 
  geom_linerange(aes(ymin = mean_span-se, ymax = mean_span+se)) +
  geom_text(aes(label=round(as.numeric(mean_span), digits = 2)), vjust = -0.5) +
  labs(title = "Average lifespan with standard error bars", y = "Average lifespan", fill = "Resource mgC/L")

lifespan %>% drop_na() %>% ggplot(.,aes(x=as.factor(temperature), y = as.numeric(span), fill=as.factor(resource))) + 
  geom_boxplot(position="dodge") + 
  labs(title="avearge lifespan", x="temperature", y ="lifespan")

avg_span <- lifespan %>% drop_na() %>% ggplot(.,aes(x=ID, y = as.numeric(span), fill=as.factor(resource))) + 
  geom_boxplot(position = "dodge") + stat_summary(geom="text", fun=quantile, 
                                                  aes(label=sprintf("%1.1f", ..y..), 
                                                      color=factor(ID)), 
                                                  position=position_nudge(y=0.35), size=3.5) + 
  labs(title = "Host lifespan", x = "Treatment", y = "Time to death in days", fill = "Resource mg C/L") +
  guides(color = "none") + proj_theme

ggsave("lifespan.png", avg_span, width = outwidth[1], height = outwidth[1]/golden, unit = "in", path = here("figures"))
