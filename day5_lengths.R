#day5_length

library(here)
library(tidyverse)
library(magrittr)

lengths <- read_csv(here("raw_data", "day5_length.csv"))

lengths %<>% mutate(mm = raw_meas*17.86/1000)
#at default magnification (5.6x), 1 unit is equal to 17.86 micron

length_summ <- lengths %>% 
  group_by(temp_id, resource) %>% 
  summarize(mean_mm = mean(mm),
            var = var(mm),
            sd = sd(mm),
            se = sd(mm)/sqrt(n()))

plot <- length_summ %>% filter(temp_id %in% c(15, 20, 25)) %>% 
  ggplot(.,aes(x=resource, y=mean_mm, color = temp_id)) +
  geom_point(size=4) +
  scale_color_manual(values = c("#FFC107", "#1E88E5", "#D81B60")) +
  geom_linerange(aes(ymin=mean_mm-se, ymax=mean_mm+se)) + 
  theme_bw()


ggsave(here("day5_lengths.png"), plot = plot)
