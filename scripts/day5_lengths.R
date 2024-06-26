#day5_length

library(here)
library(tidyverse)
library(magrittr)

theme_set(theme_bw(base_size = 12))

lengths <- read_csv(here("raw_data", "day5_length.csv")) #life table

fora_data <- readRDS(here("processed_data","foraging_raw.rds")) #foraging assay

fora_lengths <- fora_data %>% filter(temp %in% c(15, 20, 25)) %>% dplyr::select(temp, resource, mm) %>% mutate(expt = "fora")

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


#ggsave(here("day5_lengths.png"), plot = plot)


fora_summ <- fora_lengths %>% 
  filter(!is.na(mm)) %>%
  group_by(temp, resource) %>% 
  summarize(mean_mm = mean(mm),
            var = var(mm),
            sd = sd(mm),
            se = sd(mm)/sqrt(n()))

fora_lengths <- fora_summ %>% ggplot(.,aes(x=resource, y=mean_mm, color = as.factor(temp))) +
  geom_point(size=3) +
  scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
  geom_linerange(aes(ymin=mean_mm-se, ymax=mean_mm+se)) +
  scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + 
  labs(x="Resource (mgC/L)", y="Length(mm)", color="Temperature", title = "") + 
  theme_bw(base_size = 10) + 
  theme(legend.position = "none", title = element_blank())

fora_lengths <- fora_summ %>% ggplot(.,aes(x=resource, y=mean_mm, color = as.factor(temp))) +
  geom_point(size=3) +
  #scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
  #scale_color_manual(values = c("#fde725", "#21918c", "#440154")) +
  scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
  #scale_color_manual(values = c("#1E88E5", "#004D40", "#D81B60")) +
  geom_linerange(size = 1,
                 aes(ymin=mean_mm-1.96*se, ymax=mean_mm+1.96*se)) +
  scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + 
  scale_y_continuous(breaks = c(0.8, 1.0, 1.2)) + 
  guides(color = guide_legend(reverse=T)) + 
  labs(x="Resource (mgC/L)", y="Length(mm)", color="Temperature", title = "") + 
  theme_bw(base_size = 10) + 
  theme(legend.position = "none", title = element_blank())
  #theme(title = element_blank())
  

ggsave(here("workshop", "figures","fora_lengths.png"),
       fora_lengths,
       width = 3,
       height = 3,
       units = "in")


#Lengths across experiment and treatment


life_lengths <- lengths %>% 
  filter(temp_var == 0) %>%
  dplyr::select(temp_id, resource, mm) %>% 
  mutate(temp = temp_id, expt = "inf") %>% 
  dplyr::select(-c(temp_id))


data <- rbind(fora_lengths, life_lengths)

data %>% filter(!is.na(mm)) %>% 
  ggplot(., aes(x=as.factor(resource), 
                y=mm, 
                color = as.factor(temp), 
                fill = as.factor(temp),
                shape = expt)) + 
  geom_jitter(size = 2) + 
  scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
  scale_shape_manual(values = c(1, 10)) + 
  labs(x = "Resource Level (mgC/L)", y = "Body Length (mm)", fill = "Temp", color = "Temp", shape = "Expt") +
  theme_minimal()

size_boxplot_15 <- data %>% filter(!is.na(mm)) %>% 
  ggplot(., aes(x=as.factor(resource), y=mm, color = expt)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(position = position_jitterdodge(dodge.width=0.8, jitter.width=0.25)) + 
  scale_color_manual(values = c("orange", "purple")) +
  facet_wrap(~temp) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = "#619CFF")) + 
  labs(x = "Resource Level (mgC/L)", y = "Body Length (mm)", color = "Temp")


size_boxplot_20 <- data %>% filter(!is.na(mm)) %>% 
  ggplot(., aes(x=as.factor(resource), y=mm, color = expt)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(position = position_jitterdodge(dodge.width=0.8, jitter.width=0.25)) + 
  scale_color_manual(values = c("orange", "purple")) +
  facet_wrap(~temp) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = "#00BA38")) + 
  labs(x = "Resource Level (mgC/L)", y = "Body Length (mm)", color = "Temp")


size_boxplot_25 <- data %>% filter(!is.na(mm)) %>% 
  ggplot(., aes(x=as.factor(resource), y=mm, color = expt)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(position = position_jitterdodge(dodge.width=0.8, jitter.width=0.25)) + 
  scale_color_manual(values = c("orange", "purple")) +
  facet_wrap(~temp) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = "#F8766D")) + 
  labs(x = "Resource Level (mgC/L)", y = "Body Length (mm)", color = "Temp")



#ggsave(here("workshop", "figures","size_plots_15.png"),
#        size_boxplot_15,
#        width = 9,
#        height = 6,
#        units = "in")
# 
#ggsave(here("workshop", "figures","size_plots_20.png"),
#        size_boxplot_20,
#        width = 9,
#        height = 6,
#        units = "in")
# 
# ggsave(here("workshop", "figures","size_plots_25.png"),
#        size_boxplot_25,
#        width = 9,
#        height = 6,
#        units = "in")


summary(glm(mm ~ as.numeric(temp) + resource + expt, family = "gaussian", data = data))









