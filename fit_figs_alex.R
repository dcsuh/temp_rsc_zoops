
library(tidyverse)
library(magrittr)

seq_data <- readRDS("/Users/danielsuh/Documents/GitHub/temp_rsc_zoops/rate_fit.rds")

prevalence <- readRDS("/Users/danielsuh/Documents/GitHub/temp_rsc_zoops/processed_data/prevalence.rds")

theme_set(theme_bw(base_size=18))

symbol_size <- 3

#size-adjusted foraging rate
seq_data %>% dplyr::select(resource, temp, rate_len_mean, rate_len_mean_se) %>% na.omit() %>% 
  ggplot(., aes(x=as.factor(resource), y=rate_len_mean, color=as.factor(temp))) +
  geom_point(size=symbol_size) +
  geom_linerange(aes(ymin=rate_len_mean-rate_len_mean_se,
                     ymax=rate_len_mean+rate_len_mean_se)) +
  scale_color_manual(values = c("#FFC107", "#1E88E5", "#D81B60")) +
  labs(x="Resource", y="Foraging rate (mL/day)", color="Temperature") +
  guides(color = guide_legend(reverse=T))

#prevalence dot plots
prevalence %>% filter(temp %in% c(15, 20, 25) & species =="D") %>%
  ggplot(.,aes(x=resource, y = prev, group = temp, color = temp)) +
  geom_point(size=symbol_size) +
  geom_linerange(aes(ymax=conf$upper, ymin=conf$lower)) +
  scale_color_manual(values = c("#FFC107", "#1E88E5", "#D81B60")) + 
  labs(y = "Probability of Infection", x = "Resource Concentration (mgC/L)", fill = "Temperature") +
  guides(color = guide_legend(reverse=T))

#prevalence bar plots
prevalence %>% filter(temp %in% c(15, 20, 25) & species =="D") %>%
  ggplot(.,aes(x=resource, y = prev, group = temp, color = temp)) +
  geom_bar(stat = "identity", aes(fill = temp), width = 0.5, position = position_dodge(width = 0.5), color = "black") +
  geom_linerange(aes(ymax=conf$upper, ymin=conf$lower), position = position_dodge(width = 0.5), color = "black") +
  scale_fill_manual(values = c("#FFC107", "#1E88E5", "#D81B60")) + 
  labs(y = "Probability of Infection", x = "Resource Concentration (mgC/L)", fill = "Temperature") +
  guides(color = guide_legend(reverse=T))

#generates plot with resource on x-axis
seq_plot <- function(rate){
  seq_data %>% filter(temp %in% c(15,20,25)) %>%
    ggplot(.) +
    geom_line(aes(x=resource, 
                  y=!!sym(rate),
                  group=temp_factor,
                  color=as.factor(temp_factor))) +
    geom_point(aes(x=resource, 
                   y=rate_len_mean, 
                   group = as.factor(temp), 
                   color = as.factor(temp),
                   shape = as.factor(temp)), size=symbol_size) + 
    geom_linerange(aes(x=resource, 
                       ymin=rate_len_mean-rate_len_mean_se,
                       ymax=rate_len_mean+rate_len_mean_se,
                       color=as.factor(temp))) +
    scale_color_manual(values = c("#FFC107", "#1E88E5", "#D81B60")) +
    labs(x="Resource", y="Foraging rate (mL/day)", color="Temperature", shape="Temperature", title = "") +
    guides(color = guide_legend(reverse=T), shape = guide_legend(reverse=T))
}

#generates plot with temp on x-axis
seq_plot_temp <- function(rate){
  seq_data %>% filter(resource %in% c(0.1, 0.5, 1)) %>%
    ggplot(.) +
    geom_line(aes(x=temp, 
                  y=!!sym(rate),
                  group=resource,
                  color=as.factor(resource))) +
    geom_point(aes(x=temp, 
                   y=rate_len_mean, 
                   group = as.factor(resource), 
                   color = as.factor(resource),
                   shape = as.factor(resource)), size=symbol_size) +
    geom_linerange(aes(x=temp, 
                       ymin=rate_len_mean-rate_len_mean_se,
                       ymax=rate_len_mean+rate_len_mean_se,
                       color=as.factor(resource))) + 
    scale_color_manual(values = c("#CABD88", "#CA5216", "#65320D")) +
    labs(x="Temperature", y="Foraging rate (mL/day)", color="", shape="", title = "")
}

seq_plot("m2_rate") #temperature-dependent only
seq_plot("m3_rate") #resource-dependent only
seq_plot("m4_rate") #temp and resource-dependent (full model)
seq_plot("m6c_rate") #full model with handling time scaled by exponential effect of temperature

seq_plot_temp("m2_rate")
seq_plot_temp("m3_rate")
seq_plot_temp("m4_rate")
seq_plot_temp("m6c_rate")

