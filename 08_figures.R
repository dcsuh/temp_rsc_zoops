#define functions for making figures
#Daniel Suh

library(here)

source(here("base","src.R"))



#read data
f_seq <- readRDS(file=here("processed_data", "seq_data", "foraging_rate_fit_data.rds"))
u_seq <- readRDS(file=here("processed_data", "seq_data", "infection_fit_data.rds"))

data_summ <- readRDS(here("processed_data", "foraging.rds"))



h_ci_seq <- readRDS(here("processed_data", "seq_data", "h_ci_bootstrap.rds"))
f_ci_seq <- readRDS(here("processed_data", "seq_data", "f_ci_bootstrap.rds"))
u_ci_seq <- readRDS(here("processed_data", "seq_data", "u_ci_bootstrap.rds"))

mod_quantiles <- readRDS(here("processed_data", "m2E_bootstrap_quantiles.rds"))



# global formatting -------------------------------------------------------

theme_set(theme_bw(base_size = 12))


# some data prep ----------------------------------------------------------

data_summ %<>% dplyr::select(temp, resource, amt_init_mean, mm_mean, mm_se)
f_seq %<>% left_join(., data_summ)

#f_seq %<>% mutate(amt_diff = resource-amt_rem_mean,
f_seq %<>% mutate(amt_diff = amt_init_mean-amt_rem_mean,
                  m1_diff = resource - m1_R_end/15*1000,
                  m2_diff = resource - m2_R_end/15*1000,
                  m3_diff = resource - m3_R_end/15*1000,
                  m4_diff = resource - m4_R_end/15*1000,
                  m5_diff = resource - m5_R_end/15*1000)


f_seq %<>% arrange(desc(temp))


if(dir.exists(here("figures")) == FALSE) {
  message("Welcome! Let's make some room for figures.")
  dir.create(here("figures")) 
} else {
  message("/figures exists! Proceeeding to save.")
}


# Figure 1: Data summary --------------------------------------------------
length_data <- data_summ %>%
  drop_na() %>%
  ggplot(., aes(x=resource, y=mm_mean, group = as.factor(temp), color = as.factor(temp))) + 
  geom_point(size=3) +
  geom_linerange(size = 1, 
                 aes(ymin=mm_mean-1.96*mm_se,
                     ymax=mm_mean+1.96*mm_se,
                     color=as.factor(temp))) +
  scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
  scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + 
  theme(legend.position = "none", title = element_blank()) + 
  labs(x="Resource (mgC/L)", y="Length (mm)", color="Temperature", title = "")

ggsave(here("figures","01_length_data.png"),
       length_data,
       width = 3,
       height = 3,
       units = "in")

f_data <- f_seq %>% 
  dplyr::select(resource, temp, amt_diff, amt_rem_se, amt_init_mean) %>% 
  drop_na() %>%
  #ggplot(., aes(x=resource, y=amt_diff, group = as.factor(temp), color = as.factor(temp))) + 
  ggplot(., aes(x=amt_init_mean, y=amt_diff, group = as.factor(temp), color = as.factor(temp))) + 
  geom_point(size=3) +
  geom_linerange(size = 1, 
                 aes(ymin=amt_diff-1.96*amt_rem_se,
                     ymax=amt_diff+1.96*amt_rem_se,
                     color=as.factor(temp))) +
  scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
  scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + 
  theme(legend.position = "none", title = element_blank()) + 
  labs(x="Resource (mgC/L)", y="Resource difference (mgC/L)", color="Temperature", title = "")

ggsave(here("figures","01_f_data.png"),
       f_data,
       width = 3,
       height = 3,
       units = "in")


prev_data <- u_seq %>% filter(!is.na(prev)) %>% 
  ggplot(.,aes(x=resource, 
               y=prev, 
               color=as.factor(temp))) + 
  geom_point(size = 3) + 
  geom_linerange(size=1,
                 aes(ymin = conf$lower, ymax=conf$upper, color=as.factor(temp))) +
  #scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
  scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
  scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + 
  scale_y_continuous(breaks = c(0, 0.5, 1.0)) + 
  theme(legend.position = "none", title = element_blank()) + 
  labs(x="Resource (mgC/L)", y="Probability of Infection", color="Temperature", title = "")

ggsave(here("figures","01_prev_data.png"),
       prev_data,
       width = 3,
       height = 3,
       units = "in")



# plotting fn's -----------------------------------------------------------

seq_plot_alt <- function(seq_data, rate){
  seq_data %>% filter(temp %in% c(15,20,25)) %>%
    ggplot(.) +
    geom_line(aes(x=resource, 
                  y=!!sym(rate),
                  group=as.factor(temp),
                  color=as.factor(temp)),
              size = 1) +
    scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
    scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + 
    scale_y_continuous(breaks = c(0, 20, 40), limits = c(0, 45)) + 
    guides(color = guide_legend(reverse=T)) + 
    labs(x="Resource (mgC/L)", y="Foraging rate (mL/day)", color="Temperature", shape="Temperature", title = "")
}

algae_plot_diff <- function(seq_data, var){
  seq_data %>% filter(temp %in% const_temp) %>%
    ggplot(., aes(x=resource, y=!!sym(var), color=as.factor(temp), group = as.factor(temp))) +
    geom_line(size=1) +
    geom_point(aes(x=amt_init_mean,
                   y=amt_diff,
                   group = as.factor(temp),
                   color = as.factor(temp)),
               size=0.5) +
    geom_linerange(aes(x=amt_init_mean,
                       ymin=amt_diff-1.96*amt_rem_se,
                       ymax=amt_diff+1.96*amt_rem_se,
                       color=as.factor(temp))) +
    labs(x = "Resource (mgC/L)",
         y = "Resource difference (mgC/L)",
         color = "Temp") +
    scale_x_continuous(breaks = c(0.1, 0.5, 1.0), limits = c()) + 
    scale_y_continuous(breaks = c(0, 0.2, 0.4), limits = c(-0.05, 0.425)) + 
    guides(color = guide_legend(reverse=T)) + 
    scale_color_manual(values = c("#FFC107", "#21918c", "#440154"))
}

spores_plot <- function(seq_data, var){
  seq_data %>% filter(temp %in% const_temp) %>% 
    ggplot(., aes(x=resource, y=!!sym(var), color=as.factor(temp), group = as.factor(temp))) + 
    geom_line(size = 1) + 
    labs(x = "Resource (mgC/L)", 
         y = "Spores consumed",
         color = "Temp") + 
    scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + 
    scale_y_continuous(breaks = c(0, 2000, 4000, 6000), limits = c(0, 6100)) + 
    guides(color = guide_legend(reverse=T)) + 
    scale_color_manual(values = c("#FFC107", "#21918c", "#440154"))
}

u_plot <- function(seq_data, var){
  seq_data %>% 
    filter(temp %in% const_temp) %>% 
    ggplot(.,aes(x=resource, y=!!sym(var), color=as.factor(temp))) + 
    geom_line(size = 1) + 
    scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
    scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + 
    scale_y_continuous(breaks = c(0, 0.001, 0.002), limits = c(0, 0.002)) + 
    guides(color = guide_legend(reverse=T)) + 
    labs(x = "Resource (mgC/L)", y = "Per-Spore Susceptibility", color = "Temp", title = "")
}

prev_plot <- function(seq_data, var){
  seq_data %>% filter(temp %in% const_temp) %>% 
    ggplot(.,aes(x=resource, 
                 y=!!sym(var), 
                 color=as.factor(temp))) + 
    geom_line(size=1) +
    scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
    new_scale_color() + 
    geom_point(aes(x=resource, y=prev, color=as.factor(temp)), size=3) + 
    geom_linerange(aes(ymin = conf$lower, ymax=conf$upper, color=as.factor(temp))) +
    scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
    scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + 
    scale_y_continuous(breaks = c(0, 0.5, 1.0)) + 
    guides(color = guide_legend(reverse=T)) + 
    labs(x="Resource (mgC/L)", y="Infection Prevalence", color="Temperature", title = "")
}


inside_theme <- 
  theme(axis.text = element_blank(), 
        axis.title = element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm")) #top, right, bottom, left

bottom_theme <- 
  theme(axis.text.y = element_blank(), 
        axis.title = element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm")) #top, right, bottom, left

algae_theme <- theme(axis.text = element_blank(), 
                     axis.title = element_blank(), 
                     plot.title = element_blank(), 
                     legend.position = "none")

# Figure 2: foraging results ----------------------------------------------


f0 <- seq_plot_alt(f_seq, "m1_rate") + 
  theme(axis.title = element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm"))

f1 <- seq_plot_alt(f_seq, "m2_rate") +
  bottom_theme

f2 <- seq_plot_alt(f_seq, "m3_rate") +
  bottom_theme

f3 <- seq_plot_alt(f_seq, "m4_rate") +
  bottom_theme

f4 <- seq_plot_alt(f_seq, "m5_rate") + 
  theme(axis.text.y = element_blank(), 
        axis.title = element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none", 
        axis.title.y = element_blank(), 
        panel.background = element_rect(fill = "#ADD0E9"),
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm"))

clearance_plots <- f0 | f1 | f2 | f3 | f4


ggsave(here("figures", "02_clearance_plots.png"),
       clearance_plots,
       width = 8,
       height = 1.5,
       units = "in")




a0 <- algae_plot_diff(f_seq, "m1_diff") +
  theme(axis.text.x = element_blank(), 
        axis.title = element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm"))

a1 <- algae_plot_diff(f_seq, "m2_diff") +
  inside_theme

a2 <- algae_plot_diff(f_seq, "m3_diff") +
  inside_theme

a3 <- algae_plot_diff(f_seq, "m4_diff") +
  inside_theme

a4 <- algae_plot_diff(f_seq, "m5_diff") + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none", 
        axis.title.y = element_blank(), 
        panel.background = element_rect(fill = "#ADD0E9"),
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm"))


algae_plots_diff <- a0 | a1 | a2 | a3 | a4


ggsave(here("figures", "02_algae_plots.png"),
       algae_plots_diff,
       width = 8,
       height = 1.5,
       units = "in")



# Figure 3: infection results ---------------------------------------------



s0 <- spores_plot(u_seq, "spores_consumed_m1") + 
  theme(axis.text.y = element_text(size = 6, angle = 25),
        axis.text.x = element_blank(),
        axis.title = element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm"))


s1 <- spores_plot(u_seq, "spores_consumed_m2") + 
  inside_theme

s2 <- spores_plot(u_seq, "spores_consumed_m3") + 
  inside_theme

s3 <- spores_plot(u_seq, "spores_consumed_m4") + 
  inside_theme

s4 <- spores_plot(u_seq, "spores_consumed_m5") + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none", 
        axis.title.y = element_blank(), 
        panel.background = element_rect(fill = "#ADD0E9"),
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm"))

spores_plots <- s0 | s1 | s2 | s3 | s4


ggsave(here("figures", "03_inf_spores_plots.png"),
       spores_plots,
       width = 8,
       height = 1.5,
       units = "in")



u0 <- u_plot(u_seq, "m1_susc") + 
  theme(axis.text.y = element_text(size = 6, angle = 25),
        axis.text.x = element_blank(),
        axis.title = element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm"))

u1 <- u_plot(u_seq, "m2_susc") + 
  inside_theme

u2 <- u_plot(u_seq, "m3_susc") + 
  inside_theme

u3 <- u_plot(u_seq, "m4_susc") + 
  inside_theme

u4 <- u_plot(u_seq, "m5_susc") + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none", 
        axis.title.y = element_blank(), 
        panel.background = element_rect(fill = "#ADD0E9"),
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm"))

susc_plots <- u0 | u1 | u2 | u3 | u4


ggsave(here("figures", "03_susc_plots.png"),
       susc_plots,
       width = 8,
       height = 1.5,
       units = "in")


p_theme <-
  theme(axis.text.y = element_blank(), 
        axis.title = element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm"))

p0 <- prev_plot(u_seq, "m1_I_end") + 
  theme(axis.title = element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm"))

p1 <- prev_plot(u_seq, "m2_I_end") + 
  p_theme

p2 <- prev_plot(u_seq, "m3_I_end") + 
  p_theme

p3 <- prev_plot(u_seq, "m4_I_end") + 
  p_theme

p4 <- prev_plot(u_seq, "m5_I_end") + 
  theme(axis.title=element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none", 
        panel.background = element_rect(fill = "#ADD0E9"), 
        axis.text.y = element_blank(),
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm"))


prev_plots <- p0 | p1 | p2 | p3 | p4


ggsave(here("figures", "03_prev_plots.png"),
       prev_plots,
       width = 8,
       height = 1.5,
       units = "in")




# Figure 4: estimates -----------------------------------------------------




param_ci <- 
  mod_quantiles %>% 
  ggplot(aes(x = name, 
             y = est)) + 
  geom_point() + 
  geom_linerange(aes(ymin = lower.025, 
                     ymax = upper.975)) + 
  geom_hline(yintercept = 0) + 
  labs(x = "",
       y = "",
       title = "Maximum Likelihood Estimates (95% CI)") + 
  facet_wrap(~name, scales = "free",
             nrow = 2) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(here("figures", "04_param_est.png"), width = 6, height = 4)

#handling time



h_plot <-
  h_ci_seq %>%
  ggplot(aes(x = temp)) + 
  geom_ribbon(aes(ymin = h_lower,
                  ymax = h_upper),
              fill = "grey70") + 
  geom_line(aes(y = h)) + 
  labs(y = "Handling time (days/mgC)",
       x = "Temperature (C)",
       title = "Handling time with 95% CI") + 
  theme_bw(base_size = 12)


ggsave(here("figures", "04_handling.png"), width = 6, height = 4)




# Figure 5: bootstrapped CI's ---------------------------------------------


f_res_factor <-
f_ci_seq %>%
  filter(resource %in% c(0.1, 0.5, 1)) %>%
  ggplot(aes(x = temp, group = resource)) + 
  geom_ribbon(aes(ymin = f_lower,
                  ymax = f_upper,
                  fill = as.factor(resource)),
              alpha = 0.3) + 
  geom_line(aes(y = f, color = as.factor(resource))) + 
  labs(y = "Foraging Rate (mL/day)",
       x = "Temperature (C)",
       fill = "Resource\n(mgC/L)",
       color = "",
       title = "") + 
  guides(color = "none", fill = "none") + 
  scale_color_manual(values = c("#CABD88", "#CA5216", "#65320D")) +
  scale_fill_manual(values = c("#CABD88", "#CA5216", "#65320D")) +
  theme_bw(base_size = 12)

ggsave(here("figures", "05_fora_ci_temp.png"), 
       plot = set_panel_size(p=f_res_factor, width=unit(5, "in"), height=unit(3, "in")),
       width = 7,
       height = 4,
       units = "in")


f_temp_factor <-
f_ci_seq %>%
  filter(temp %in% c(15, 20, 25)) %>%
  ggplot(aes(x = resource, group = temp)) + 
  geom_ribbon(aes(ymin = f_lower,
                  ymax = f_upper,
                  fill = as.factor(temp)),
              alpha = 0.3) + 
  geom_line(aes(y = f, color = as.factor(temp))) + 
  labs(y = "Foraging Rate (mL/day)",
       x = "Resource (mgC/L)",
       fill = "Temp (C)",
       title = "") + 
  guides(color = "none", fill = guide_legend(reverse=T)) + 
  guides(color = "none", fill = "none") + 
  scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
  scale_fill_manual(values = c("#FFC107", "#21918c", "#440154")) +
  theme(legend.position = "none") + 
  scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + 
  theme_bw(base_size = 12)


ggsave(here("figures", "05_fora_ci_rsc.png"), 
       plot = set_panel_size(p=f_temp_factor, width=unit(5, "in"), height=unit(3, "in")),
       width = 7,
       height = 4,
       units = "in")




inf_res_factor <-
u_ci_seq %>%
  filter(resource %in% c(0.1, 0.5, 1)) %>%
  ggplot(aes(x = temp, group = resource)) + 
  geom_ribbon(aes(ymin = u_lower/100000,
                  ymax = u_upper/100000,
                  fill = as.factor(resource)),
              alpha = 0.3) + 
  geom_line(aes(y = u, color = as.factor(resource))) + 
  labs(y = "Per-Parasite Susceptibility",
       x = "Temperature (C)",
       fill = "Resource\n(mgC/L)",
       color = "",
       title = "") + 
  scale_color_manual(values = c("#CABD88", "#CA5216", "#65320D")) +
  scale_fill_manual(values = c("#CABD88", "#CA5216", "#65320D")) +
  theme(legend.position = "none") + 
  guides(color = "none", fill = "none") + 
  theme_bw(base_size = 12)


ggsave(here("figures", "05_susc_ci_temp.png"), 
       plot = set_panel_size(p=inf_res_factor, width=unit(5, "in"), height=unit(3, "in")),
       width = 7,
       height = 4,
       units = "in")

inf_temp_factor <-
u_ci_seq %>%
  filter(temp %in% c(15, 20, 25)) %>%
  ggplot(aes(x = resource, group = temp)) + 
  geom_ribbon(aes(ymin = u_lower/100000,
                  ymax = u_upper/100000,
                  fill = as.factor(temp)),
              alpha = 0.3) + 
  geom_line(aes(y = u, color = as.factor(temp))) + 
  labs(y = "Per-Parasite Susceptibility",
       x = "Resource (mgC/L)",
       fill = "Temp (C)",
       color = "",
       title = "") + 
  scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
  scale_fill_manual(values = c("#FFC107", "#21918c", "#440154")) +
  theme(legend.position = "none") + 
  guides(color = "none", fill = guide_legend(reverse=T)) + 
  guides(color = "none", fill = "none") + 
  scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + 
  theme_bw(base_size = 12)


ggsave(here("figures", "05_susc_ci_rsc.png"), 
       plot = set_panel_size(p=inf_temp_factor, width=unit(5, "in"), height=unit(3, "in")),
       width = 7,
       height = 4,
       units = "in")


