#define functions for making figures
#Daniel Suh

library(here)

source(here("base","src.R"))


#read data
f_seq <- readRDS(file=here("processed_data", "fit_data","foraging_rate_fit_data.rds"))
u_seq <- readRDS(file=here("processed_data", "fit_data","infection_fit_data.rds"))

data_summ <- readRDS(here("processed_data", "foraging.rds"))

mod_boot <- readRDS(file = here("processed_data", "m2E_bootstraps.rds"))

m5_fit <- readRDS(here("processed_data", "mle", "m5_combined_fit.rds"))


# global formatting -------------------------------------------------------

theme_set(theme_bw(base_size = 12))


# some data prep ----------------------------------------------------------

data_summ %<>% dplyr::select(temp, resource, amt_init_mean, mm_mean, mm_mean_se)
f_seq %<>% left_join(., data_summ)

#f_seq %<>% mutate(amt_diff = resource-amt_rem_mean,
f_seq %<>% mutate(amt_diff = amt_init_mean-amt_rem_mean,
                  m1_diff = resource - m1_R_end/15*1000,
                  m2_diff = resource - m2_R_end/15*1000,
                  m3_diff = resource - m3_R_end/15*1000,
                  m4_diff = resource - m4_R_end/15*1000,
                  m6c_diff = resource - m6c_R_end/15*1000)


f_seq %<>% arrange(desc(temp))




# Figure 1: Data summary --------------------------------------------------
length_data <- data_summ %>%
  drop_na() %>%
  ggplot(., aes(x=resource, y=mm_mean, group = as.factor(temp), color = as.factor(temp))) + 
  geom_point(size=3) +
  geom_linerange(size = 1, 
                 aes(ymin=mm_mean-1.96*mm_mean_se,
                     ymax=mm_mean+1.96*mm_mean_se,
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

f4 <- seq_plot_alt(f_seq, "m6c_rate") + 
  theme(axis.text.y = element_blank(), 
        axis.title = element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none", 
        axis.title.y = element_blank(), 
        panel.background = element_rect(fill = "#ADD0E9"),
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm"))




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

a4 <- algae_plot_diff(f_seq, "m6c_diff") + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none", 
        axis.title.y = element_blank(), 
        panel.background = element_rect(fill = "#ADD0E9"),
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm"))



foraging <- 
  grid.arrange(grobs=lapply(list(f0, f1, f2, f3, f4,
                                 a0, a1, a2, a3, a4,
                                 fs0, fs1, fs2, fs3, fs4), 
                            set_panel_size,
                            width=unit(1.3, "in"), 
                            height=unit(1, "in")),
               ncol = 5)


ggsave(here("figures", "foraging.png"),
       foraging,
       width = 8.5,
       height = 4.5,
       units = "in")


# figure 3: infection results ---------------------------------------------



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



u0 <- u_plot(u_seq, "m1_u") + 
  theme(axis.text.y = element_text(size = 6, angle = 25),
        axis.text.x = element_blank(),
        axis.title = element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm"))

u1 <- u_plot(u_seq, "m2_u") + 
  inside_theme

u2 <- u_plot(u_seq, "m3_u") + 
  inside_theme

u3 <- u_plot(u_seq, "m4_u") + 
  inside_theme

u4 <- u_plot(u_seq, "m5_u") + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none", 
        axis.title.y = element_blank(), 
        panel.background = element_rect(fill = "#ADD0E9"),
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm"))




p0 <- prev_plot(u_seq, "m1_I_end") + 
  theme(axis.title = element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm"))

p1 <- prev_plot(u_seq, "m2_I_end") + 
  fs_theme

p2 <- prev_plot(u_seq, "m3_I_end") + 
  fs_theme

p3 <- prev_plot(u_seq, "m4_I_end") + 
  fs_theme

p4 <- prev_plot(u_seq, "m5_I_end") + 
  theme(axis.title=element_blank(), 
        plot.title = element_blank(), 
        legend.position = "none", 
        panel.background = element_rect(fill = "#ADD0E9"), 
        axis.text.y = element_blank(),
        plot.margin = unit(c(0, 0.5, 0, 0.5), "mm"))


infection <- 
  grid.arrange(grobs=lapply(list(s0, s1, s2, s3, s4,
                                 u0, u1, u2, u3, u4,
                                 p0, p1, p2, p3, p4), 
                            set_panel_size,
                            width=unit(1.3, "in"), 
                            height=unit(1, "in")),
               ncol = 5)


ggsave(here("figures", "infection.png"),
       infection,
       width = 8.5,
       height = 4.5,
       units = "in")





# Figure 4: estimates -----------------------------------------------------

#format bootstrap data

tib_length <- length(mod_boot)

mod_boot_coefs <- tibble(f = rep_len(0, tib_length),
                         u = rep_len(0, tib_length),
                         arr = rep_len(0, tib_length),
                         arr_u = rep_len(0, tib_length),
                         h = rep_len(0, tib_length),
                         w = rep_len(0, tib_length),
                         rho = rep_len(0, tib_length),
                         phi = rep_len(0, tib_length),
                         sd_est = rep_len(0, tib_length),
                         warning = NA)


for(i in 1:length(mod_boot)){
  if(mod_boot[[i]]@details[4] == 52){
    print(paste(i, "convergence_failure", sep="_"))
    mod_boot_coefs$f[i] <- NA
    mod_boot_coefs$u[i] <- NA
    mod_boot_coefs$arr[i] <- NA
    mod_boot_coefs$arr_u[i] <- NA
    mod_boot_coefs$h[i] <- NA
    mod_boot_coefs$w[i] <- NA
    mod_boot_coefs$rho[i] <- NA
    mod_boot_coefs$phi[i] <- NA
    mod_boot_coefs$sd_est[i] <- NA
    mod_boot_coefs$warning[i] <- "52"
  }
  else if(mod_boot[[i]]@details[4] == 10){
    print(paste(i, "convergence_failure", sep="_"))
    mod_boot_coefs$f[i] <- NA
    mod_boot_coefs$u[i] <- NA
    mod_boot_coefs$arr[i] <- NA
    mod_boot_coefs$arr_u[i] <- NA
    mod_boot_coefs$h[i] <- NA
    mod_boot_coefs$w[i] <- NA
    mod_boot_coefs$rho[i] <- NA
    mod_boot_coefs$phi[i] <- NA
    mod_boot_coefs$sd_est[i] <- NA
    mod_boot_coefs$warning[i] <- "10"
  }
  else if(mod_boot[[i]]@min == 0){
    print(paste(i, "zero_ll", sep="_"))
    mod_boot_coefs$f[i] <- NA
    mod_boot_coefs$u[i] <- NA
    mod_boot_coefs$arr[i] <- NA
    mod_boot_coefs$arr_u[i] <- NA
    mod_boot_coefs$h[i] <- NA
    mod_boot_coefs$w[i] <- NA
    mod_boot_coefs$rho[i] <- NA
    mod_boot_coefs$phi[i] <- NA
    mod_boot_coefs$sd_est[i] <- NA
    mod_boot_coefs$warning[i] <- "zero_ll"
  }
  else if(mod_boot[[i]]@details[4] == 1){
    print(paste(i, "max_iter_reached", sep="_"))
    boot_coef <- coef(mod_boot[[i]])
    mod_boot_coefs$f[i] <- as.numeric(boot_coef[1])
    mod_boot_coefs$u[i] <- as.numeric(boot_coef[2])
    mod_boot_coefs$arr[i] <- as.numeric(boot_coef[3])
    mod_boot_coefs$arr_u[i] <- as.numeric(boot_coef[4])
    mod_boot_coefs$h[i] <- as.numeric(boot_coef[5])
    mod_boot_coefs$w[i] <- as.numeric(boot_coef[6])
    mod_boot_coefs$rho[i] <- as.numeric(boot_coef[7])
    mod_boot_coefs$phi[i] <- as.numeric(boot_coef[8])
    mod_boot_coefs$sd_est[i] <- as.numeric(boot_coef[9])
    mod_boot_coefs$warning[i] <- "1"
  }
  else{
    boot_coef <- coef(mod_boot[[i]])
    mod_boot_coefs$f[i] <- as.numeric(boot_coef[1])
    mod_boot_coefs$u[i] <- as.numeric(boot_coef[2])
    mod_boot_coefs$arr[i] <- as.numeric(boot_coef[3])
    mod_boot_coefs$arr_u[i] <- as.numeric(boot_coef[4])
    mod_boot_coefs$h[i] <- as.numeric(boot_coef[5])
    mod_boot_coefs$w[i] <- as.numeric(boot_coef[6])
    mod_boot_coefs$rho[i] <- as.numeric(boot_coef[7])
    mod_boot_coefs$phi[i] <- as.numeric(boot_coef[8])
    mod_boot_coefs$sd_est[i] <- as.numeric(boot_coef[9])
    mod_boot_coefs$warning[i] <- "0"
  }
}

mod_quantiles <- tibble(id = c("est", "mean", "lower.025", "upper.975"),
                        f = NA,
                        u = NA,
                        arr = NA,
                        arr_u = NA,
                        h = NA,
                        w = NA,
                        rho = NA,
                        phi = NA,
                        sd_est = NA)

mod_quantiles$f[1] <- coef(m5_fit)[1]
mod_quantiles$f[2] <- mean(mod_boot_coefs$f)
mod_quantiles$f[3] <- quantile(mod_boot_coefs$f, probs=seq(0.025, 0.975, 0.95))[1]
mod_quantiles$f[4] <- quantile(mod_boot_coefs$f, probs=seq(0.025, 0.975, 0.95))[2]


mod_quantiles$u[1] <- coef(m5_fit)[2]/100000
mod_quantiles$u[2] <- mean(mod_boot_coefs$u)/100000
mod_quantiles$u[3] <- quantile(mod_boot_coefs$u, probs=seq(0.025, 0.975, 0.95))[1]/100000
mod_quantiles$u[4] <- quantile(mod_boot_coefs$u, probs=seq(0.025, 0.975, 0.95))[2]/100000


mod_quantiles$arr[1] <- coef(m5_fit)[3]
mod_quantiles$arr[2] <- mean(mod_boot_coefs$arr)
mod_quantiles$arr[3] <- quantile(mod_boot_coefs$arr, probs=seq(0.025, 0.975, 0.95))[1]
mod_quantiles$arr[4] <- quantile(mod_boot_coefs$arr, probs=seq(0.025, 0.975, 0.95))[2]

mod_quantiles$arr_u[1] <- coef(m5_fit)[4]
mod_quantiles$arr_u[2] <- mean(mod_boot_coefs$arr_u)
mod_quantiles$arr_u[3] <- quantile(mod_boot_coefs$arr_u, probs=seq(0.025, 0.975, 0.95))[1]
mod_quantiles$arr_u[4] <- quantile(mod_boot_coefs$arr_u, probs=seq(0.025, 0.975, 0.95))[2]

mod_quantiles$h[1] <- coef(m5_fit)[5]
mod_quantiles$h[2] <- mean(mod_boot_coefs$h)
mod_quantiles$h[3] <- quantile(mod_boot_coefs$h, probs=seq(0.025, 0.975, 0.95))[1]
mod_quantiles$h[4] <- quantile(mod_boot_coefs$h, probs=seq(0.025, 0.975, 0.95))[2]

mod_quantiles$w[1] <- coef(m5_fit)[6]
mod_quantiles$w[2] <- mean(mod_boot_coefs$w)
mod_quantiles$w[3] <- quantile(mod_boot_coefs$w, probs=seq(0.025, 0.975, 0.95))[1]
mod_quantiles$w[4] <- quantile(mod_boot_coefs$w, probs=seq(0.025, 0.975, 0.95))[2]


mod_quantiles$rho[1] <- coef(m5_fit)[7]
mod_quantiles$rho[2] <- mean(mod_boot_coefs$rho)
mod_quantiles$rho[3] <- quantile(mod_boot_coefs$rho, probs=seq(0.025, 0.975, 0.95))[1]
mod_quantiles$rho[4] <- quantile(mod_boot_coefs$rho, probs=seq(0.025, 0.975, 0.95))[2]

mod_quantiles$phi[1] <- coef(m5_fit)[8]
mod_quantiles$phi[2] <- mean(mod_boot_coefs$phi)
mod_quantiles$phi[3] <- quantile(mod_boot_coefs$phi, probs=seq(0.025, 0.975, 0.95))[1]
mod_quantiles$phi[4] <- quantile(mod_boot_coefs$phi, probs=seq(0.025, 0.975, 0.95))[2]

mod_quantiles$sd_est[1] <- coef(m5_fit)[9]
mod_quantiles$sd_est[2] <- mean(mod_boot_coefs$sd_est)
mod_quantiles$sd_est[3] <- quantile(mod_boot_coefs$sd_est, probs=seq(0.025, 0.975, 0.95))[1]
mod_quantiles$sd_est[4] <- quantile(mod_boot_coefs$sd_est, probs=seq(0.025, 0.975, 0.95))[2]


theme_set(theme_bw(base_size = 8))

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
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  facet_wrap(~name, scales = "free",
             nrow = 2)

ggsave(here("figures", "param_est.png"), width = 6, height = 4)

#handling time

h_seq <- tibble(temp = seq(15, 25, by = 0.1),
                h = NA,
                h_lower = NA,
                h_upper = NA)

h_est <- mod_quantiles$est[5]
w_est <- mod_quantiles$est[6]

get_handling_time <- function(temp, bound){
  tmp <- mod_boot_coefs %>% 
    mutate(h_time = h * exp(temp*w)) %>%
    summarize(min = quantile(h_time, probs = seq(0.025, 0.975, 0.95))[1],
              max = quantile(h_time, probs = seq(0.025, 0.975, 0.95))[2])
  if (bound == "upper"){
    return(tmp$max)  
  }
  else if (bound == "lower") {
    return(tmp$min) 
  }
}


h_seq %<>% mutate(h = h_est * exp(temp*w_est),
                  h_lower = mapply(get_handling_time, temp = temp, bound = "lower"),
                  h_upper = mapply(get_handling_time, temp = temp, bound = "upper"))


h_plot <-
  h_seq %>%
  ggplot(aes(x = temp)) + 
  geom_ribbon(aes(ymin = h_lower,
                  ymax = h_upper),
              fill = "grey70") + 
  geom_line(aes(y = h)) + 
  labs(y = "Handling time (days/mgC)",
       x = "Temperature (C)",
       title = "Handling time with 95% CI") + 
  theme_bw(base_size = 12)


ggsave(here("figures", "handling.png"), width = 6, height = 4)




# figure 5: bootstrapped CI's ---------------------------------------------

f_seq <- seq_data %>% 
  mutate(f = NA,
         f_lower = NA,
         f_upper = NA)


f_est <- mod_quantiles$est[1]
arr_est <- mod_quantiles$est[3]



get_fora <- function(temp, resource, bound, mm = mm){
  tmp <- mod_boot_coefs %>% 
    mutate(fora = (f*exp(arr*(1/ref_t - 1/temp))*(mm^gamma))/
             (1+f*exp(arr*(1/ref_t - 1/temp))*
                (mm^gamma)*
                h*exp(w*temp)*
                resource/1000)) %>%
    summarize(min = quantile(fora, probs = seq(0.025, 0.975, 0.95))[1],
              max = quantile(fora, probs = seq(0.025, 0.975, 0.95))[2])
  if (bound == "upper"){
    return(tmp$max)  
  }
  else if (bound == "lower") {
    return(tmp$min) 
  }
}


f_seq %<>% mutate(f = (f_est*exp(arr_est*(1/ref_t - 1/temp))*(mm^gamma))/
                    (1+f_est*exp(arr_est*(1/ref_t - 1/temp))*
                       (mm^gamma)*
                       h_est*exp(w_est*temp)*
                       resource/1000),
                  f_lower = mapply(get_fora, temp = temp, resource = resource, mm = mm, bound = "lower"),
                  f_upper = mapply(get_fora, temp = temp, resource = resource, mm = mm, bound = "upper"))

f_res_factor <-
f_seq %>%
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

ggsave(here("figures", "fora_ci_temp.png"), 
       plot = set_panel_size(p=f_res_factor, width=unit(5, "in"), height=unit(3, "in")),
       width = 7,
       height = 4,
       units = "in")


f_temp_factor <-
f_seq %>%
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
  theme_bw(base_size = 12)


ggsave(here("figures", "fora_ci_rsc.png"), 
       plot = set_panel_size(p=f_temp_factor, width=unit(5, "in"), height=unit(3, "in")),
       width = 7,
       height = 4,
       units = "in")


u_seq <- seq_data %>% 
  mutate(u = NA,
         u_lower = NA,
         u_upper = NA)



u_est <- mod_quantiles$est[2]
arr_u_est <- mod_quantiles$est[4]
rho_est <- mod_quantiles$est[7]
phi_est <- mod_quantiles$est[8]


get_susc <- function(temp, resource, bound){
  tmp <- mod_boot_coefs %>% 
    mutate(per_spore_susc = u*exp(arr_u*(1/ref_t - 1/temp))*exp(rho*resource)*exp(resource*temp*phi)) %>%
    summarize(min = quantile(per_spore_susc, probs=seq(0.025, 0.975, 0.95))[1],
              max = quantile(per_spore_susc, probs=seq(0.025, 0.975, 0.95))[2])
  if (bound == "upper"){
    return(tmp$max)  
  }
  else if (bound == "lower") {
    return(tmp$min) 
  }
}


u_seq %<>% mutate(u = u_est*exp(arr_u_est*(1/ref_t - 1/temp))*exp(rho_est*resource)*exp(resource*temp*phi_est),
                  u_lower = mapply(get_susc, temp = temp, resource = resource, bound = "lower"),
                  u_upper = mapply(get_susc, temp = temp, resource = resource, bound = "upper"),
                  u_coef = exp(arr_u_est*(1/ref_t - 1/temp))*exp(rho_est*resource)*exp(resource*temp*phi_est),
                  arr_eff = exp(arr_u_est*(1/ref_t - 1/temp)),
                  rho_eff = exp(rho_est*resource),
                  phi_eff = exp(resource*temp*phi_est))

inf_res_factor <-
u_seq %>%
  filter(resource %in% c(0.1, 0.5, 1)) %>%
  ggplot(aes(x = temp, group = resource)) + 
  geom_ribbon(aes(ymin = u_lower/100000,
                  ymax = u_upper/100000,
                  fill = as.factor(resource)),
              alpha = 0.3) + 
  geom_line(aes(y = u/100000, color = as.factor(resource))) + 
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


ggsave(here("figures", "susc_ci_temp.png"), 
       plot = set_panel_size(p=inf_res_factor, width=unit(5, "in"), height=unit(3, "in")),
       width = 7,
       height = 4,
       units = "in")

inf_temp_factor <-
u_seq %>%
  filter(temp %in% c(15, 20, 25)) %>%
  ggplot(aes(x = resource, group = temp)) + 
  geom_ribbon(aes(ymin = u_lower/100000,
                  ymax = u_upper/100000,
                  fill = as.factor(temp)),
              alpha = 0.3) + 
  geom_line(aes(y = u/100000, color = as.factor(temp))) + 
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
  theme_bw(base_size = 12)


ggsave(here("figures", "susc_ci_rsc.png"), 
       plot = set_panel_size(p=inf_temp_factor, width=unit(5, "in"), height=unit(3, "in")),
       width = 7,
       height = 4,
       units = "in")


