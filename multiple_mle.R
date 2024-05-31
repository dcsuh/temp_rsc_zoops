library(here)
source(here("model_fitting", "combined", "final_prep.R")) #data, global variables and functions

library(patchwork)

#model 5 results with different starting conditions.
#est_1 uses starting conditions from profile m5a profile
#m5a_profile <- readRDS(file = here("mle", "final", "m5a_profile_fit.rds"))
#est_1.1 uses 1.1x starting conditions
#est_0.5 uses 0.5x starting conditions

est_1 <- readRDS(file = here("model_fitting", "combined", "m5_1_fit.rds"))
est_1.1 <- readRDS(file = here("model_fitting", "combined", "m5_1.1_fit.rds"))
est_0.5 <- readRDS(file = here("model_fitting", "combined", "m5_0.5_fit.rds"))

coef(est_1)
coef(est_1.1)
coef(est_0.5)

getAIC(est_1)
getAIC(est_1.1)
getAIC(est_0.5)



temp_range <- seq(15,25,by=0.1)
resource_range <- seq(0.1,1,by=0.01)

seq_data <- tibble(temp=c(),
                   resourece=c())

for (i in 1:length(temp_range)){
  tmp <- tibble(temp = temp_range[i],
                resource = seq(0.1,1,by=0.01))
  seq_data %<>% bind_rows(., tmp)
}


lm(mm ~ resource + temp, data = data)

length_coef <- coef(lm(mm ~ resource + temp, data = data)) #lengths from foraging rate experiment
#length_ceof <- coef(lm(mm ~ resource + as.numeric(temp_id), data = lengths)) #5 day lengths from life table


interpolate_length <- function(resource, temp){
  length <- length_coef[1] + resource*length_coef[2] + temp*length_coef[3]
}

seq_data %<>% mutate(length = mapply(interpolate_length, 
                                     resource = resource,
                                     temp = temp))


fora_summ <- data %>% 
  mutate(ID = treatment_ID) %>%
  group_by(ID) %>% 
  summarize(end_r_mean = mean(amt_rem)*fora_vol/1000,
            end_r_se = (sd(amt_rem)/sqrt(n()))*fora_vol/1000)

prevalence <- readRDS(here("processed_data", "prevalence.rds"))
prevalence %<>% filter(temp %in% c(15, 20, 25)) %>% filter(species=="D") %>% mutate(ID = paste(temp, resource, sep="_")) %>% dplyr::select(-c(temp, resource, species))


viz_summ <- left_join(fora_summ, prevalence)
viz_summ %<>% left_join(., data_summ)
viz_summ %<>% left_join(., length_summ)


seq_data %<>% left_join(., viz_summ, by=c("temp", "resource"))


r_end_mean <- 
  data %>% group_by(temp, resource) %>% 
  summarize(amt_rem_mean = mean(amt_rem),
            amt_rem_se = sd(amt_rem)/sqrt(n()),
            amt_consumed_mean = mean(amt_consumed),
            amt_consumed_se = sd(amt_consumed)/sqrt(n())) %>%
  ungroup()

seq_data %<>% left_join(., r_end_mean)


coef(est_1)

m5b_u_f_fit <- est_1

m5b_coef <- coef(m5b_u_f_fit)
m5b_u_prime <- as.numeric(m5b_coef[1])
m5b_f <- as.numeric(m5b_coef[2])
m5b_arr <- as.numeric(m5b_coef[3])
m5b_arr_u <- as.numeric(m5b_coef[4])
m5b_h <- as.numeric(m5b_coef[5])
m5b_w <- as.numeric(m5b_coef[6])
m5b_rho_15 <- as.numeric(m5b_coef[7])
m5b_rho_20 <- as.numeric(m5b_coef[8])
m5b_rho_25 <- as.numeric(m5b_coef[9])

seq_data %<>% mutate(est_1_rate = (m5b_f*exp(m5b_arr*(1/ref_t - 1/temp))*(length^gamma))/
                       (1+m5b_f*exp(m5b_arr*(1/ref_t - 1/temp))*(length^gamma)*m5b_h*exp(m5b_w*temp)*resource/1000),
                     est_1_u = case_when(temp == 15 ~ m5b_u_prime*exp(m5b_arr_u*(1/ref_t - 1/temp))*exp(m5b_rho_15*resource),
                                       temp == 20 ~ m5b_u_prime*exp(m5b_arr_u*(1/ref_t - 1/temp))*exp(m5b_rho_20*resource),
                                       temp == 25 ~ m5b_u_prime*exp(m5b_arr_u*(1/ref_t - 1/temp))*exp(m5b_rho_25*resource)),
                     est_1_arr = m5b_u_prime*exp(m5b_arr_u*(1/ref_t - 1/temp)),
                     est_1_rho = case_when(temp == 15 ~ m5b_u_prime*exp(m5b_rho_15*resource),
                                         temp == 20 ~ m5b_u_prime*exp(m5b_rho_20*resource),
                                         temp == 25 ~ m5b_u_prime*exp(m5b_rho_25*resource)))

seq_data %<>% mutate(m5_1_end = mapply(m4_sim, 
                                       R=resource*life_vol/1000, 
                                       time=1, 
                                       f=m5b_f/life_vol, 
                                       u = case_when(temp == 15 ~ m5b_u_prime*exp(m5b_rho_15*resource),
                                                     temp == 20 ~ m5b_u_prime*exp(m5b_rho_20*resource),
                                                     temp == 25 ~ m5b_u_prime*exp(m5b_rho_25*resource)), 
                                       length=length, 
                                       gamma=gamma, 
                                       Z=200*life_vol,
                                       ref_t = ref_t,
                                       arr_t_f = m5b_arr,
                                       arr_t_u = m5b_arr_u,
                                       h = m5b_h*exp(m5b_w*temp),
                                       temp = temp))

coef(est_1.1)

m5b_u_f_fit <- est_1.1

m5b_coef <- coef(m5b_u_f_fit)
m5b_u_prime <- as.numeric(m5b_coef[1])
m5b_f <- as.numeric(m5b_coef[2])
m5b_arr <- as.numeric(m5b_coef[3])
m5b_arr_u <- as.numeric(m5b_coef[4])
m5b_h <- as.numeric(m5b_coef[5])
m5b_w <- as.numeric(m5b_coef[6])
m5b_rho_15 <- as.numeric(m5b_coef[7])
m5b_rho_20 <- as.numeric(m5b_coef[8])
m5b_rho_25 <- as.numeric(m5b_coef[9])

seq_data %<>% mutate(est_1.1_rate = (m5b_f*exp(m5b_arr*(1/ref_t - 1/temp))*(length^gamma))/
                       (1+m5b_f*exp(m5b_arr*(1/ref_t - 1/temp))*(length^gamma)*m5b_h*exp(m5b_w*temp)*resource/1000),
                     est_1.1_u = case_when(temp == 15 ~ m5b_u_prime*exp(m5b_arr_u*(1/ref_t - 1/temp))*exp(m5b_rho_15*resource),
                                       temp == 20 ~ m5b_u_prime*exp(m5b_arr_u*(1/ref_t - 1/temp))*exp(m5b_rho_20*resource),
                                       temp == 25 ~ m5b_u_prime*exp(m5b_arr_u*(1/ref_t - 1/temp))*exp(m5b_rho_25*resource)),
                     est_1.1_arr = m5b_u_prime*exp(m5b_arr_u*(1/ref_t - 1/temp)),
                     est_1.1_rho = case_when(temp == 15 ~ m5b_u_prime*exp(m5b_rho_15*resource),
                                           temp == 20 ~ m5b_u_prime*exp(m5b_rho_20*resource),
                                           temp == 25 ~ m5b_u_prime*exp(m5b_rho_25*resource)))

seq_data %<>% mutate(m5_1.1_end = mapply(m4_sim, 
                                         R=resource*life_vol/1000, 
                                         time=1, 
                                         f=m5b_f/life_vol, 
                                         u = case_when(temp == 15 ~ m5b_u_prime*exp(m5b_rho_15*resource),
                                                       temp == 20 ~ m5b_u_prime*exp(m5b_rho_20*resource),
                                                       temp == 25 ~ m5b_u_prime*exp(m5b_rho_25*resource)), 
                                         length=length, 
                                         gamma=gamma, 
                                         Z=200*life_vol,
                                         ref_t = ref_t,
                                         arr_t_f = m5b_arr,
                                         arr_t_u = m5b_arr_u,
                                         h = m5b_h*exp(m5b_w*temp),
                                         temp = temp))

coef(est_0.5)

m5b_u_f_fit <- est_0.5

m5b_coef <- coef(m5b_u_f_fit)
m5b_u_prime <- as.numeric(m5b_coef[1])
m5b_f <- as.numeric(m5b_coef[2])
m5b_arr <- as.numeric(m5b_coef[3])
m5b_arr_u <- as.numeric(m5b_coef[4])
m5b_h <- as.numeric(m5b_coef[5])
m5b_w <- as.numeric(m5b_coef[6])
m5b_rho_15 <- as.numeric(m5b_coef[7])
m5b_rho_20 <- as.numeric(m5b_coef[8])
m5b_rho_25 <- as.numeric(m5b_coef[9])

seq_data %<>% mutate(est_0.5_rate = (m5b_f*exp(m5b_arr*(1/ref_t - 1/temp))*(length^gamma))/
                       (1+m5b_f*exp(m5b_arr*(1/ref_t - 1/temp))*(length^gamma)*m5b_h*exp(m5b_w*temp)*resource/1000),
                     est_0.5_u = case_when(temp == 15 ~ m5b_u_prime*exp(m5b_arr_u*(1/ref_t - 1/temp))*exp(m5b_rho_15*resource),
                                       temp == 20 ~ m5b_u_prime*exp(m5b_arr_u*(1/ref_t - 1/temp))*exp(m5b_rho_20*resource),
                                       temp == 25 ~ m5b_u_prime*exp(m5b_arr_u*(1/ref_t - 1/temp))*exp(m5b_rho_25*resource)),
                     est_0.5_arr = m5b_u_prime*exp(m5b_arr_u*(1/ref_t - 1/temp)),
                     est_0.5_rho = case_when(temp == 15 ~ m5b_u_prime*exp(m5b_rho_15*resource),
                                           temp == 20 ~ m5b_u_prime*exp(m5b_rho_20*resource),
                                           temp == 25 ~ m5b_u_prime*exp(m5b_rho_25*resource)))

seq_data %<>% mutate(m5_0.5_end = mapply(m4_sim, 
                                         R=resource*life_vol/1000, 
                                         time=1, 
                                         f=m5b_f/life_vol, 
                                         u = case_when(temp == 15 ~ m5b_u_prime*exp(m5b_rho_15*resource),
                                                       temp == 20 ~ m5b_u_prime*exp(m5b_rho_20*resource),
                                                       temp == 25 ~ m5b_u_prime*exp(m5b_rho_25*resource)), 
                                         length=length, 
                                         gamma=gamma, 
                                         Z=200*life_vol,
                                         ref_t = ref_t,
                                         arr_t_f = m5b_arr,
                                         arr_t_u = m5b_arr_u,
                                         h = m5b_h*exp(m5b_w*temp),
                                         temp = temp))
#f
seq_data %>% ggplot(., aes(x=est_1_rate, y=est_1.1_rate)) + geom_point() + geom_abline()
seq_data %>% ggplot(., aes(x=est_1_rate, y=est_0.5_rate)) + geom_point() + geom_abline()
seq_data %>% ggplot(., aes(x=est_1.1_rate, y=est_0.5_rate)) + geom_point() + geom_abline()

#arr
seq_data %>% ggplot(., aes(x=est_1_arr, y=est_1.1_arr)) + geom_point() + geom_abline()
seq_data %>% ggplot(., aes(x=est_1_arr, y=est_0.5_arr)) + geom_point() + geom_abline()
seq_data %>% ggplot(., aes(x=est_1.1_arr, y=est_0.5_arr)) + geom_point() + geom_abline()

#rho
seq_data %>% ggplot(., aes(x=est_1_rho, y=est_1.1_rho)) + geom_point() + geom_abline()
seq_data %>% ggplot(., aes(x=est_1_rho, y=est_0.5_rho)) + geom_point() + geom_abline()
seq_data %>% ggplot(., aes(x=est_1.1_rho, y=est_0.5_rho)) + geom_point() + geom_abline()

#arr X rho
seq_data %>% ggplot(., aes(x=est_1_arr*est_1_rho, y=est_1.1_arr*est_1.1_rho)) + geom_point() + geom_abline()
seq_data %>% ggplot(., aes(x=est_1_arr*est_1_rho, y=est_0.5_arr*est_0.5_rho)) + geom_point() + geom_abline()
seq_data %>% ggplot(., aes(x=est_1.1_arr*est_1.1_rho, y=est_0.5_arr*est_0.5_rho)) + geom_point() + geom_abline()

#prev
seq_data %>% ggplot(., aes(x=m5_1_end, y=m5_1.1_end)) + geom_point() + geom_abline()
seq_data %>% ggplot(., aes(x=m5_1_end, y=m5_0.5_end)) + geom_point() + geom_abline()
seq_data %>% ggplot(., aes(x=m5_1.1_end, y=m5_0.5_end)) + geom_point() + geom_abline()

#u
seq_data %>% ggplot(., aes(x=est_1_u, y=est_1.1_u, color = temp)) + geom_point() + geom_abline() + scale_color_viridis_c()
seq_data %>% ggplot(., aes(x=est_1_u, y=est_0.5_u, color = temp)) + geom_point() + geom_abline() + scale_color_viridis_c()
seq_data %>% ggplot(., aes(x=est_1.1_u, y=est_0.5_u, color = temp)) + geom_point() + geom_abline() + scale_color_viridis_c()



est_m4 <- readRDS(file = here("mle", "final", "m4_fit.rds"))
est_m4_1 <- readRDS(file = here("model_fitting", "m4_1_fit.rds"))
est_m4_1.1 <- readRDS(file = here("model_fitting", "m4_1.1_fit.rds"))
est_m4_0.5 <- readRDS(file = here("model_fitting", "m4_0.5_fit.rds"))

coef(est_m4_1)
coef(est_m4_1.1)
coef(est_m4_0.5)

getAIC(est_m4)
getAIC(est_m4_1)
getAIC(est_m4_1.1)
getAIC(est_m4_0.5)


m4_coef <- coef(est_m4)
m4_f <- as.numeric(m4_coef[1])
m4_u <- as.numeric(m4_coef[2])
m4_arr <- as.numeric(m4_coef[3])
m4_arr_u <- as.numeric(m4_coef[4])
m4_h <- as.numeric(m4_coef[5])
m4_w <- as.numeric(m4_coef[6])
m4_rho <- as.numeric(m4_coef[7])

m1_sim_Z <- function(R, time, f, u, length, gamma, Z, ref_t, arr_t_f, h, temp){
  times <- seq(0, time, by=0.01)
  xstart <- c(R=R,
              S=1,
              I=0,
              Z=Z)
  params <- c(f=f,
              u=0,
              length = length,
              gamma = gamma,
              ref_t = ref_t,
              arr_t_f = arr_t_f,
              h = h,
              temp = temp)
  output <- as.data.frame(lsoda(y=xstart, times, m1_num_sol, params))
  
  end_data <- slice_min(output, time)[,5] - slice_max(output, time)[,5]
  
  
  return(end_data)
}



seq_data %<>% mutate(m4_end = mapply(m4_sim, 
                                       R=resource*life_vol/1000, 
                                       time=1, 
                                       f=m4_f/life_vol, 
                                       u = m4_u*exp(m4_rho*resource), 
                                       length=length, 
                                       gamma=gamma, 
                                       Z=200*life_vol,
                                       ref_t = ref_t,
                                       arr_t_f = m4_arr,
                                       arr_t_u = m4_arr_u,
                                       h = m4_h*exp(m4_w*temp),
                                       temp = temp))

seq_data %<>% mutate(m4_rate = (m4_f*exp(m4_arr*(1/ref_t - 1/temp))*(length^gamma))/
                       (1+m4_f*exp(m4_arr*(1/ref_t - 1/temp))*(length^gamma)*m4_h*exp(m4_w*temp)*resource/1000),
                     m4_u = m4_u*exp(m4_arr_u*(1/ref_t - 1/temp))*exp(m4_rho*resource),
                     m4_res = m4_h*exp(m4_w*temp),
                     m4_rho = exp(m4_rho*resource),
                     m4_arr = exp(m4_arr_u*(1/ref_t - 1/temp)))

seq_data %<>% mutate(m4_sp = mapply(m1_sim_Z, 
                                   R=resource*life_vol/1000, 
                                   time=1, 
                                   f=m4_f/life_vol, 
                                   u = 0, 
                                   length=length, 
                                   gamma=gamma, 
                                   Z=200*life_vol,
                                   ref_t = ref_t,
                                   arr_t_f = m4_arr,
                                   h = m4_h*exp(m4_w*temp),
                                   temp = temp))

seq_data %<>% mutate(m4_R = mapply(m4_sim, 
                                     R=resource*life_vol/1000, 
                                     time=1, 
                                     f=m4_f/life_vol, 
                                     u = m4_u*exp(m4_rho*resource), 
                                     length=length, 
                                     gamma=gamma, 
                                     Z=0,
                                     ref_t = ref_t,
                                     arr_t_f = m4_arr,
                                     arr_t_u = m4_arr_u,
                                     h = m4_h*exp(m4_w*temp),
                                     temp = temp))


m4_coef <- coef(est_m4_1)
m4_f <- as.numeric(m4_coef[1])
m4_u <- as.numeric(m4_coef[2])
m4_arr <- as.numeric(m4_coef[3])
m4_arr_u <- as.numeric(m4_coef[4])
m4_h <- as.numeric(m4_coef[5])
m4_w <- as.numeric(m4_coef[6])
m4_rho <- as.numeric(m4_coef[7])


seq_data %<>% mutate(m4_1_end = mapply(m4_sim, 
                                       R=resource*life_vol/1000, 
                                       time=1, 
                                       f=m4_f/life_vol, 
                                       u = m4_u*exp(m4_rho*resource), 
                                       length=length, 
                                       gamma=gamma, 
                                       Z=200*life_vol,
                                       ref_t = ref_t,
                                       arr_t_f = m4_arr,
                                       arr_t_u = m4_arr_u,
                                       h = m4_h*exp(m4_w*temp),
                                       temp = temp))

seq_data %<>% mutate(m4_1_rate = (m4_f*exp(m4_arr*(1/ref_t - 1/temp))*(length^gamma))/
                       (1+m4_f*exp(m4_arr*(1/ref_t - 1/temp))*(length^gamma)*m4_h*exp(m4_w*temp)*resource/1000),
                     m4_1_u = m4_u*exp(m4_arr_u*(1/ref_t - 1/temp))*exp(m4_rho*resource),
                     m4_1_res = m4_h*exp(m4_w*temp),
                     m4_1_rho = exp(m4_rho*resource),
                     m4_1_arr = exp(m4_arr_u*(1/ref_t - 1/temp)))

seq_data %<>% mutate(m4_1_sp = mapply(m1_sim_Z, 
                                    R=resource*life_vol/1000, 
                                    time=1, 
                                    f=m4_f/life_vol, 
                                    u = 0, 
                                    length=length, 
                                    gamma=gamma, 
                                    Z=200*life_vol,
                                    ref_t = ref_t,
                                    arr_t_f = m4_arr,
                                    h = m4_h*exp(m4_w*temp),
                                    temp = temp))

seq_data %<>% mutate(m4_1_R = mapply(m4_sim, 
                                   R=resource*life_vol/1000, 
                                   time=1, 
                                   f=m4_f/life_vol, 
                                   u = m4_u*exp(m4_rho*resource), 
                                   length=length, 
                                   gamma=gamma, 
                                   Z=0,
                                   ref_t = ref_t,
                                   arr_t_f = m4_arr,
                                   arr_t_u = m4_arr_u,
                                   h = m4_h*exp(m4_w*temp),
                                   temp = temp))

m4_coef <- coef(est_m4_1.1)
m4_f <- as.numeric(m4_coef[1])
m4_u <- as.numeric(m4_coef[2])
m4_arr <- as.numeric(m4_coef[3])
m4_arr_u <- as.numeric(m4_coef[4])
m4_h <- as.numeric(m4_coef[5])
m4_w <- as.numeric(m4_coef[6])
m4_rho <- as.numeric(m4_coef[7])


seq_data %<>% mutate(m4_1.1_end = mapply(m4_sim, 
                                       R=resource*life_vol/1000, 
                                       time=1, 
                                       f=m4_f/life_vol, 
                                       u = m4_u*exp(m4_rho*resource), 
                                       length=length, 
                                       gamma=gamma, 
                                       Z=200*life_vol,
                                       ref_t = ref_t,
                                       arr_t_f = m4_arr,
                                       arr_t_u = m4_arr_u,
                                       h = m4_h*exp(m4_w*temp),
                                       temp = temp))

seq_data %<>% mutate(m4_1.1_rate = (m4_f*exp(m4_arr*(1/ref_t - 1/temp))*(length^gamma))/
                       (1+m4_f*exp(m4_arr*(1/ref_t - 1/temp))*(length^gamma)*m4_h*exp(m4_w*temp)*resource/1000),
                     m4_1.1_u = m4_u*exp(m4_arr_u*(1/ref_t - 1/temp))*exp(m4_rho*resource),
                     m4_1.1_res = m4_h*exp(m4_w*temp),
                     m4_1.1_rho = exp(m4_rho*resource),
                     m4_1.1_arr = exp(m4_arr_u*(1/ref_t - 1/temp)))


seq_data %<>% mutate(m4_1.1_sp = mapply(m1_sim_Z, 
                                    R=resource*life_vol/1000, 
                                    time=1, 
                                    f=m4_f/life_vol, 
                                    u = 0, 
                                    length=length, 
                                    gamma=gamma, 
                                    Z=200*life_vol,
                                    ref_t = ref_t,
                                    arr_t_f = m4_arr,
                                    h = m4_h*exp(m4_w*temp),
                                    temp = temp))

seq_data %<>% mutate(m4_1.1_R = mapply(m4_sim, 
                                   R=resource*life_vol/1000, 
                                   time=1, 
                                   f=m4_f/life_vol, 
                                   u = m4_u*exp(m4_rho*resource), 
                                   length=length, 
                                   gamma=gamma, 
                                   Z=0,
                                   ref_t = ref_t,
                                   arr_t_f = m4_arr,
                                   arr_t_u = m4_arr_u,
                                   h = m4_h*exp(m4_w*temp),
                                   temp = temp))


m4_coef <- coef(est_m4_0.5)
m4_f <- as.numeric(m4_coef[1])
m4_u <- as.numeric(m4_coef[2])
m4_arr <- as.numeric(m4_coef[3])
m4_arr_u <- as.numeric(m4_coef[4])
m4_h <- as.numeric(m4_coef[5])
m4_w <- as.numeric(m4_coef[6])
m4_rho <- as.numeric(m4_coef[7])


seq_data %<>% mutate(m4_0.5_end = mapply(m4_sim, 
                                       R=resource*life_vol/1000, 
                                       time=1, 
                                       f=m4_f/life_vol, 
                                       u = m4_u*exp(m4_rho*resource), 
                                       length=length, 
                                       gamma=gamma, 
                                       Z=200*life_vol,
                                       ref_t = ref_t,
                                       arr_t_f = m4_arr,
                                       arr_t_u = m4_arr_u,
                                       h = m4_h*exp(m4_w*temp),
                                       temp = temp))

seq_data %<>% mutate(m4_0.5_rate = (m4_f*exp(m4_arr*(1/ref_t - 1/temp))*(length^gamma))/
                       (1+m4_f*exp(m4_arr*(1/ref_t - 1/temp))*(length^gamma)*m4_h*exp(m4_w*temp)*resource/1000),
                     m4_0.5_u = m4_u*exp(m4_arr_u*(1/ref_t - 1/temp))*exp(m4_rho*resource),
                     m4_0.5_res = m4_h*exp(m4_w*temp),
                     m4_0.5_rho = exp(m4_rho*resource),
                     m4_0.5_arr = exp(m4_arr_u*(1/ref_t - 1/temp)))

seq_data %<>% mutate(m4_0.5_sp = mapply(m1_sim_Z, 
                                    R=resource*life_vol/1000, 
                                    time=1, 
                                    f=m4_f/life_vol, 
                                    u = 0, 
                                    length=length, 
                                    gamma=gamma, 
                                    Z=200*life_vol,
                                    ref_t = ref_t,
                                    arr_t_f = m4_arr,
                                    h = m4_h*exp(m4_w*temp),
                                    temp = temp))

seq_data %<>% mutate(m4_0.5_R = mapply(m4_sim, 
                                   R=resource*life_vol/1000, 
                                   time=1, 
                                   f=m4_f/life_vol, 
                                   u = m4_u*exp(m4_rho*resource), 
                                   length=length, 
                                   gamma=gamma, 
                                   Z=0,
                                   ref_t = ref_t,
                                   arr_t_f = m4_arr,
                                   arr_t_u = m4_arr_u,
                                   h = m4_h*exp(m4_w*temp),
                                   temp = temp))

#f
#seq_data %>% ggplot(., aes(x=m4_rate, y=m4_1_rate)) + geom_point() + geom_abline()
f_1_1.1 <- seq_data %>% ggplot(., aes(x=m4_1_rate, y=m4_1.1_rate)) + geom_point() + geom_abline()
f_1_0.5 <- seq_data %>% ggplot(., aes(x=m4_1_rate, y=m4_0.5_rate)) + geom_point() + geom_abline()
#seq_data %>% ggplot(., aes(x=m4_1.1_rate, y=m4_0.5_rate)) + geom_point() + geom_abline()

#prev
prev_1_1.1 <- seq_data %>% ggplot(., aes(x=m4_1_end, y=m4_1.1_end, color = temp)) + geom_point() + geom_abline()+ scale_color_viridis_c()
prev_1_0.5 <- seq_data %>% ggplot(., aes(x=m4_1_end, y=m4_0.5_end, color = temp)) + geom_point() + geom_abline()+ scale_color_viridis_c()
#seq_data %>% ggplot(., aes(x=m4_1.1_end, y=m4_0.5_end, color = temp)) + geom_point() + geom_abline()+ scale_color_viridis_c()

#u
seq_data %>% ggplot(., aes(x=m4_1_u, y=m4_1.1_u, color = temp)) + geom_point() + geom_abline() + scale_color_viridis_c()
seq_data %>% ggplot(., aes(x=m4_1_u, y=m4_0.5_u, color = temp)) + geom_point() + geom_abline() + scale_color_viridis_c()
seq_data %>% ggplot(., aes(x=m4_1.1_u, y=m4_0.5_u, color = temp)) + geom_point() + geom_abline() + scale_color_viridis_c()

#temp on resource interaction
#seq_data %>% ggplot(., aes(x=m4_res, y=m4_1_res, color = temp)) + geom_point() + scale_color_viridis_c()
seq_data %>% ggplot(., aes(x=m4_1_res, y=m4_1.1_res, color = temp)) + geom_point() + scale_color_viridis_c() + geom_abline()
seq_data %>% ggplot(., aes(x=m4_1_res, y=m4_0.5_res, color = temp)) + geom_point() + scale_color_viridis_c() + geom_abline()
#seq_data %>% ggplot(., aes(x=m4_1.1_res, y=m4_0.5_res, color = temp)) + geom_point() + scale_color_viridis_c()

#temp effect on u
seq_data %>% ggplot(., aes(x=m4_1_arr, y=m4_1.1_arr, color = temp)) + geom_point() + scale_color_viridis_c() + geom_abline()
seq_data %>% ggplot(., aes(x=m4_1_arr, y=m4_0.5_arr, color = temp)) + geom_point() + scale_color_viridis_c() + geom_abline()

#resource effect on u
seq_data %>% ggplot(., aes(x=m4_1_rho, y=m4_1.1_rho, color = temp)) + geom_point() + scale_color_viridis_c() + geom_abline()
seq_data %>% ggplot(., aes(x=m4_1_rho, y=m4_0.5_rho, color = temp)) + geom_point() + scale_color_viridis_c() + geom_abline()


#spores consumed
seq_data %>% ggplot(., aes(x=m4_1_sp, y=m4_1.1_sp)) + geom_point() + geom_abline()
seq_data %>% ggplot(., aes(x=m4_1_sp, y=m4_0.5_sp)) + geom_point() + geom_abline()
seq_data %>% ggplot(., aes(x=m4_1.1_sp, y=m4_0.5_sp)) + geom_point() + geom_abline()


seq_plot <- function(seq_data, rate){
  seq_data %>% filter(temp %in% c(25,20,15)) %>%
    ggplot(.) +
    geom_line(aes(x=resource, 
                  y=!!sym(rate),
                  group=as.factor(temp),
                  color=as.factor(temp)),
              size = 1) +
    scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
    guides(color = guide_legend(reverse=T)) + 
    labs(x="Resource (mgC/L)", y="Foraging rate (mL/day)", color="Temperature", shape="Temperature", title = "") 
}

seq_plot(seq_data, "m4_rate")
seq_plot(seq_data, "m4_1_rate")
seq_plot(seq_data, "m4_1.1_rate")
seq_plot(seq_data, "m4_0.5_rate")



spores_plot <- function(seq_data, var){
  seq_data %>% filter(temp %in% const_temp) %>% 
    ggplot(., aes(x=resource, y=!!sym(var), color=as.factor(temp), group = as.factor(temp))) + 
    geom_line(size = 1) + 
    labs(x = "Resource (mgC/L)", 
         y = "Spores consumed",
         color = "Temp") + 
    scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + 
    guides(color = guide_legend(reverse=T)) + 
    #scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) + 
    scale_color_manual(values = c("#FFC107", "#21918c", "#440154"))
}

spores_plot(seq_data, "m4_sp")
spores_plot(seq_data, "m4_1_sp")
spores_plot(seq_data, "m4_1.1_sp")
spores_plot(seq_data, "m4_0.5_sp")


seq_data %<>%
  mutate(amt_diff = amt_init_mean-amt_rem_mean,
         m4_diff = resource - m4_R/50*1000,
         m4_1_diff = resource - m4_1_R/50*1000,
         m4_1.1_diff = resource - m4_1.1_R/50*1000,
         m4_0.5_diff = resource - m4_0.5_R/50*1000)


algae_plot_diff <- function(seq_data, var){
  seq_data %>% filter(temp %in% const_temp) %>%
    ggplot(., aes(x=resource, y=!!sym(var), color=as.factor(temp), group = as.factor(temp))) +
    geom_line(size=1) +
    labs(x = "Resource (mgC/L)",
         y = "Resource difference (mgC/L)",
         color = "Temp") +
    scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + 
    guides(color = guide_legend(reverse=T)) + 
    #scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D"))
    scale_color_manual(values = c("#FFC107", "#21918c", "#440154"))
}

algae_plot_diff(seq_data, "m4_diff")
algae_plot_diff(seq_data, "m4_1_diff")
algae_plot_diff(seq_data, "m4_1.1_diff")
algae_plot_diff(seq_data, "m4_0.5_diff")


prev_plot <- function(seq_data, var){
  seq_data %>% filter(temp %in% const_temp) %>% 
    ggplot(.,aes(x=resource, 
                 y=!!sym(var), 
                 color=as.factor(temp))) + 
    geom_line(size=1) +
    #scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
    scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
    new_scale_color() + 
    geom_point(aes(x=resource, y=prev, color=as.factor(temp)), size=3) + 
    geom_linerange(aes(ymin = conf$lower, ymax=conf$upper, color=as.factor(temp))) +
    #scale_color_manual(values = c("#619CFF", "#00BA38", "#F8766D")) +
    scale_color_manual(values = c("#FFC107", "#21918c", "#440154")) +
    scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + 
    scale_y_continuous(breaks = c(0, 0.5, 1.0)) + 
    guides(color = guide_legend(reverse=T)) + 
    labs(x="Resource (mgC/L)", y="Infection Prevalence", color="Temperature", title = "")
}

prev_plot(seq_data, "m4_end")
prev_plot(seq_data, "m4_1_end")
prev_plot(seq_data, "m4_1.1_end")
prev_plot(seq_data, "m4_0.5_end")




