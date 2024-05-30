library(here)
source(here("model_fitting", "combined", "final_prep.R")) #data, global variables and functions

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
