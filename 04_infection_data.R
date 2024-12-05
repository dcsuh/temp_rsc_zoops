#Generate transmission model results
#Daniel Suh

library(here)


if(file.exists(here("processed_data", "seq_data", "infection_fit_data.rds")) == TRUE) {
  
  message("Simulated data already exist.")
  
} else {
  



source(here("base","src.R"))
source(here("03_infection.R"))



m1 <- readRDS(here("processed_data", "mle", "m1_combined_fit.rds"))
m2 <- readRDS(here("processed_data", "mle", "m2_combined_fit.rds"))
m3 <- readRDS(here("processed_data", "mle", "m3_combined_fit.rds"))
m4 <- readRDS(here("processed_data", "mle", "m4_combined_fit.rds"))
m5 <- readRDS(here("processed_data", "mle", "m5_combined_fit.rds"))


temp_range <- seq(15,25,by=0.1)
resource_range <- seq(0.1,1,by=0.01)

seq_data <- tibble(temp=c(),
                   resource=c())

for (i in 1:length(temp_range)){
  tmp <- tibble(temp = temp_range[i],
                resource = seq(0.1,1,by=0.01))
  seq_data %<>% bind_rows(., tmp)
}


lm(mm ~ resource + temp, data = data)

length_coef <- coef(lm(mm ~ resource + temp, data = data)) #lengths from foraging rate experiment
#length_coef <- coef(lm(mm ~ resource + as.numeric(temp_id), data = lengths)) #5 day lengths from life table


interpolate_length <- function(resource, temp){
  length <- length_coef[1] + resource*length_coef[2] + temp*length_coef[3]
}

seq_data %<>% mutate(length = mapply(interpolate_length, 
                                     resource = resource,
                                     temp = temp))

prevalence <- readRDS(here("processed_data", "prevalence.rds"))
prevalence %<>%  
  mutate(resource = as.numeric(as.character(resource))) %>%
  mutate(ID = paste(temp, resource, sep="_"))

seq_data %<>% left_join(., prevalence, by=c("temp", "resource"))


# solve for terminal infection prevalence ---------------------------------



m1_coef <- coef(m1)
m1_u <- as.numeric(m1_coef[1])
m1_f <- as.numeric(m1_coef[2])
m1_arr <- as.numeric(m1_coef[3])
m1_h <- as.numeric(m1_coef[4])
m1_w <- as.numeric(m1_coef[5])

m2_coef <- coef(m2)
m2_f <- as.numeric(m2_coef[1])
m2_u <- as.numeric(m2_coef[2])
m2_arr <- as.numeric(m2_coef[3])
m2_arr_u <- as.numeric(m2_coef[4])
m2_h <- as.numeric(m2_coef[5])
m2_w <- as.numeric(m2_coef[6])

m3_coef <- coef(m3)
m3_f <- as.numeric(m3_coef[1])
m3_u <- as.numeric(m3_coef[2])/10000
m3_arr <- as.numeric(m3_coef[3])
m3_h <- as.numeric(m3_coef[4])
m3_w <- as.numeric(m3_coef[5])
m3_rho <- as.numeric(m3_coef[6])

m4_coef <- coef(m4)
m4_f <- as.numeric(m4_coef[1])
m4_u <- as.numeric(m4_coef[2])/10000
m4_arr <- as.numeric(m4_coef[3])
m4_arr_u <- as.numeric(m4_coef[4])
m4_h <- as.numeric(m4_coef[5])
m4_w <- as.numeric(m4_coef[6])
m4_rho <- as.numeric(m4_coef[7])

m5_coef <- coef(m5)
m5_f <- as.numeric(m5_coef[1])
m5_u <- as.numeric(m5_coef[2])/100000
m5_arr <- as.numeric(m5_coef[3])
m5_arr_u <- as.numeric(m5_coef[4])
m5_h <- as.numeric(m5_coef[5])
m5_w <- as.numeric(m5_coef[6])
m5_rho <- as.numeric(m5_coef[7])
m5_phi <- as.numeric(m5_coef[8])


seq_data %<>% mutate(m1_I_end = mapply(m1_sim, 
                                       R=resource*life_vol/1000, 
                                       time=1, 
                                       f=m1_f/life_vol, 
                                       u = m1_u, 
                                       length=length, 
                                       gamma=gamma, 
                                       Z=200*life_vol,
                                       ref_t = ref_t,
                                       arr_t_f = m1_arr,
                                       h = m1_h*exp(m1_w*temp),
                                       temp = temp),
                     m2_I_end = mapply(m2_sim, 
                                       R=resource*life_vol/1000,
                                       time=1,
                                       f=m2_f/life_vol,
                                       u = m2_u,
                                       length=length,
                                       gamma=gamma,
                                       Z=200*life_vol,
                                       ref_t = ref_t,
                                       arr_t_f = m2_arr,
                                       arr_t_u = m2_arr_u,
                                       h = m2_h*exp(m2_w*temp),
                                       temp = temp),
                     m3_I_end = mapply(m3_sim, 
                                       R=resource*life_vol/1000,
                                       time=1,
                                       f=m3_f/life_vol,
                                       u = m3_u*exp(m3_rho*resource),
                                       length=length,
                                       gamma=gamma,
                                       Z=200*life_vol,
                                       ref_t = ref_t,
                                       arr_t_f = m3_arr,
                                       h = m3_h*exp(m3_w*temp),
                                       temp = temp),
                     m4_I_end = mapply(m4_sim, 
                                       R=resource*life_vol/1000, 
                                       time=1, 
                                       f=m4_f/life_vol, 
                                       u = m4_u, 
                                       length=length, 
                                       gamma=gamma, 
                                       Z=200*life_vol,
                                       ref_t = ref_t,
                                       arr_t_f = m4_arr,
                                       arr_t_u = m4_arr_u,
                                       h = m4_h,
                                       temp = temp,
                                       resource = resource,
                                       w = m4_w,
                                       rho = m4_rho),
                     m5_I_end = mapply(m5_sim, 
                                       R=resource*life_vol/1000, 
                                       time=1, 
                                       f=m5_f/life_vol, 
                                       u = m5_u, 
                                       length=length, 
                                       gamma=gamma, 
                                       Z=200*life_vol,
                                       ref_t = ref_t,
                                       arr_t_f = m5_arr,
                                       arr_t_u = m5_arr_u,
                                       h = m5_h,
                                       temp = temp,
                                       resource = resource,
                                       w = m5_w,
                                       rho = m5_rho,
                                       phi = m5_phi))



# estimate rate and susceptibility ----------------------------------------



seq_data %<>% mutate(m1_rate = (m1_f*exp(m1_arr*(1/ref_t - 1/temp))*(length^gamma))/
                       (1+m1_f*exp(m1_arr*(1/ref_t - 1/temp))*(length^gamma)*m1_h*exp(m1_w*temp)*resource/1000),
                     m2_rate = (m2_f*exp(m2_arr*(1/ref_t - 1/temp))*(length^gamma))/
                       (1+m2_f*exp(m2_arr*(1/ref_t - 1/temp))*(length^gamma)*m2_h*exp(m2_w*temp)*resource/1000),
                     m3_rate = (m3_f*exp(m3_arr*(1/ref_t - 1/temp))*(length^gamma))/
                       (1+m3_f*exp(m2_arr*(1/ref_t - 1/temp))*(length^gamma)*m3_h*exp(m3_w*temp)*resource/1000),
                     m4_rate = (m4_f*exp(m4_arr*(1/ref_t - 1/temp))*(length^gamma))/
                       (1+m4_f*exp(m4_arr*(1/ref_t - 1/temp))*(length^gamma)*m4_h*exp(m4_w*temp)*resource/1000),
                     m5_rate = (m5_f*exp(m5_arr*(1/ref_t - 1/temp))*(length^gamma))/
                       (1+m5_f*exp(m5_arr*(1/ref_t - 1/temp))*
                          (length^gamma)*
                          m5_h*exp(m5_w*temp)*
                          resource/1000),
                     m1_susc = m1_u,
                     m2_susc = m2_u*exp(m2_arr_u*(1/ref_t - 1/temp)),
                     m3_susc = m3_u*exp(m3_rho*resource),
                     m4_susc = m4_u*exp(m4_arr_u*(1/ref_t - 1/temp))*exp(m4_rho*resource),
                     m5_susc = m5_u*
                       exp(m5_arr_u*(1/ref_t - 1/temp))*
                       exp(m5_rho*resource)*
                       exp(resource*temp*m5_phi))



# spores consumed ---------------------------------------------------------

best_f <- readRDS(here("processed_data", "mle", "m5_f_fit.rds")) # full model with exponential effect of temp on handling time

best_f_coef <- coef(best_f)

f_est <- best_f_coef[1] #ml/day
arr_f_est <- best_f_coef[3] #Celsius
h_est <- best_f_coef[5] #mg dry weight carbon/day
#h is a huge number because it is milligrams carbon/day but daphnia are eating micrograms per day
w_est <- best_f_coef[6] #per-degree Celsius (i.e. 1/Celsius)

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



seq_data %<>% mutate(spores_consumed_m1 = mapply(m1_sim_Z, 
                                                 R=resource*life_vol/1000, 
                                                 time=1, 
                                                 f=m1_f/life_vol, 
                                                 u = 0, 
                                                 length=length, 
                                                 gamma=gamma, 
                                                 Z=200*life_vol,
                                                 ref_t = ref_t,
                                                 arr_t_f = m1_arr,
                                                 h = m1_h*exp(m1_w*temp),
                                                 temp = temp),
                     spores_consumed_m2 = mapply(m1_sim_Z, 
                                                 R=resource*life_vol/1000, 
                                                 time=1, 
                                                 f=m2_f/life_vol, 
                                                 u = 0, 
                                                 length=length, 
                                                 gamma=gamma, 
                                                 Z=200*life_vol,
                                                 ref_t = ref_t,
                                                 arr_t_f = m2_arr,
                                                 h = m2_h*exp(m2_w*temp),
                                                 temp = temp),
                     spores_consumed_m3 = mapply(m1_sim_Z, 
                                                 R=resource*life_vol/1000, 
                                                 time=1, 
                                                 f=m3_f/life_vol, 
                                                 u = 0, 
                                                 length=length, 
                                                 gamma=gamma, 
                                                 Z=200*life_vol,
                                                 ref_t = ref_t,
                                                 arr_t_f = m3_arr,
                                                 h = m3_h*exp(m3_w*temp),
                                                 temp = temp),
                     spores_consumed_m4 = mapply(m1_sim_Z, 
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
                                                 temp = temp),
                     spores_consumed_m5 = mapply(m1_sim_Z, 
                                                 R=resource*life_vol/1000, 
                                                 time=1, 
                                                 f=m5_f/life_vol, 
                                                 u = 0, 
                                                 length=length, 
                                                 gamma=gamma, 
                                                 Z=200*life_vol,
                                                 ref_t = ref_t,
                                                 arr_t_f = m5_arr,
                                                 h = m5_h*exp(m5_w*temp),
                                                 temp = temp))



# save data ---------------------------------------------------------------

if(dir.exists(here("processed_data", "seq_data")) == FALSE) {
  message("Welcome! Let's make some room for simulated model data.")
  dir.create(here("processed_data", "seq_data")) 
} else {
  message("/processed_data/seq_data exists! Proceeeding to save.")
}


saveRDS(seq_data, file=here("processed_data", "seq_data", "infection_fit_data.rds"))

}
