#generate foraging rate model results
#Daniel Suh

library(here)



if(file.exists(here("processed_data", "seq_data", "foraging_rate_fit_data.rds")) == TRUE) {
  
  message("Simulated data already exist.")
  
} else {
  


source(here("base","src.R"))
source(here("01_foraging.R"))



m1 <- readRDS(here("processed_data", "mle", "m1_f_fit.rds"))
m2 <- readRDS(here("processed_data", "mle", "m2_f_fit.rds"))
m3 <- readRDS(here("processed_data", "mle", "m3_f_fit.rds"))
m4 <- readRDS(here("processed_data", "mle", "m4_f_fit.rds"))
m5 <- readRDS(here("processed_data", "mle", "m5_f_fit.rds"))



# cleaning ----------------------------------------------------------------

#change units
data_summ %<>% mutate(rate_len_mean = rate_len_mean*60*24,
                      rate_len_mean_se = rate_len_mean_se*60*24)



# define local variables --------------------------------------------------



gamma <- 2
ref_t <- 15
vol <- 15
mean_time <- mean(data$time/60/24)

exposure <- 200
life_vol <- 50




# make skeleton data frame ------------------------------------------------

temp_range <- seq(15,25,by=0.1)
resource_range <- seq(min(data$amt_init), max(data$amt_init),by=0.01)
resource_range <- append(resource_range, values = c(0.1, 0.5, 1))

seq_data <- tibble(temp=c(),
                   resourece=c())

for (i in 1:length(temp_range)){
  tmp <- tibble(temp = temp_range[i],
                resource = resource_range)
  seq_data %<>% bind_rows(., tmp)
}


seq_data %<>% mutate(temp_factor = as.factor(case_when(temp==25 ~ 25,
                                                       temp==20 ~ 20,
                                                       temp==15 ~ 15)))

#interpolate length
lm(mm ~ resource + temp, data = data)

length_coef <- coef(lm(mm ~ resource + temp, data = data))


interpolate_length <- function(resource, temp){
  length <- length_coef[1] + resource*length_coef[2] + temp*length_coef[3]
}

seq_data %<>% mutate(length = mapply(interpolate_length, 
                                     resource = resource,
                                     temp = temp))

#interpolate initial resource concentration
lm(amt_init ~ resource, data = data)

proj_res_coef <- coef(lm(amt_init ~ resource, data = data))

interpolate_resources <- function(resource){
  output <- proj_res_coef[1] + resource*proj_res_coef[2]
}

seq_data %<>% mutate(proj_res = mapply(interpolate_resources,
                                       resource = resource))

data_summ_se <- data_summ %>% 
  dplyr::select(temp, resource, rate_len_mean, rate_len_mean_se, amt_init_mean)

seq_data %<>% left_join(., data_summ_se)

r_end_mean <- 
  data %>% 
  group_by(temp, resource) %>%
  summarize(amt_rem_mean = mean(amt_rem),
            amt_rem_se = sd(amt_rem)/sqrt(n()))

seq_data %<>% left_join(., r_end_mean)



# get foraging rate -------------------------------------------------------


m1_coef <- coef(m1)
m1_f <- as.numeric(m1_coef[1])

m2_coef <- coef(m2)
m2_f <- as.numeric(m2_coef[1])
m2_arr <- as.numeric(m2_coef[2])

m3_coef <- coef(m3)
m3_f <- as.numeric(m3_coef[1])
m3_h <- as.numeric(m3_coef[2])

m4_coef <- m4@coef
m4_f <- as.numeric(m4_coef[1])
m4_arr <- as.numeric(m4_coef[2])
m4_h <- as.numeric(m4_coef[3])

m5_coef <- m5@coef
m5_f <- as.numeric(m5_coef[1])
m5_arr <- as.numeric(m5_coef[2])
m5_h <- as.numeric(m5_coef[3])
m5_w <- as.numeric(m5_coef[4])


get_rate <- function(f, arr, temp, h, resource, length){
  ((f*length^gamma*exp(arr*(1/ref_t - 1/temp)))/
     (1+f*length^gamma*exp(arr*(1/ref_t - 1/temp))*h*resource))
}


seq_data %<>% mutate(m1_rate = length^gamma*m1_f,
                     m2_rate = (length^gamma)*m2_f*exp(m2_arr*(1/ref_t - 1/temp)),
                     m3_rate=((m3_f*length^gamma)/(1+m3_f*length^gamma*m3_h*resource/1000)),
                     m4_rate = mapply(get_rate, 
                                      f=m4_f, 
                                      arr=m4_arr, 
                                      temp=seq_data$temp, 
                                      h=m4_h, 
                                      resource=seq_data$resource/1000, 
                                      length = length),
                     m5_rate = mapply(get_rate, 
                                       f=m5_f, 
                                       arr=m5_arr, 
                                       h=m5_h*exp(m5_w*seq_data$temp), 
                                       temp=seq_data$temp, 
                                       resource=seq_data$resource/1000, 
                                       length = length))



# solve for final resource quantity ---------------------------------------

seq_data %<>% mutate(m1_R_end = mapply(m1_sim,
                                       R = resource*vol/1000,
                                       time = mean_time,
                                       f=m1_f/vol,
                                       length = length,
                                       gamma = 2),
                     m2_R_end = mapply(m2_sim,
                                       R = resource*vol/1000,
                                       time = mean_time,
                                       f=m2_f/vol,
                                       length = length,
                                       gamma = 2,
                                       ref_t = 15,
                                       arr_t = m2_arr,
                                       temp = temp),
                     m3_R_end = mapply(m3_sim,
                                       R = resource*vol/1000,
                                       time = mean_time,
                                       f=m3_f/vol,
                                       length = length,
                                       gamma = 2,
                                       h = m3_h),
                     m4_R_end = mapply(m4_sim,
                                       R = resource*vol/1000,
                                       time = mean_time,
                                       f=m4_f/vol,
                                       length = length,
                                       gamma = 2,
                                       ref_t = 15,
                                       arr_t = m4_arr,
                                       temp = temp,
                                       h = m4_h),
                     m5_R_end = mapply(m4_sim,
                                        R = resource*vol/1000,
                                        time = mean_time,
                                        f=m5_f/vol,
                                        length = length,
                                        gamma = 2,
                                        ref_t = 15,
                                        arr_t = m5_arr,
                                        temp = temp,
                                        h = m5_h*exp(m5_w*temp)),
                     m1_consumed = amt_init_mean*15/1000 - m1_R_end,
                     m2_consumed = amt_init_mean*15/1000 - m2_R_end,
                     m3_consumed = amt_init_mean*15/1000 - m3_R_end,
                     m4_consumed = amt_init_mean*15/1000 - m4_R_end,
                     m5_consumed = amt_init_mean*15/1000 - m5_R_end)



# simulate spore exposure -------------------------------------------------


m1_num_sol_z <- function(t, x, params){
  R <- x[1]
  Z <- x[2]
  with(as.list(params),{
    dR <- -(length^gamma)*f*R
    dZ <- -(length^gamma)*f*Z
    res <- c(dR, dZ)
    list(res)}
  )
}

m1_sim_z <- function(R, time, f, length, gamma, Z){
  times <- seq(0, time, by=0.01)
  xstart <- c(R=R,
              Z=Z)
  params <- c(f=f,
              length=length,
              gamma=gamma)
  output <- as.data.frame(lsoda(y=xstart, times, m1_num_sol_z, params))
  
  end <- slice_min(output, time)[,3] - slice_max(output, time)[,3]
  
  
  return(end)
}

m2_num_sol_z <- function(t, x, params){
  R <- x[1]
  Z <- x[2]
  with(as.list(params),{
    dR <- -(length^gamma)*f*exp(arr_t*(1/ref_t - 1/temp))*R
    dZ <- -(length^gamma)*f*exp(arr_t*(1/ref_t - 1/temp))*Z
    res <- c(dR, dZ)
    list(res)}
  )
}

m2_sim_z <- function(R, time, f, length, gamma, ref_t, arr_t, temp, Z){
  times <- seq(0, time, by=0.01)
  xstart <- c(R=R,
              Z=Z)
  params <- c(f=f,
              length=length,
              gamma=gamma,
              arr_t=arr_t,
              ref_t=ref_t,
              temp=temp)
  output <- as.data.frame(lsoda(y=xstart, times, m2_num_sol_z, params))
  
  end <- slice_min(output, time)[,3] - slice_max(output, time)[,3]
  
  
  return(end)
}

m3_num_sol_z <- function(t, x, params){
  R <- x[1]
  Z <- x[2]
  with(as.list(params),{
    dR <- -((f*length^gamma*R)/(1+f*(length^gamma)*h*R))
    dZ <- -((f*length^gamma*Z)/(1+f*(length^gamma)*h*R))
    res <- c(dR, dZ)
    list(res)}
  )
}

m3_sim_z <- function(R, time, f, length, gamma, h, Z){
  times <- seq(0, time, by=0.01)
  xstart <- c(R=R,
              Z=Z)
  params <- c(f=f,
              length=length,
              gamma=gamma,
              h=h)
  output <- as.data.frame(lsoda(y=xstart, times, m3_num_sol_z, params))
  
  end <- slice_min(output, time)[,3] - slice_max(output, time)[,3]
  
  
  return(end)
}

m4_num_sol_z <- function(t, x, params){
  R <- x[1]
  Z <- x[2]
  with(as.list(params),{
    dR <- -((f*exp(arr_t*(1/ref_t - 1/temp))*(length^gamma)*R)/
              (1+f*exp(arr_t*(1/ref_t - 1/temp))*(length^gamma)*h*R))
    dZ <- -((f*exp(arr_t*(1/ref_t - 1/temp))*(length^gamma))/
              (1+f*exp(arr_t*(1/ref_t - 1/temp))*(length^gamma)*h*R))*Z
    res <- c(dR, dZ)
    list(res)}
  )
}

m4_sim_z <- function(R, time, f, length, gamma, arr_t, ref_t, temp, h, Z){
  times <- seq(0, time, by=0.01)
  xstart <- c(R=R,
              Z=Z)
  params <- c(f=f,
              length=length,
              gamma=gamma,
              arr_t=arr_t,
              ref_t=ref_t,
              temp=temp,
              h=h)
  output <- as.data.frame(lsoda(y=xstart, times, m4_num_sol_z, params))
  
  end <- slice_min(output, time)[,3] - slice_max(output, time)[,3]
  
  
  return(end)
}


seq_data %<>% mutate(m1_Z_end = mapply(m1_sim_z,
                                       R = resource*life_vol/1000,
                                       time = 1,
                                       f=m1_f/life_vol,
                                       length = length,
                                       gamma = 2,
                                       Z = exposure*life_vol),
                     m2_Z_end = mapply(m2_sim_z,
                                       R = resource*life_vol/1000,
                                       time = 1,
                                       f=m2_f/life_vol,
                                       length = length,
                                       gamma = 2,
                                       ref_t = 15,
                                       arr_t = m2_arr,
                                       temp = temp,
                                       Z = exposure*life_vol),
                     m3_Z_end = mapply(m3_sim_z,
                                       R = resource*life_vol/1000,
                                       time = 1,
                                       f=m3_f/life_vol,
                                       length = length,
                                       gamma = 2,
                                       h = m3_h,
                                       Z = exposure*life_vol),
                     m4_Z_end = mapply(m4_sim_z,
                                       R = resource*life_vol/1000,
                                       time = 1,
                                       f=m4_f/life_vol,
                                       length = length,
                                       gamma = 2,
                                       ref_t = 15,
                                       arr_t = m4_arr,
                                       temp = temp,
                                       h = m4_h,
                                       Z = exposure*life_vol),
                     m5_Z_end = mapply(m4_sim_z,
                                       R = resource*life_vol/1000,
                                       time = 1,
                                       f=m5_f/life_vol,
                                       length = length,
                                       gamma = 2,
                                       ref_t = 15,
                                       arr_t = m5_arr,
                                       temp = temp,
                                       h = m5_h*exp(m5_w*temp),
                                       Z = exposure*life_vol))



# save data ---------------------------------------------------------------

if(dir.exists(here("processed_data", "seq_data")) == FALSE) {
  message("Welcome! Let's make some room for simulated model data.")
  dir.create(here("processed_data", "seq_data")) 
} else {
  message("/processed_data/seq_data exists! Proceeeding to save.")
}

saveRDS(seq_data, file=here("processed_data", "seq_data", "foraging_rate_fit_data.rds"))

}

