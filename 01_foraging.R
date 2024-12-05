#Build foraging rate model
#Daniel Suh

library(here)

source(here("base","src.R"))

#read data
data <- readRDS(here("processed_data","foraging_raw.rds"))
data_summ <- readRDS(here("processed_data", "foraging.rds"))



# Define local variables --------------------------------------------------


treatment_IDs <- unique(data$treatment_ID)
resource_IDs <- unique(data$resource)
temp_IDs <- unique(data$temp)

vol <- 15



# Define model fn's -------------------------------------------------------


# model 1A - size-dependent model -----------------------------------------



m1_num_sol <- function(t, x, params){
  R <- x[1]
  with(as.list(params),{
    dR <- -(length^gamma)*f*R
    res <- c(dR)
    list(res)}
  )
}

m1_sim <- function(R, time, f, length, gamma){
  times <- seq(0, time, by=0.01)
  xstart <- c(R=R)
  params <- c(f=f,
              length=length,
              gamma=gamma)
  output <- as.data.frame(lsoda(y=xstart, times, m1_num_sol, params))
  
  r_end <- slice_max(output, time)[,2]
  
  return(r_end)
}



m1_ll <- function(f, sd_est){
  m1_end <- as.data.frame(mapply(m1_sim, R=data$amt_init*data$vol/1000, time=data$time/60/24, f=f/data$vol,
                                 length=data$mm, gamma=2))
  colnames(m1_end) <- "m1_end"
  data$end <- m1_end$m1_end
  
  data %<>% mutate(resid = sqrt(end) - sqrt(amt_rem*vol/1000),
                   sqrt_model_end = sqrt(end),
                   sqrt_data_end = sqrt(amt_rem*vol/1000))
  nll_sum <- 0
  for(i in 1:length(treatment_IDs)){
    treatment_data <- data %>% filter(treatment_ID == treatment_IDs[i])
    nll <- dnorm(treatment_data$sqrt_data_end,
                 mean = mean(treatment_data$sqrt_model_end),
                 sd = sd_est,
                 log = T)
    nll_sum <- -sum(nll) + nll_sum
  }
  
  nll_sum
}



# model 1B - temperature-only model ---------------------------------------

m2_num_sol <- function(t, x, params){
  R <- x[1]
  with(as.list(params),{
    dR <- -(length^gamma)*f*exp(arr_t*(1/ref_t - 1/temp))*R
    res <- c(dR)
    list(res)}
  )
}

m2_sim <- function(R, time, f, length, gamma, ref_t, arr_t, temp){
  times <- seq(0, time, by=0.01)
  xstart <- c(R=R)
  params <- c(f=f,
              length=length,
              gamma=gamma,
              arr_t=arr_t,
              ref_t=ref_t,
              temp=temp)
  output <- as.data.frame(lsoda(y=xstart, times, m2_num_sol, params))
  
  r_end <- slice_max(output, time)[,2]
  
  return(r_end)
}

m2_ll <- function(f, arr_t, sd_est){
  
  m2_end <- as.data.frame(mapply(m2_sim, R=data$amt_init*data$vol/1000, time=data$time/60/24, f=f/data$vol,
                                 length=data$mm, gamma=2, arr_t=arr_t, ref_t=15, temp=data$temp))
  colnames(m2_end) <- "m2_end"
  data$end <- m2_end$m2_end
  
  data %<>% mutate(resid = sqrt(end) - sqrt(amt_rem*vol/1000),
                   sqrt_model_end = sqrt(end),
                   sqrt_data_end = sqrt(amt_rem*vol/1000))
  
  
  nll_sum <- 0
  for(i in 1:length(treatment_IDs)){
    treatment_data <- data %>% filter(treatment_ID == treatment_IDs[i])
    nll <- dnorm(treatment_data$sqrt_data_end,
                 mean = mean(treatment_data$sqrt_model_end),
                 sd = sd_est,
                 log = T)
    nll_sum <- -sum(nll) + nll_sum
  }
  
  nll_sum
}

# model 1C - resource-only model ------------------------------------------

m3_num_sol <- function(t, x, params){
  R <- x[1]
  with(as.list(params),{
    dR <- -((f*length^gamma*R)/(1+f*(length^gamma)*h*R))
    res <- c(dR)
    list(res)}
  )
}

m3_sim <- function(R, time, f, length, gamma, h){
  times <- seq(0, time, by=0.01)
  xstart <- c(R=R)
  params <- c(f=f,
              length=length,
              gamma=gamma,
              h=h)
  output <- as.data.frame(lsoda(y=xstart, times, m3_num_sol, params))
  
  r_end <- slice_max(output, time)[,2]

  return(r_end)
}

m3_ll <- function(f, h, sd_est){
  
  output <- as.data.frame(mapply(m3_sim, R=data$amt_init*data$vol/1000, time=data$time/60/24, f=f/data$vol,
                                 length=data$mm, gamma=2, h=h))
  colnames(output) <- "end"
  data$end <- output$end
  
  data %<>% mutate(resid = sqrt(end) - sqrt(amt_rem*vol/1000),
                   sqrt_model_end = sqrt(end),
                   sqrt_data_end = sqrt(amt_rem*vol/1000))
  
  
  nll_sum <- 0
  for(i in 1:length(treatment_IDs)){
    treatment_data <- data %>% filter(treatment_ID == treatment_IDs[i])
    nll <- dnorm(treatment_data$sqrt_data_end,
                 mean = mean(treatment_data$sqrt_model_end),
                 sd = sd_est,
                 log = T)
    nll_sum <- -sum(nll) + nll_sum
  }
  
  nll_sum
}


# model 1D - additive model -----------------------------------------------

m4_num_sol <- function(t, x, params){
  R <- x[1]
  with(as.list(params),{
    dR <- -((f*exp(arr_t*(1/ref_t - 1/temp))*(length^gamma)*R)/(1+f*exp(arr_t*(1/ref_t - 1/temp))*(length^gamma)*h*R))
    res <- c(dR)
    list(res)}
  )
}

m4_sim <- function(R, time, f, length, gamma, arr_t, ref_t, temp, h){
  times <- seq(0, time, by=0.01)
  xstart <- c(R=R)
  params <- c(f=f,
              length=length,
              gamma=gamma,
              arr_t=arr_t,
              ref_t=ref_t,
              temp=temp,
              h=h)
  output <- as.data.frame(lsoda(y=xstart, times, m4_num_sol, params))
  
  r_end <- slice_max(output, time)[,2]
  
  return(r_end)
}

m4_ll <- function(f, arr_t, h, sd_est){
  
  m4_end <- as.data.frame(mapply(m4_sim, R=data$amt_init*data$vol/1000, time=data$time/60/24, f=f/data$vol,
                                 length=data$mm, gamma=2, arr_t=arr_t, ref_t=15, temp=data$temp, h=h))
  colnames(m4_end) <- "m4_end"
  data$end <- m4_end$m4_end

  data %<>% mutate(resid = sqrt(end) - sqrt(amt_rem*vol/1000),
                   sqrt_model_end = sqrt(end),
                   sqrt_data_end = sqrt(amt_rem*vol/1000))
  
  
  nll_sum <- 0
  for(i in 1:length(treatment_IDs)){
    treatment_data <- data %>% filter(treatment_ID == treatment_IDs[i])
    nll <- dnorm(treatment_data$sqrt_data_end,
                 mean = mean(treatment_data$sqrt_model_end),
                 sd = sd_est,
                 log = T)
    nll_sum <- -sum(nll) + nll_sum
  }
  
  nll_sum
}

# model 1E - interactive model --------------------------------------------

m5_ll <- function(f, arr_t, h, w, sd_est){
  
  
  output <- as.data.frame(mapply(m4_sim, 
                                 R=data$amt_init*data$vol/1000, 
                                 time=data$time/60/24, 
                                 f=f/data$vol,
                                 length=data$mm, 
                                 gamma=2, 
                                 arr_t=arr_t, 
                                 ref_t=15, 
                                 temp=data$temp, 
                                 h=h*exp(w*data$temp)))
  colnames(output) <- "endpoint"
  data$end <- output$endpoint
  #data %<>% mutate(resid = log(end) - log(amt_rem*vol/1000))
  data %<>% mutate(resid = sqrt(end) - sqrt(amt_rem*vol/1000))
  data %<>% mutate(sqrt_model_end = sqrt(end),
                   sqrt_data_end = sqrt(amt_rem*vol/1000),
                   log_model_end = log(end),
                   log_data_end = log(amt_rem*vol/1000))
  
  
  nll_sum <- 0
  for(i in 1:length(treatment_IDs)){
    treatment_data <- data %>% filter(treatment_ID == treatment_IDs[i])
    nll <- dnorm(treatment_data$sqrt_data_end,
                 mean = mean(treatment_data$sqrt_model_end),
                 sd = sd_est,
                 log = T)
    nll_sum <- -sum(nll) + nll_sum
  }
  
  nll_sum
  
}


