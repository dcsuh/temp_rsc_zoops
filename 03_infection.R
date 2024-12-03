#Build transmission model
#Daniel Suh

library(here)

source(here("base","src.R"))


#read infection data
mort <-  read_csv(here("raw_data/infection.csv")) #mortality data
lengths <- read_csv(here("raw_data", "day5_length.csv")) #get average lengths

#read foraging data
data <- readRDS(here("processed_data","foraging_raw.rds")) #foraging rate data
data_summ <- readRDS(here("processed_data", "foraging.rds")) #foraging rate data summary



data %<>% mutate(treatment_ID = paste(temp, resource, sep = "_"))

data %<>% mutate(resource_tot = resource*vol) #resource_tot is the total amount of resources in our tube. resource concentration * total volume of tube




# define local variables --------------------------------------------------


treatment_IDs <- unique(data$treatment_ID)
resource_IDs <- unique(data$resource)
temp_IDs <- unique(data$temp)

fora_vol <- 15 #mL
life_vol <- 50 #mL
spore_conc <- 200 #spores/mL
gamma <- 2 
ref_t <- 15 #celsius




# cleaning ----------------------------------------------------------------

#infection data
dataset <- 
  mort %>%
  mutate(final_date = ifelse(is.na(mortality_day), end_data_date, mortality_day))


dataset %<>% filter(species == "D") %>% #only daphnia for this analysis
  filter(temp %in% const_temp) #only constant temp for this analysis

dataset %<>% mutate(birthdate = ifelse(species == "daphnia", "4/5/22", "4/6/22"), 
                    lifespan = as.numeric(mdy(final_date) - mdy(birthdate)),
                    temp = as.numeric(temp),
                    treatment = paste(temp, resource, sep = "_")) %>% 
  filter(is.na(male)) %>% #remove males
  filter(is.na(missing)) %>% #remove missing
  filter(!is.na(inf)) #remove NAs for inf for estimating beta. sometimes they died too young to tell


dataset %<>% mutate(inf_status = inf, 
                    spore_exposure = 200, #spores/mL
                    uninf = 1-inf_status,
                    time = 1, #duration of exposure in days
                    trt = paste(temp, resource, species, sep = "_"))


#lengths from infection assay
lengths %<>% mutate(mm = raw_meas*17.86/1000)
#at default magnification (5.6x), 1 unit is equal to 17.86 micron

length_summ <- lengths %>% 
  group_by(temp_id, resource) %>% 
  summarize(life_mm = mean(mm),
            var = var(mm),
            sd = sd(mm),
            se = sd(mm)/sqrt(n())) %>%
  ungroup() %>%
  dplyr::select(temp_id, resource, life_mm)

length_summ %<>% 
  add_row(temp_id=as.numeric(25), 
          resource=as.numeric(0.5), 
          #midpoint between 1.0 and 0.1
          life_mm=1.16) %>%
  rename(temp = temp_id)

dataset %<>% left_join(., length_summ)

length_summ %<>% mutate(ID = paste(temp, resource, sep="_")) %>% dplyr::select(-c(temp, resource))

#interpolate lengths from foraging assay
mean_length_summ <- data_summ %>% dplyr::select(temp, resource, mm_mean)
data %<>% left_join(., mean_length_summ)
data %<>% mutate(mm = if_else(!is.na(mm), mm, mm_mean))

mean_length_summ %<>% rename(fora_mm = "mm_mean")
dataset %<>% left_join(., mean_length_summ)




# define model fn's -------------------------------------------------------



# model 2A - independent --------------------------------------------------


m1_num_sol <- function(t, x, params){
  R <- x[1]
  S <- x[2]
  I <- x[3]
  Z <- x[4]
  
  with(as.list(params),{
    dR <- -((f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma))/
              (1+f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma)*h*R))*R*(S+I)
    dS <- -u*((f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma))/
                (1+f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma)*h*R))*Z*S
    dI <- u*((f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma))/
               (1+f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma)*h*R))*Z*S
    dZ <- -((f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma))/
              (1+f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma)*h*R))*Z*(S+I)
    res <- c(dR, dS, dI, dZ)
    list(res)}
  )
}

m1_sim <- function(R, time, f, u, length, gamma, Z, ref_t, arr_t_f, h, temp){
  times <- seq(0, time, by=0.01)
  xstart <- c(R=R,
              S=1,
              I=0,
              Z=Z)
  params <- c(f=f,
              u=u,
              length = length,
              gamma = gamma,
              Z = Z,
              ref_t = ref_t,
              arr_t_f = arr_t_f,
              h = h,
              temp = temp)
  output <- as.data.frame(lsoda(y=xstart, times, m1_num_sol, params))
  
  if(Z == 0) {
    end_data <- slice_max(output, time)[,2]
  }
  else {
    end_data <- slice_max(output, time)[,4]  
  }
  
  return(end_data)
}

m1_ll <- function(u, f, arr_t_f, h, w, sd_est){
  
  R_end <- as.data.frame(mapply(m1_sim, 
                                R=data$amt_init*fora_vol/1000, 
                                time=data$time/60/24, 
                                f=f/fora_vol,
                                u = 0, 
                                length=data$mm, 
                                gamma=gamma, 
                                Z=0,
                                ref_t=ref_t,
                                arr_t_f=arr_t_f,
                                h=h*exp(w*data$temp),
                                temp=data$temp))
  
  colnames(R_end) <- "R_end"
  data$end <- R_end$R_end
  
  data %<>% mutate(resid = sqrt(end) - sqrt(amt_rem*fora_vol/1000))
  data %<>% mutate(sqrt_model_end = sqrt(end),
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
  
  exp_data <- dataset %>% filter(exposed==TRUE)
  
  I_end <- as.data.frame(mapply(m1_sim, 
                                R=exp_data$resource*life_vol/1000, 
                                time=exp_data$time, 
                                f=f/life_vol, 
                                u=u, 
                                length = exp_data$life_mm, 
                                gamma=gamma, 
                                Z=spore_conc*life_vol,
                                ref_t=ref_t,
                                arr_t_f=arr_t_f,
                                h=h*exp(w*exp_data$temp),
                                temp=exp_data$temp))
  colnames(I_end) <- "I_end"
  
  exp_data %<>% cbind(I_end)
  
  exp_data %<>% 
    mutate(ll = mapply(inf_out, inf_status = inf_status, I_end=I_end))
  nll_sum <- nll_sum + -sum(exp_data$ll)
  return(nll_sum)
}



# model 2B - Temperature-only ---------------------------------------------


m2_sim <-  function(R, time, f, u, length, gamma, Z, ref_t, arr_t_f, arr_t_u, h, temp){
  times <- seq(0, time, by=0.01)
  xstart <- c(R=R,
              S=1,
              I=0,
              Z=Z)
  params <- c(f=f,
              u=u*exp(arr_t_u*(1/ref_t - 1/temp)),
              length = length,
              gamma = gamma,
              Z = Z,
              ref_t = ref_t,
              arr_t_f = arr_t_f,
              h = h,
              temp = temp)
  output <- as.data.frame(lsoda(y=xstart, times, m1_num_sol, params))
  
  if(Z == 0) {
    end_data <- slice_max(output, time)[,2]
  }
  else {
    end_data <- slice_max(output, time)[,4]  
  }
  
  return(end_data)
}

m2_ll <- function(f, u, arr_t_f, arr_t_u, h, w, sd_est){

    R_end <- as.data.frame(mapply(m2_sim, 
                                R=data$amt_init*fora_vol/1000, 
                                time=data$time/60/24, 
                                f=f/fora_vol,
                                u = 0, 
                                length=data$mm, 
                                gamma=gamma, 
                                Z=0,
                                arr_t_f=arr_t_f,
                                arr_t_u=arr_t_u,
                                ref_t=ref_t,
                                h=h*exp(w*data$temp),
                                temp=data$temp))
  colnames(R_end) <- "R_end"
  data$end <- R_end$R_end
  
  data %<>% mutate(resid = sqrt(end) - sqrt(amt_rem*fora_vol/1000))
  data %<>% mutate(sqrt_model_end = sqrt(end),
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
  
  exp_data <- dataset %>% filter(exposed==TRUE)
  
  I_end <- as.data.frame(mapply(m2_sim, 
                                R=exp_data$resource*life_vol/1000, 
                                time=exp_data$time, 
                                f=f/life_vol, 
                                u=u, 
                                length = exp_data$life_mm, 
                                gamma=gamma, 
                                Z=spore_conc*life_vol,
                                arr_t_f=arr_t_f,
                                arr_t_u=arr_t_u,
                                ref_t=ref_t,
                                h=h*exp(w*exp_data$temp),
                                temp=exp_data$temp))
  
  colnames(I_end) <- "I_end"
  
  exp_data %<>% cbind(I_end)
  
  exp_data %<>% 
    mutate(ll = mapply(inf_out, inf_status = inf_status, I_end=I_end))
  nll_sum <- nll_sum + -sum(exp_data$ll)
  return(nll_sum)
}

# model 2C - Resource-only ------------------------------------------------

m3_sim <- function(R, time, f, u, length, gamma, Z, ref_t, arr_t_f, h, temp){
  times <- seq(0, time, by=0.01)
  xstart <- c(R=R,
              S=1,
              I=0,
              Z=Z)
  params <- c(f=f,
              u=u,
              length = length,
              gamma = gamma,
              Z = Z,
              ref_t = ref_t,
              arr_t_f = arr_t_f,
              h = h,
              temp = temp)
  output <- as.data.frame(lsoda(y=xstart, times, m1_num_sol, params))
  if(Z == 0) {
    end_data <- slice_max(output, time)[,2]
  }
  else {
    end_data <- slice_max(output, time)[,4]  
  }
  
  return(end_data)
}

m3_ll <- function(f, u, arr_t_f, h, w, rho, sd_est){
  
  #exit conditions if returning nonsense values
  if(f<0 | h<0 | u<0 | arr_t_f<0 |
     (u/10000*exp(0.1/1000*rho)) <0 |
     (u/10000*exp(1.0/1000*rho)) <0) {
    
    return(NA)
    
  }
  
  else {
  
  R_end <- as.data.frame(mapply(m3_sim, 
                                R=data$amt_init*fora_vol/1000, 
                                time=data$time/60/24, 
                                f=f/fora_vol,
                                u = 0, 
                                length=data$mm, 
                                gamma=gamma, 
                                Z=0,
                                arr_t_f=arr_t_f,
                                ref_t=ref_t,
                                h=h*exp(w*data$temp),
                                temp=data$temp))
  colnames(R_end) <- "R_end"
  data$end <- R_end$R_end

  data %<>% mutate(resid = sqrt(end) - sqrt(amt_rem*fora_vol/1000))
  data %<>% mutate(sqrt_model_end = sqrt(end),
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
  
  
  exp_data <- dataset %>% filter(exposed==TRUE)
  
  
  
  I_end <- as.data.frame(mapply(m3_sim, 
                                R=exp_data$resource*life_vol/1000, 
                                time=exp_data$time, 
                                f=f/life_vol, 
                                u=u*exp(exp_data$resource*rho)/10000, 
                                length = exp_data$life_mm, 
                                gamma=gamma, 
                                Z=spore_conc*life_vol,
                                arr_t_f=arr_t_f,
                                ref_t=ref_t,
                                h=h*exp(w*exp_data$temp),
                                temp=exp_data$temp))
  
  colnames(I_end) <- "I_end"
  
  exp_data %<>% cbind(I_end)
  
  exp_data %<>% 
    mutate(ll = mapply(inf_out, inf_status = inf_status, I_end=I_end))
  nll_sum <- nll_sum + -sum(exp_data$ll)
  return(nll_sum)
  }
}

# model 2D - Additive -----------------------------------------------------

m4_num_sol <- function(t, x, params){
  R <- x[1]
  S <- x[2]
  I <- x[3]
  Z <- x[4]
  
  with(as.list(params),{
    dR <- -R*(S+I)*(f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma))/
      (1+f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma)*h*R*exp(w*temp))
    
    dS <- -Z*S*u*exp(arr_t_u*(1/ref_t - 1/temp))*exp(resource*rho)*
      (f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma))/
      (1+f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma)*h*R*exp(w*temp))
    
    dI <- Z*S*u*exp(arr_t_u*(1/ref_t - 1/temp))*exp(resource*rho)*
      (f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma))/
      (1+f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma)*h*R*exp(w*temp))
    
    dZ <- -Z*(S+I)*(f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma))/
      (1+f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma)*h*R*exp(w*temp))
    res <- c(dR, dS, dI, dZ)
    list(res)}
  )
}


m4_sim <- function(R, time, f, u, length, gamma, Z, ref_t, arr_t_f, arr_t_u, h, temp, resource, w, rho){
  times <- seq(0, time, by=0.01)
  xstart <- c(R=R,
              S=1,
              I=0,
              Z=Z)
  params <- c(f=f,
              u=u,
              length = length,
              gamma = gamma,
              Z = Z,
              ref_t = ref_t,
              arr_t_f = arr_t_f,
              arr_t_u = arr_t_u,
              h = h,
              temp = temp,
              resource = resource,
              w = w,
              rho = rho)
  output <- as.data.frame(lsoda(y=xstart, times, m4_num_sol, params))
  
  if(Z == 0) {
    end_data <- slice_max(output, time)[,2]
  }
  else {
    end_data <- slice_max(output, time)[,4]  
  }
  
  return(end_data)
}

m4_ll <- function(f, u, arr_t_f, arr_t_u, h, w, rho, sd_est){
  
  #exit conditions if returning nonsense values
  if(f<0 | h<0 | u<0 | arr_t_f<0 | arr_t_u<0 |
     (u/10000*exp(arr_t_u*((1/15)-(1/15)))*exp(0.1/1000*rho)) <0 |
     (u/10000*exp(arr_t_u*((1/15)-(1/25)))*exp(1.0/1000*rho)) <0) {
    
    return(NA)
    
  }
  
  else {
    
    
    R_end <- as.data.frame(mapply(m4_sim, 
                                  R=data$amt_init*fora_vol/1000, 
                                  time=data$time/60/24, 
                                  f=f/fora_vol,
                                  u = 0, 
                                  length=data$mm, 
                                  gamma=gamma, 
                                  Z=0,
                                  ref_t=ref_t,
                                  arr_t_f=arr_t_f,
                                  arr_t_u=arr_t_u,
                                  h=h,
                                  temp = data$temp,
                                  resource = data$resource,
                                  w = w,
                                  rho = rho))
    colnames(R_end) <- "R_end"
    data$end <- R_end$R_end
    
    data %<>% mutate(resid = sqrt(end) - sqrt(amt_rem*fora_vol/1000))
    data %<>% mutate(sqrt_model_end = sqrt(end),
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
    
    
    exp_data <- dataset %>% filter(exposed==TRUE)
    
    
    I_end <- as.data.frame(mapply(m4_sim, 
                                  R=exp_data$resource*life_vol/1000, 
                                  time=exp_data$time, 
                                  f=f/life_vol, 
                                  u=(u/10000), 
                                  length = exp_data$life_mm, 
                                  gamma=gamma, 
                                  Z=spore_conc*life_vol,
                                  ref_t=ref_t,
                                  arr_t_f=arr_t_f,
                                  arr_t_u=arr_t_u,
                                  h=h,
                                  temp=exp_data$temp,
                                  resource=exp_data$resource,
                                  w = w,
                                  rho = rho))
    
    colnames(I_end) <- "I_end"
    
    exp_data %<>% cbind(I_end)
    
    exp_data %<>% 
      mutate(ll = mapply(inf_out, inf_status = inf_status, I_end=I_end))
    nll_sum <- nll_sum + -sum(exp_data$ll)
    return(nll_sum)
    
  }
}

# model 2E - interactive --------------------------------------------------


m5_num_sol <- function(t, x, params){
  R <- x[1]
  S <- x[2]
  I <- x[3]
  Z <- x[4]
  
  with(as.list(params),{
    dR <- -R*(S+I)*(f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma))/
      (1+f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma)*h*R*exp(w*temp))
    
    dS <- -Z*S*u*exp(arr_t_u*(1/ref_t - 1/temp))*exp(resource*rho)*exp(resource*temp*phi)*
      (f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma))/
      (1+f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma)*h*R*exp(w*temp))
    
    dI <- Z*S*u*exp(arr_t_u*(1/ref_t - 1/temp))*exp(resource*rho)*exp(resource*temp*phi)*
      (f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma))/
      (1+f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma)*h*R*exp(w*temp))
    
    dZ <- -Z*(S+I)*(f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma))/
      (1+f*exp(arr_t_f*(1/ref_t - 1/temp))*(length^gamma)*h*R*exp(w*temp))
    res <- c(dR, dS, dI, dZ)
    list(res)}
  )
}


m5_sim <- function(R, time, f, u, length, gamma, Z, ref_t, arr_t_f, arr_t_u, h, temp, resource, w, rho, phi){
  times <- seq(0, time, by=0.01)
  xstart <- c(R=R,
              S=1,
              I=0,
              Z=Z)
  params <- c(f=f,
              u=u,
              length = length,
              gamma = gamma,
              Z = Z,
              ref_t = ref_t,
              arr_t_f = arr_t_f,
              arr_t_u = arr_t_u,
              h = h,
              temp = temp,
              resource = resource,
              w = w,
              rho = rho,
              phi = phi)
  output <- as.data.frame(lsoda(y=xstart, times, m5_num_sol, params))
  
  if(Z == 0) {
    end_data <- slice_max(output, time)[,2]
  }
  else {
    end_data <- slice_max(output, time)[,4]  
  }
  
  return(end_data)
}

m5_ll <- function(f, u, arr_t_f, arr_t_u, h, w, rho, phi, sd_est){
  
  if(f<0 | h<0 | u<0 | arr_t_f<0 | arr_t_u<0 |
     (u/100000*exp(arr_t_u*((1/15)-(1/15)))*exp(0.1/1000*rho)*exp(0.1/1000*15*phi)) <0 |
     (u/100000*exp(arr_t_u*((1/15)-(1/25)))*exp(1.0/1000*rho)*exp(1.0/1000*25*phi)) <0) {
    
    return(NA)
    
  } 
  
  else {  
    
    
    
    R_end <- as.data.frame(mapply(m5_sim, 
                                  R=data$amt_init*fora_vol/1000, 
                                  time=data$time/60/24, 
                                  f=f/fora_vol,
                                  u = 0, 
                                  length=data$mm, 
                                  gamma=gamma, 
                                  Z=0,
                                  ref_t=ref_t,
                                  arr_t_f=arr_t_f,
                                  arr_t_u=arr_t_u,
                                  h=h,
                                  temp = data$temp,
                                  resource = data$resource,
                                  w = w,
                                  rho = rho,
                                  phi = phi))
    colnames(R_end) <- "R_end"
    data$end <- R_end$R_end
    
    data %<>% mutate(resid = sqrt(end) - sqrt(amt_rem*fora_vol/1000))
    data %<>% mutate(sqrt_model_end = sqrt(end),
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
    
    
    exp_data <- dataset %>% filter(exposed==TRUE)
    
    
    I_end <- as.data.frame(mapply(m5_sim, 
                                  R=exp_data$resource*life_vol/1000, 
                                  time=exp_data$time, 
                                  f=f/life_vol, 
                                  u=u/100000, 
                                  length = exp_data$life_mm, 
                                  gamma=gamma, 
                                  Z=spore_conc*life_vol,
                                  ref_t=ref_t,
                                  arr_t_f=arr_t_f,
                                  arr_t_u=arr_t_u,
                                  h=h,
                                  temp=exp_data$temp,
                                  resource=exp_data$resource,
                                  w = w,
                                  rho = rho,
                                  phi = phi))
    
    colnames(I_end) <- "I_end"
    
    exp_data %<>% cbind(I_end)
    
    exp_data %<>% 
      mutate(ll = mapply(inf_out, inf_status = inf_status, I_end=I_end))
    nll_sum <- nll_sum + -sum(exp_data$ll)
    return(nll_sum)
    
  }
}


