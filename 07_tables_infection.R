#generate tables for foraging model results
#Daniel Suh

library(here)

source(here("base","src.R"))
source(here("03_infection.R"))

#read data

m1 <- readRDS(here("processed_data", "mle", "m1_combined_fit.rds"))
m2 <- readRDS(here("processed_data", "mle", "m2_combined_fit.rds"))
m3 <- readRDS(here("processed_data", "mle", "m3_combined_fit.rds"))
m4 <- readRDS(here("processed_data", "mle", "m4_combined_fit.rds"))
m5 <- readRDS(here("processed_data", "mle", "m5_combined_fit.rds"))



# component likelihood functions ------------------------------------------

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
  
  f_nll <- nll_sum
  
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
  u_nll <- -sum(exp_data$ll)
  return(c(f_nll, u_nll))
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
  
  f_nll <- nll_sum
  
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
  u_nll <- -sum(exp_data$ll)
  return(c(f_nll, u_nll))
}


m3_ll <- function(f, u, arr_t_f, h, w, rho, sd_est){
  
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
  
  f_nll <- nll_sum
  
  
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
  u_nll <- -sum(exp_data$ll)
  return(c(f_nll, u_nll))
}



m4_ll <- function(f, u, arr_t_f, arr_t_u, h, w, rho, sd_est){
  
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
    
    f_nll <- nll_sum
    
    
    exp_data <- dataset %>% filter(exposed==TRUE)
    #  exp_data <- dataset_tmp %>% filter(exposed==TRUE)
    
    
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
    u_nll <- -sum(exp_data$ll)
    return(c(f_nll, u_nll))
  }
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
    
    f_nll <- nll_sum
    
    
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
    u_nll <- -sum(exp_data$ll)
    return(c(f_nll, u_nll))
  }
}


m1_output <- m1_ll(u = coef(m1)[1], 
                   f=coef(m1)[2], 
                   arr_t_f=coef(m1)[3], 
                   h=coef(m1)[4], 
                   w=coef(m1)[5],
                   sd_est=coef(m1)[6])

m2_output <- m2_ll(f = coef(m2)[1], 
                   u=coef(m2)[2], 
                   arr_t_f=coef(m2)[3], 
                   arr_t_u=coef(m2)[4], 
                   h=coef(m2)[5], 
                   w=coef(m2)[6],
                   sd_est=coef(m2)[7])

m3_output <- m3_ll(f = coef(m3)[1], 
                   u=coef(m3)[2], 
                   arr_t_f=coef(m3)[3], 
                   h=coef(m3)[4], 
                   w=coef(m3)[5],
                   rho=coef(m3)[6],
                   sd_est=coef(m3)[7])

m4_output <- m4_ll(f = coef(m4)[1], 
                   u=coef(m4)[2], 
                   arr_t_f=coef(m4)[3], 
                   arr_t_u=coef(m4)[4], 
                   h=coef(m4)[5], 
                   w=coef(m4)[6],
                   rho=coef(m4)[7],
                   sd_est=coef(m4)[8])

m5_output <- m5_ll(f = coef(m5)[1], 
                   u=coef(m5)[2], 
                   arr_t_f=coef(m5)[3], 
                   arr_t_u=coef(m5)[4], 
                   h=coef(m5)[5], 
                   w=coef(m5)[6],
                   rho=coef(m5)[7],
                   phi=coef(m5)[8],
                   sd_est=coef(m5)[9])



# corrected AIC -----------------------------------------------------------



fora_n <- nrow(data)
inf_n <- nrow(dataset %>% filter(exposed==T))


fora_params <- 5

get_combined_corrected_AIC <- function(model, output){
  (2*fora_params + 2*output[1]) + (2*fora_params^2 + 2*fora_params)/(fora_n - fora_params - 1) + 
    
    (2*(length(coef(model))-fora_params) + 2*output[2]) + 
    (2*(length(coef(m1))-fora_params)^2 + 2*(length(coef(m1))-fora_params))/
    (inf_n - (length(coef(m1))-fora_params) - 1)
}


m1_AICc <- get_combined_corrected_AIC(m1, m1_output)

m2_AICc <- get_combined_corrected_AIC(m2, m2_output)

m3_AICc <- get_combined_corrected_AIC(m3, m3_output)

m4_AICc <- get_combined_corrected_AIC(m4, m4_output)

m5_AICc <- get_combined_corrected_AIC(m5, m5_output)

m1_AICc <- as.numeric(round(m1_AICc, digits = 1))
m2_AICc <- as.numeric(round(m2_AICc, digits = 1))
m3_AICc <- as.numeric(round(m3_AICc, digits = 1))
m4_AICc <- as.numeric(round(m4_AICc, digits = 1))
m5_AICc <- as.numeric(round(m5_AICc, digits = 1))

m1_df <- length(coef(m1)) 
m2_df <- length(coef(m2)) 
m3_df <- length(coef(m3)) 
m4_df <- length(coef(m4)) 
m5_df <- length(coef(m5)) 

m1_ll <- round(-sum(m1_output), digits = 1)
m2_ll <- round(-sum(m2_output), digits = 1)
m3_ll <- round(-sum(m3_output), digits = 1)
m4_ll <- round(-sum(m4_output), digits = 1)
m5_ll <- round(-sum(m5_output), digits = 1)

m1_dll <- m1_ll - m1_ll
m2_dll <- m2_ll - m1_ll
m3_dll <- m3_ll - m1_ll
m4_dll <- m4_ll - m1_ll
m5_dll <- m5_ll - m1_ll

m5_dAICc <- m5_AICc-m5_AICc
m4_dAICc <- m4_AICc-m5_AICc
m3_dAICc <- m3_AICc-m5_AICc
m2_dAICc <- m2_AICc-m5_AICc
m1_dAICc <- m1_AICc-m5_AICc

sum_deltaAIC <- 
  exp(-0.5*(m1_dAICc)) + 
  exp(-0.5*(m2_dAICc)) + 
  exp(-0.5*(m3_dAICc)) + 
  exp(-0.5*(m4_dAICc)) + 
  exp(-0.5*(m5_dAICc))

m1_weight <- exp(-0.5*(m1_AICc-m5_AICc))/sum_deltaAIC
m2_weight <- exp(-0.5*(m2_AICc-m5_AICc))/sum_deltaAIC
m3_weight <- exp(-0.5*(m3_AICc-m5_AICc))/sum_deltaAIC
m4_weight <- exp(-0.5*(m4_AICc-m5_AICc))/sum_deltaAIC
m5_weight <- exp(-0.5*(m5_AICc-m5_AICc))/sum_deltaAIC

mod_names <- c("(2A) Independent", 
               "(2B) Temperature-only", 
               "(2C) Resource-only", 
               "(2D) Additive", 
               "(2E) Interactive")

mod_results_combined_df <- tibble(model = mod_names,
                                  logLik = c(m1_ll, m2_ll, m3_ll, m4_ll, m5_ll),
                                  AICc = c(m1_AICc, m2_AICc, m3_AICc, m4_AICc, m5_AICc),
                                  dLogLik = c(m1_dll, m2_dll, m3_dll, m4_dll, m5_dll),
                                  dAICc = c(m1_dAICc, m2_dAICc, m3_dAICc, m4_dAICc, m5_dAICc),
                                  df = c(m1_df, m2_df, m3_df, m4_df, m5_df),
                                  weight = round(c(m1_weight, m2_weight, m3_weight, m4_weight, m5_weight), digits = 4))

mod_results_combined_flex <- 
  flextable(mod_results_combined_df %>% 
              arrange(., desc(weight))) %>% 
  colformat_double(j = "AICc", big.mark = "")

if(dir.exists(here("figures")) == FALSE) {
  message("Welcome! Let's make some room for figures.")
  dir.create(here("figures")) 
} else {
  message("/figures exists! Proceeeding to save.")
}

save_as_docx("aicc_combined" = mod_results_combined_flex, path = here("figures", "combined_table_results.docx"))

