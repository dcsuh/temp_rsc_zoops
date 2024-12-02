#generate bootstrapped confidence intervals
#Daniel Suh

library(here)

source(here("base","src.R"))
source(here("03_infection.R"))

#read data
m5 <- readRDS(here("processed_data", "mle", "m5_combined_fit.rds"))


dataset %<>% filter(exposed==T)
life_data <- dataset
fora_data <- data

m5_coef <- coef(m5)
m5_f <- as.numeric(m5_coef[1])
m5_u <- as.numeric(m5_coef[2])/100000
m5_arr <- as.numeric(m5_coef[3])
m5_arr_u <- as.numeric(m5_coef[4])
m5_h <- as.numeric(m5_coef[5])
m5_w <- as.numeric(m5_coef[6])
m5_rho <- as.numeric(m5_coef[7])
m5_phi <- as.numeric(m5_coef[8])
m5_sd_est <- as.numeric(m5_coef[9])

#1005 iterations were completed
iterations <- 1

boot_01 <- list()

for(i in 1:iterations){
  boot_data <- fora_data %>% 
    group_by(treatment_ID) %>%
    slice_sample(., prop=1, replace=T) %>%
    ungroup()
  
  boot_dataset <- life_data %>% 
    group_by(treatment) %>%
    slice_sample(., prop=1, replace=T) %>%
    ungroup()
  
  boot_01[[i]] <-mle2(m5_ll, start=list(u = m5_u*100000,
                                        f=m5_f, 
                                        arr_t_f=m5_arr, 
                                        arr_t_u=m5_arr_u,
                                        h=m5_h,
                                        w=m5_w,
                                        rho=m5_rho,
                                        phi=m5_phi,
                                        sd_est=m5_sd_est), 
                      control=list(parscale = c(u = m5_u*100000,
                                                f=m5_f, 
                                                arr_t_f=m5_arr, 
                                                arr_t_u=m5_arr_u,
                                                h=m5_h,
                                                w=m5_w,
                                                rho=m5_rho,
                                                phi=m5_phi,
                                                sd_est=m5_sd_est),
                                   maxit=5000,
                                   reltol=0.00001),
                      skip.hessian=F, 
                      method="Nelder-Mead")
  print(i)
  end_time <- Sys.time()
  print(end_time-start_time)
}


#format bootstrap data

#mod_boot <- boot_01
mod_boot <- readRDS(here("processed_data", "m2E_bootstraps.rds"))

m5_fit <- readRDS(here("processed_data", "mle", "m5_combined_fit.rds"))

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

mod_quantiles %<>% pivot_longer(cols = f:sd_est) %>% pivot_wider(., names_from="id")





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


f_seq <- seq_data %>%
  mutate(f = NA,
         f_lower = NA,
         f_upper = NA)


f_est <- mod_quantiles$est[1]
arr_est <- mod_quantiles$est[3]
ref_t <- 15


get_fora <- function(temp, resource, bound, mm){
  tmp <- mod_boot_coefs %>% 
    mutate(fora = (f*exp(arr*(1/ref_t - 1/temp))*(mm^2))/
             (1+f*exp(arr*(1/ref_t - 1/temp))*
                (mm^2)*
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


f_seq %<>% mutate(f = (f_est*exp(arr_est*(1/ref_t - 1/temp))*(length^2))/
                    (1+f_est*exp(arr_est*(1/ref_t - 1/temp))*
                       (length^2)*
                       h_est*exp(w_est*temp)*
                       resource/1000),
                  f_lower = mapply(get_fora, temp = temp, resource = resource, mm = length, bound = "lower"),
                  f_upper = mapply(get_fora, temp = temp, resource = resource, mm = length, bound = "upper"))



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


saveRDS(boot_01, file = here("processed_data", "m2E_bootstraps.rds"))
saveRDS(mod_quantiles, file = here("processed_data", "m2E_bootstrap_quantiles.rds"))
saveRDS(h_seq, here("processed_data", "seq_data", "h_ci_bootstrap.rds"))
saveRDS(f_seq, here("processed_data", "seq_data", "f_ci_bootstrap.rds"))
saveRDS(u_seq, here("processed_data", "seq_data", "u_ci_bootstrap.rds"))




