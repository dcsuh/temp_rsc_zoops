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

saveRDS(boot_01, file = here("processed_data", "seq_data", "bootstrap_results.rds"))
