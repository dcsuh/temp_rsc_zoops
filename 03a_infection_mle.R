#estimate parameters for transmission via maximum likelihood
#Daniel Suh

library(here)


if(file.exists(here("processed_data", "mle", "m1_combined_fit.rds")) == TRUE |
   file.exists(here("processed_data", "mle", "m2_combined_fit.rds")) == TRUE |
   file.exists(here("processed_data", "mle", "m3_combined_fit.rds")) == TRUE |
   file.exists(here("processed_data", "mle", "m4_combined_fit.rds")) == TRUE |
   file.exists(here("processed_data", "mle", "m5_combined_fit.rds")) == TRUE) {
  
  message("Model parameter estimates already exist.")
  
} else {
  



source(here("base","src.R"))
source(here("03_infection.R"))



# maximum likelihood estimation -------------------------------------------

m1_fit <- 
  mle2(m1_ll, start=list(u = 0.0001, 
                         f=18, 
                         arr_t_f=30, 
                         h=17000, 
                         w=-0.2,
                         sd_est = 0.01),
       control=list(parscale = c(u = 0.0001, 
                                 f=1, 
                                 arr_t_f=10, 
                                 h=10000, 
                                 w=0.1,
                                 sd_est = 0.01)),
       skip.hessian=F, method="BFGS")


m2_fit <-
  mle2(m2_ll, start=list(u = 0.0001, 
                         f=18,
                         arr_t_f = 30,
                         arr_t_u = 65,
                         h=16000,
                         w=-0.1,
                         sd_est = 0.01), 
       control=list(parscale = c(u = 0.0001, 
                                 f=1,
                                 arr_t_f = 10,
                                 arr_t_u = 10,
                                 h=10000,
                                 w=0.1,
                                 sd_est = 0.01)),
       skip.hessian=F, method="BFGS")


m3_fit <- 
  mle2(m3_ll, start=list(u = 1, 
                         f=30, 
                         arr_t_f = 30,
                         h=16000,
                         w=0.1,
                         rho = 0.1,
                         sd_est = 0.01), 
       control=list(parscale = c(u = 1, 
                                 f=1, 
                                 arr_t_f = 1,
                                 h=10000,
                                 w=0.1,
                                 rho = 0.1,
                                 sd_est = 0.01), maxit = 10000),
       skip.hessian=F, method="BFGS")


m4_fit <-
  mle2(m4_ll, start=list(f = 10,
                         u=1, 
                         arr_t_f=30, 
                         arr_t_u=30,
                         h=10000,
                         w=0.1,
                         rho=0.1,
                         sd_est = 0.01),
       control=list(parscale = 
                      c(f = 10,
                        u=1, 
                        arr_t_f=30, 
                        arr_t_u=30,
                        h=10000,
                        w=0.1,
                        rho=0.1,
                        sd_est = 0.01), maxit = 10000),
       skip.hessian=F, method="BFGS")

m5_fit <-
  mle2(m5_ll, start=list(f = 10,
                         u=10, 
                         arr_t_f=30, 
                         arr_t_u=30,
                         h=100,
                         w=0.1,
                         rho=0.1,
                         phi=0.1,
                         sd_est = 0.01), 
       control=list(parscale = c(f = 10,
                                 u=10, 
                                 arr_t_f=10, 
                                 arr_t_u=10,
                                 h=100,
                                 w=0.1,
                                 rho=0.1,
                                 phi=0.1,
                                 sd_est = 0.01),
                    maxit=5000),
       skip.hessian=F, 
       method="BFGS")



# save fits ---------------------------------------------------------------

if(dir.exists(here("processed_data", "mle")) == FALSE) {
  message("Welcome! Let's make some room for the model estimates.")
  dir.create(here("processed_data", "mle")) 
} else {
  message("/processed_data/mle exists! Proceeeding to save.")
}

saveRDS(m1_fit, file = here("processed_data", "mle", "m1_combined_fit.rds"))
saveRDS(m2_fit, file = here("processed_data", "mle", "m2_combined_fit.rds"))
saveRDS(m3_fit, file = here("processed_data", "mle", "m3_combined_fit.rds"))
saveRDS(m4_fit, file = here("processed_data", "mle", "m4_combined_fit.rds"))
saveRDS(m5_fit, file = here("processed_data", "mle", "m5_combined_fit.rds"))

}
