#estimate parameters for foraging rate via maximum likelihood
#Daniel Suh

library(here)




if(file.exists(here("processed_data", "mle", "m1_f_fit.rds")) == TRUE |
   file.exists(here("processed_data", "mle", "m2_f_fit.rds")) == TRUE |
   file.exists(here("processed_data", "mle", "m3_f_fit.rds")) == TRUE |
   file.exists(here("processed_data", "mle", "m4_f_fit.rds")) == TRUE |
   file.exists(here("processed_data", "mle", "m5_f_fit.rds")) == TRUE) {
  
  message("Model parameter estimates already exist.")
  
  } else {
    


source(here("base","src.R"))
source(here("01_foraging.R"))

# maximum likelihood estimation -------------------------------------------


m1_f_fit <- mle2(minuslogl=m1_ll, 
                 skip.hessian=T, 
                 start=list(f=10,
                            sd_est = 0.01),
                 method = "L-BFGS-B",
                 lower = c(0, 0),
                 upper = c(Inf, Inf),
                 control=list(parscale=c(f=10, sd_est = 0.001), maxit=10000))

m2_f_fit <- mle2(m2_ll,
                 skip.hessian = T,
                 start=list(f=5,
                            arr_t=50,
                            sd_est=0.001),
                 method = "L-BFGS-B",
                 lower = c(0, 0, 0),
                 upper = c(Inf, Inf, Inf),
                 control=list(parscale=c(f=5, arr_t=50, sd_est=0.001), maxit=10000))

m3_f_fit <- mle2(m3_ll,
                 skip.hessian = T,
                 start=list(f=40,
                            h=100,
                            sd_est=0.001),
                 method = "L-BFGS-B",
                 lower = c(0, 0, 0),
                 upper = c(Inf, Inf, Inf),
                 control=list(parscale=c(f=40, h=100, sd_est=0.001), maxit=10000))

m4_f_fit <- mle2(m4_ll,
                 skip.hessian = T,
                 start=list(f=5, arr_t=50, h=100, sd_est=0.001),
                 method = "L-BFGS-B",
                 lower = c(0, 0, 0, 0),
                 upper = c(Inf, Inf, Inf, Inf),
                 control=list(parscale=c(f=5, arr_t=50, h=100, sd_est=0.001), maxit=10000))

m5_f_fit <- mle2(m5_ll,
                 skip.hessian = T,
                 start=list(f=10, arr_t=50, h=10, w=0.1, sd_est = 0.001),
                 method = "L-BFGS-B",
                 lower = c(f=0, arr_t=0, h=0, w=-Inf, sd_est = 0),
                 upper = c(Inf, Inf, Inf, Inf, Inf),
                 control=list(parscale=c(f=10, arr_t=50, h=10, w=0.1, sd_est = 0.001), maxit=10000))

# save outputs ------------------------------------------------------------

if(dir.exists(here("processed_data", "mle")) == FALSE) {
  message("Welcome! Let's make some room for the model estimates.")
  dir.create(here("processed_data", "mle")) 
} else {
  message("/processed_data/mle exists! Proceeeding to save.")
}

saveRDS(m1_f_fit, file = here("processed_data", "mle", "m1_f_fit.rds"))
saveRDS(m2_f_fit, file = here("processed_data", "mle", "m2_f_fit.rds"))
saveRDS(m3_f_fit, file = here("processed_data", "mle", "m3_f_fit.rds"))
saveRDS(m4_f_fit, file = here("processed_data", "mle", "m4_f_fit.rds"))
saveRDS(m5_f_fit, file = here("processed_data", "mle", "m5_f_fit.rds"))

}
