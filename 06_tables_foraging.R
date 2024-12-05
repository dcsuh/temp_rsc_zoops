#generate tables for foraging model results
#Daniel Suh

library(here)

source(here("base","src.R"))

#read data

fora_data <- readRDS(here("processed_data","foraging_raw.rds"))

m1a <- readRDS(here("processed_data", "mle", "m1_f_fit.rds"))
m1b <- readRDS(here("processed_data", "mle", "m2_f_fit.rds"))
m1c <- readRDS(here("processed_data", "mle", "m3_f_fit.rds"))
m1d <- readRDS(here("processed_data", "mle", "m4_f_fit.rds"))
m1e <- readRDS(here("processed_data", "mle", "m5_f_fit.rds"))


#AICc table for foraing data

fora_model_list <- 
  list(m1a, m1b, m1c, m1d, m1e)
mod_names <- c("(1A) Size-only", 
               "(1B) Temperature-only", 
               "(1C) Resource-only", 
               "(1D) Additive", 
               "(1E) Interactive")

aic_c_table_final <- AICctab(fora_model_list, 
                             logLik=T, 
                             base=T, 
                             mnames=mod_names, 
                             weights=T, 
                             nobs=nrow(fora_data))

aicc_flex <- 
  aic_c_table_final %>% 
  print() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "model") %>% 
  flextable()

if(dir.exists(here("figures")) == FALSE) {
  message("Welcome! Let's make some room for figures.")
  dir.create(here("figures")) 
} else {
  message("/figures exists! Proceeeding to save.")
}

save_as_docx("aicc" = aicc_flex, path = here("figures", "fora_table_results.docx"))