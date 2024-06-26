#data prep for final fits


# packages ----------------------------------------------------------------


library(here)
source(here("base","src.R"))


# read data ---------------------------------------------------------------


fitness <- read.csv(here("raw_data/main_fitness_edit.csv")) #must be read.csv
fitness %<>% dplyr::select(tube, temp_id, REMOVED, KBP, final_date) #remove unnecessary variables
mort <-  read_csv(here("raw_data/main_mort_edit.csv")) #mortality data
lengths <- read_csv(here("raw_data", "day5_length.csv")) #get average lengths
lengths %<>% filter(temp_id %in% const_temp)


dataset <- left_join(mort,fitness) #life table data



data <- readRDS(here("processed_data","foraging_raw.rds")) #foraging rate data
data_summ <- readRDS(here("processed_data", "foraging.rds")) #foraging rate data summary

data %<>% mutate(treatment_ID = paste(temp, resource, sep = "_"))

data %<>% mutate(resource_tot = resource*vol) #resource_tot is the total amount of resources in our tube. resource concentration * total volume of tube



# global variables --------------------------------------------------------



treatment_IDs <- unique(data$treatment_ID)
resource_IDs <- unique(data$resource)
temp_IDs <- unique(data$temp)


fora_vol <- 15 #mL
life_vol <- 50 #mL
spore_conc <- 200 #spores/mL
gamma <- 2 
ref_t <- 15 #celsius



# data prep ---------------------------------------------------------------




dataset %<>% filter(species == "D") %>% #only daphnia for this analysis
  filter(temp %in% const_temp) #only constant temp for this analysis

dataset %<>% mutate(birthdate = ifelse(species == "daphnia", "4/5/22", "4/6/22"), 
                    lifespan = as.numeric(mdy(final_date) - mdy(birthdate)),
                    temp = as.numeric(temp),
                    treatment = paste(temp, resource, sep = "_")) %>% 
  filter(is.na(male)) %>% #remove males
  filter(is.na(missing)) %>% #remove missing
  filter(!is.na(inf)) #remove NAs for inf for estimating beta. sometimes they died too young to tell

dataset %<>% mutate(inf_status = inf, dead = ifelse(is.na(REMOVED) & is.na(KBP), 1, 0))


dataset %<>% mutate(spore_exposure = 200, #spores/mL
                    uninf = 1-inf_status,
                    time = 1, #duration of exposure in days
                    trt = paste(temp, resource, species, sep = "_"))




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

length_summ %<>% add_row(temp_id=as.character(25), resource=as.numeric(0.5), life_mm=1.16) #midpoint between 1.0 and 0.1

dataset %<>% left_join(length_summ)

length_summ %<>% mutate(ID = paste(temp_id, resource, sep="_"))




mean_length_summ <- data_summ %>% dplyr::select(temp, resource, mm_mean)
data %<>% left_join(., mean_length_summ)
data %<>% mutate(mm = if_else(!is.na(mm), mm, mm_mean))

mean_length_summ %<>% rename(fora_mm = "mm_mean")
dataset %<>% left_join(., mean_length_summ)




# global functions --------------------------------------------------------


inf_out <- function(inf_status, I_end){
  dbinom(x = inf_status, size = 1, prob=I_end, log=T)
}



getAIC <- function(model){
  return(2*length(coef(model))-(2*summary(model)@m2logL/-2))
}



# solver functions --------------------------------------------------------

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


m4_sim <- function(R, time, f, u, length, gamma, Z, ref_t, arr_t_f, arr_t_u, h, temp){
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


m5_ll <- function(f, u_15, u_20, u_25, arr_t_f, h, w, rho_15, rho_20, rho_25){
  data <- tmp
  R_end <- as.data.frame(mapply(m3_sim, 
                                R=data$amt_init*fora_vol/1000, 
                                time=data$time/60/24, 
                                f=f/fora_vol,
                                u = 0, 
                                length=data$mm, gamma=gamma, Z=0,
                                ref_t=ref_t,
                                arr_t_f=arr_t_f,
                                h=h*exp(w*data$temp),
                                temp = data$temp))
  colnames(R_end) <- "R_end"
  data$end <- R_end$R_end
  data %<>% mutate(resid = log(end) - log(amt_rem*fora_vol/1000))
  
  nll <- dnorm(data$resid, 
               mean = 0, 
               sd = sd(data$resid, na.rm = T), 
               log = T)
  nll_sum <- -sum(nll, na.rm = T)
  
  tmp_dataset %<>% filter(exposed==TRUE)
  
  u_vec <- c(u_15, u_20, u_25)
  rho_vec <- c(rho_15, rho_20, rho_25)
  
  for(i in 1:3){
    exp_data <- tmp_dataset %>% filter(temp == temp_IDs[i])
    
    I_end <- as.data.frame(mapply(m3_sim, 
                                  R=exp_data$resource*life_vol/1000, 
                                  time=exp_data$time, 
                                  f=f/life_vol, 
                                  u=u_vec[i]*exp(exp_data$resource*rho_vec[i]), 
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
  }
  return(nll_sum)
}
