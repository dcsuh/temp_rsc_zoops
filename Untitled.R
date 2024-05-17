
f <- coef(m2_f_fit)[1]
arr_t <- coef(m2_f_fit)[2]

data <- tmp
m2_end <- as.data.frame(mapply(m2_sim, R=data$amt_init*data$vol/1000, time=data$time/60/24, f=f/data$vol,
                               length=data$mm, gamma=2, arr_t=arr_t, ref_t=15, temp=data$temp))
colnames(m2_end) <- "m2_end"
data$end <- m2_end$m2_end
data %<>% mutate(resid = log(end) - log(amt_rem*vol/1000))

nll <- dnorm(data$resid, 
             mean = 0, 
             sd = sd(data$resid, na.rm = T), 
             log = T)
-sum(nll, na.rm = T)


data %<>% mutate(actual_end = amt_rem*vol/1000)

#data %>% ggplot(., aes(x=end, y=actual_end)) + geom_point()


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


df <- data_summ %>% dplyr::select(temp, resource, mm_mean, time_mean, conc_mean)


df %<>% mutate(r_end = mapply(m2_sim,
                               R = resource*vol/1000,
                               time = time_mean,
                               f=f/vol,
                               length = mm_mean,
                               gamma = 2,
                               ref_t = 15,
                               arr_t = arr_t,
                               temp = temp))

df %<>% mutate(r_consumed = (resource*vol/1000)-r_end)

df %>% ggplot(., aes(x=resource, y=r_consumed, color = as.factor(temp))) + 
  geom_point(shape=1) + 
  geom_point(aes(x = resource, y = (resource-conc_mean)*vol/1000, color = as.factor(temp)))
