#local functions
#Daniel Suh

#get AIC from mle outputs
getAIC <- function(model){
  return(2*length(coef(model)) + (summary(model)@m2logL))
}


#get infection prevalence
inf_out <- function(inf_status, I_end){
  dbinom(x = inf_status, size = 1, prob=I_end, log=T)
}