#local functions
#Daniel Suh

#get AIC from mle outputs
getAIC <- function(model){
  return(2*length(coef(model))-(2*summary(model)@m2logL/-2))
}
