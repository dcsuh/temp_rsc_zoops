## script for defining functions used elsewhere in project

#calculate prevalence using threshold spore yield

get_threshold <- function(threshold, data = mort, factors = treatment_factors){
  data %<>% mutate(inf_adj = ifelse(spore_yield>threshold, 1, 0))
  prev_adj <- data %>% drop_na(inf_adj) #remove observations with NA for inf_adj
  prev_adj %<>% group_by(ID) %>% summarize(n = n(), prev = sum(inf_adj)/n()) 
  prev_adj %<>% left_join(.,treatment_factors)
  return(prev_adj)
}

