## Base script for installing/loading packages and setting plot theme

# this function loads/installs packages
install_if_necessary <- function(x) {
  if(!require(x, character.only = TRUE)){
    install.packages(x)
    library(x, character.only = TRUE)
  } else {
    library(x, character.only = TRUE)  
  }
}


# all the necessary packages go here
pckg_names <- c("magrittr", 
                "here", 
                "lubridate", 
                "epitools",
                "ggpubr",
                "rstatix",
                "lhs",
                "deSolve",
                "bbmle",
                "ggnewscale",
                "flextable",
                "AICcmodavg",
                "egg",
                "patchwork",
                "tidyverse")

# load/install packages
lapply(pckg_names, 
       install_if_necessary)



# setting the plot theme for the project
proj_theme <- theme_set(theme_bw(base_size = 20))


## Intended target dims
golden <- 1.618
outwidth <- c(8.7, 11.4, 17.8)
unit <- "in"

#Some useful vectors
const_temp <- c("15", "20", "25")
