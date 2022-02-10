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
pckg_names <- c("tidyverse", "magrittr", "here", "lubridate")

# load/install packages
lapply(pckg_names, 
       install_if_necessary)



# setting the plot theme for the project
theme_set(theme_minimal(base_size = 13))


## Intended target dims
outwidth <- c('8.7cm', '11.4cm', '17.8cm')


# load packages
