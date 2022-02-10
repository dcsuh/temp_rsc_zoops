library(here)

source(here("base","src.R"))

treatments <- read_csv(here("raw_data", "treatments.csv"))
fitness <- read_csv(here("raw_data", "fitness.csv"))
mort <- read_csv(here("raw_data", "mortality.csv"))

