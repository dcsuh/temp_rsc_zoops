## length correction

library(here)

source(here("base","src.R"))

#lengths were measured at different magnifications
#this script converts raw measurements at different magnifications into standard mm lengths

#conversions as follows
#at default magnification (5.6x), 1 unit is equal to 17.86 micron
#at 5x, 1 = 20
#at 4x, 1 = 25

mort <- read_csv(here("raw_data", "main_mort_data.csv"))
mort$length <- c(0)

for (i in 1:nrow(mort)){
  if (!is.na(mort$length_mag_correction[i])){
    if (mort$lm_corr_factor[i] == -1){
      mort$length[i] <- mort$length_RAW[i] * 20 / 1000
    }
    else if (mort$lm_corr_factor[i] == -2){
      mort$length[i] <- mort$length_RAW[i] * 25 / 1000
    }
  }
  else
    mort$length[i] <- mort$length_RAW[i] * 17.86 / 1000
  print(i)
}

saveRDS(mort, file = here("raw_data", "main_mort_edit.rds"))
