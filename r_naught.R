
library(here)

source(here("base","src.R"))

lt_full <- readRDS(here("processed_data", "lt_full_summary.rds"))


lt.summary_factors %>% filter(species=="D", temp_id %in% const_temp) %>% 
  ggplot(.,aes(x=trt, y=S.R_naught)) + 
  geom_point() + 
  geom_linerange(aes(ymin=S.R_naught.025, ymax=S.R_naught.975)) + 
  theme(axis.text.x = element_text(angle=90)) + 
  geom_hline(yintercept=1)

