Code and Data for "Temperature and resources interact to affect transmission via host foraging rate and susceptibility" by Daniel C. Suh, Katie Schroeder, and Alexander T. Strauss

Scripts 00 through 05 perform all data cleaning and analysis.
Scripts 06 through 09 generate figures and statistical results.

Necessary packages and some functions are specified in scripts in the "base" folder
Figure components and tables are in the "figures" folder. Manuscript versions of figures with powerpoint file used for formatting are in the subfolder "manuscript_versions".


See below for a brief overview of each script:

00_data_processing.R
Converts raw data into usable format for future analyses
Outputted data are saved in the "processed_data" folder

01_foraging.R
Defines model functions for foraging rate analysis

01a_foraging_mle.R
Estimates model parameters via maximum likelihood
Point estimates are saved in the "mle" folder under "processed_data"
This scripts sources 01_foraging.R

02_foraging_data.R
Simulates data using ML estimates 
Data ares saved in the "seq_data" folder under "processed_data"
This scripts sources 01_foraging.R

03_infection.R
Defines model functions for foraging and transmission

03a_foraging.R
Estimates model parameters via maximum likelihood
Point estimates are saved in the "mle" folder under "processed_data"
This script sources 03_infection.R

04_infection_data.R
Simulates data using ML estimates 
Data ares saved in the "seq_data" folder under "processed_data"
This scripts sources 03_foraging.R

05_bootstrap.R
Defines function and generates bootstrapped confidence interval estimates
Simulates confidence intervals for handling time, foraging rate, and per parasite susceptibility
Generates Table S2 in supplementary materials
This scripts sources 03_foraging.R

06_tables_foraging.R
Generates AIC table for foraging models
Table 1A in manuscript

07_tables_infection.R
Generates AIC table for transmission models
Table 1B in manuscript
This scripts sources 03_foraging.R

08_figures.R
Generates figure components for constructing manuscript figures
Final manuscript versions of figures were formatted using Microsoft Powerpoint
Figures 1-5 in manuscript

09_stats.R
Performs traditional statistics
Reported in text and in figure 1 in manuscript
