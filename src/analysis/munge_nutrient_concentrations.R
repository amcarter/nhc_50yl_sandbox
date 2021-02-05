# Compare nutrients measured in Hall to nutrients in NHC today

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/code/")
library(tidyverse)
library(lubridate)

# 1. Load Hall Data ####
# table 26 has TP data

t26 <- read_csv("../data/hall/hall_table_26.csv") %>%
  dplyr::rename(site = location)
t13 <- read_csv("../data/hall/hall_table_13_p.csv")
t14 <- read_csv("../data/hall/hall_table_14_nitrogen.csv") %>%
  mutate(no2_mgl = ifelse(no2_mgl == "nd", 0, as.numeric(no2_mgl)), 
         no3_mgl = ifelse(no3_mgl =="nd", 0, as.numeric(no3_mgl)))
hall <- full_join(t26, t13, by = c("date", "site", "TP_mgl")) %>%
  full_join(t14, by = c("date", "site"))

dat_mar <- read_csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/longitudinal_sampling/longitudinal_sensor_data/NHCLongitudinalDO_20190308.csv") %>%
  mutate(depth_cm = as.numeric(depth_cm))
dat_oct <- read_csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/longitudinal_sampling/longitudinal_sensor_data/NHCLongitudinalDO_20191009.csv")

now <- bind_rows(dat_mar, dat_oct) %>%
  select(date, distance_m, depth_cm, temp_C, DO_mgL, SpC_uscm, 
         NO3.N_mgl, NPOC_mgl, TN_mgl, NH4.N_mgl, PO4.P_mgl)

summary(hall)
summary(now)
