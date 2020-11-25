# plot precip data for NHC watershed

library(ggplot2)
library(tidyverse)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/code/")

nldas <- read_csv("data/nldas.csv") %>%
  dplyr::select(datetime = DateTime, value = '1', variable) %>% 
  filter(!variable %in% c("wind_speed",
                          "max_relative_humidity",
                          "min_relative_humidity",
                          "surface_downwelling_shortwave_flux_in_air"))

ggplot(nldas, aes(x=datetime, y = value)) +
  geom_line() +
  facet_wrap(.~variable, scales = "free_y" )

  
precip = read_csv('data/prism/prism_raw.csv') %>%
  as_tibble() %>%
  dplyr::select(datetime = DateTime, precip_mm = '1') %>%
  group_by(year = as.numeric(substr(datetime, 1, 4))) %>%
  filter(year < 2013) %>%
  summarize(precip_mm = sum(precip_mm, na.rm = TRUE)) %>%
  ungroup()

ggplot(precip, aes(year, precip_mm)) +
  geom_point() +
  geom_smooth()
