# plot precip data for NHC watershed

library(ggplot2)
library(tidyverse)
library(lubridate)

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
  pivot_longer(cols = precip_mm, names_to = "variable", values_to = "value")
  # group_by(year = as.numeric(substr(datetime, 1, 4))) %>%
  # filter(year < 2013) %>%
  # summarize(precip_mm = sum(precip_mm, na.rm = TRUE)) %>%
  # ungroup()
all <- bind_rows(nldas, precip) %>% 
  pivot_wider(names_from = "variable", values_from = "value")

plot(all$precip_mm, all$precipitation_amount)

p90 <- quantile(all$precipitation_amount, .9, na.rm = T)
pp <- all %>% 
  select(datetime, precip_mmd = precipitation_amount) %>%
  filter(!is.na(precip_mmd)) %>%
  group_by(year = year(datetime)) %>%
  summarize(zero_days = length(which(precip_mmd == 0)),
            cumulative_precip = sum(precip_mmd), 
            cum90 = sum(precip_mmd[which(precip_mmd >= p90)]),
            max_precip = max(precip_mmd, na.rm = T), 
            percent_extreme = cum90/cumulative_precip) %>%
  select(-cum90, -max_precip) %>%
  pivot_longer(cols= -year, names_to = "variable", values_to = "value")

ggplot(pp, aes(x = year, y = value)) +
  geom_point() +
  facet_wrap(.~variable, scales = "free_y", dir = "v", switch = "y") +
  geom_smooth(method = lm, lwd = 1, col = "black") 
  theme_minimal()


