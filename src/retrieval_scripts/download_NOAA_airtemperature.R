# Download NOAA data
# Plot trend over time, compare to NHC water temperatures
# Calculate degree days (above 90 F, 32 C), mean annual temp, peak temperature

# Created 10.05.2020
# AMC

library(tidyverse)
library(streamMetabolizer)
library(lubridate)
setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl")

# Get air temperature data from the NHC watershed from 1968-2019:
# This sp function uses the closest NOAA data to my site
### Mike has a function that will interpolate all of the NOAA data to give me
#     a more exact estimate at a point, might be worth doing that at some point
### It would also be useful to make a custom version of this funciton to get
#     other variables at custom intervals and also to return the distance 
#     to the NOAA station selected

### change this in the future so that you download the rest of the noaa data 

# site_dat <- read_csv("data/NHCsite_metadata.csv")
# long <- site_dat$longitude[site_dat$sitecode=="NHC"]
# lat <- site_dat$latitude[site_dat$sitecode=="NHC"]
# enddate <- site_dat$enddate.UTC[site_dat$sitecode=="NHC"] %>%
#   as.POSIXct(format = "%m/%d/%Y %H:%M", tz = "UTC")
# startdate <- as.POSIXct("1968-01-01 00:00:00", tz="UTC")
# noaa_dat <- StreamPULSE:::FindandCollect_airpres(lat, long,
#                                                  startdate,
#                                                  enddate)
# write_rds(noaa_dat, "data/noaa_air_temp_raw.rds")
noaa_dat <- read_rds("data/noaa_air_temp_raw.rds")

noaa_dat <- noaa_dat %>%
  mutate(datetime = with_tz(DateTime_UTC, tzone="EST"),
         date = as.Date(datetime, tz = "EST")) %>% 
  select(-DateTime_UTC)
# summarize this into daily values
daily <- noaa_dat %>%
  select(-datetime) %>%
  group_by(date) %>%
  summarize(temp_mean = mean(air_temp, na.rm = T),
            temp_min = min(air_temp, na.rm = T))
  
# time series decompositon that plots the seasonal component separated from 
#   the trend over time.

plot(daily$date, daily$temp_mean)
dd <- decompose(ts(daily$temp_mean, deltat = 1/365))
plot(dd)

decomp_noaa <- data.frame(air_trend = dd$trend,
                          date = daily$date)

# Import data from NHC, compare the water temperature to the air temperatures
### make sure to redownload when Steve has uploaded more NHC data. 
#StreamPULSE::query_available_data("NC","NHC")

# nhc_mega <- StreamPULSE::request_data("NC_NHC", variables=c("WaterTemp_C", "DO_mgL"))
# write_rds(nhc_mega, "data/NHC_watertemp.rds")
nhc_mega <- read_rds("data/NHC_watertemp.rds")
nhc <- nhc_mega$data %>%
  mutate(datetime = with_tz(DateTime_UTC, tzone="EST")) %>% 
  filter(! flagtype %in% c("Bad Data", "Questionable")) %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  select(datetime, water_temp = WaterTemp_C, DO_mgL, flagtype) %>%
  as_tibble()

#plot(nhc$datetime, nhc$water_temp, type="l", ylim = c(0,30))

# pair NHC and NOAA 15 minute dataframes

alldat <- full_join(nhc, noaa_dat, by="datetime")
plot(alldat$air_temp, alldat$water_temp, pch=20)
points(alldat$air_temp[alldat$flagtype == "Interesting"],
       alldat$water_temp[alldat$flagtype == "Interesting"], 
       pch = 20, col = "red")


# compare daily temperature means
all_daily <- alldat %>%
  select(-datetime, -flagtype, -air_kPa) %>%
  group_by(date) %>%
  # summarize_all(list(mean = ~mean(., na.rm = TRUE),
  #                    max = ~max(., na.rm = TRUE),
  #                    min = ~min(., na.rm = TRUE))) %>%
  summarize_all(list(mean = mean, max = max, min = min)) %>%
  ungroup() %>%
  mutate(year = year(date)) %>%
  filter(! year %in% c(1968, 2020))

# find days that have hypoxia
temp_thresh <- 30

hypox_days <- all_daily %>% 
  #filter(DO_mgL_min <= 2) %>%
  filter(air_temp_max >= temp_thresh) %>%
  select(DO_mgL_min, air_temp_max, water_temp_max)

pr_hpx <- length(which(all_daily$DO_mgL_min <= 2))/length(which(!is.na(all_daily$DO_mgL_min)))
pr_hpx_30deg <- length(which(hypox_days$DO_mgL_min <= 2))/length(which(!is.na(hypox_days$DO_mgL_min)))
pr_hpx_32deg <- length(which(hypox_days$DO_mgL_min <= 2))/length(which(!is.na(hypox_days$DO_mgL_min)))

hist(hypox_days$air_temp_max)


peak_temp_days <- all_daily %>%
  filter(air_temp_max > temp_thresh) %>%
  group_by(year) %>%
  summarize(days_above_temp = length(air_temp_max)) %>%
  ungroup()

annual <- all_daily %>%
  group_by(year) %>%
  summarize(mean_water_temp = mean(water_temp_mean, na.rm = T),
            peak_water_temp = max(water_temp_mean, na.rm=T),
            mean_air_temp = mean(air_temp_mean, na.rm = T),
            peak_air_temp = max(air_temp_mean, na.rm=T)) %>%
  ungroup()
annual <- left_join(annual, peak_temp_days)

long_annual <- annual %>%
  select(-mean_water_temp, -peak_water_temp) %>%
  rename(days_above_30C = days_above_temp) %>%
  pivot_longer(cols = -year, names_to = "trend", values_to = "temperature")


png(width=7, height=6, units='in', type='cairo', res=300,
    filename='figures/air_temp_trends.png')
  
  ggplot(long_annual, mapping = aes(year, temperature)) +
    geom_line() +
    facet_grid(long_annual$trend~., scales = "free")

dev.off()

plot(all_daily$air_temp, all_daily$water_temp, pch=20)
ll <- lm(all_daily$water_temp ~ all_daily$air_temp)
abline(a = ll$coefficients[1], b = ll$coefficients[2])

plot(annual$year, annual$mean_air_temp, type = "l", ylim = c(13,33))
lines(annual$year, annual$max_air_temp, lty = 2)
# reconstruct relationship removing days where airtemp <=-5 degrees, because 
# this is where the watertemp hits 0 according to above relationship

# all_daily <- all_daily[all_daily$air_temp > -5,]
# plot(all_daily$air_temp, all_daily$water_temp, pch=20)
# ll <- lm(all_daily$water_temp ~ all_daily$air_temp)
# abline(a = ll$coefficients[1], b = ll$coefficients[2])

# remove seasonal component from the time series, are the residuals correlated?

dd_wat <- decompose(ts(all_daily$water_temp, deltat = 1/365))
dd_air <- decompose(ts(all_daily$air_temp, deltat = 1/365))
plot(dd_wat)

detrended_dat <- data.frame(date = all_daily$date,
                            wat_trend = dd_wat$trend,
                            air_trend = dd_air$trend)
detrended_dat <- detrended_dat[which(!is.na(detrended_dat$wat_trend)),]

plot(detrended_dat$air_trend, detrended_dat$wat_trend)
trendmod <- lm(detrended_dat$wat_trend ~ detrended_dat$air_trend)
abline(a = trendmod$coefficients[1], b = trendmod$coefficients[2])

# png(width=7, height=6, units='in', type='cairo', res=300,
#     filename='figures/air_water_temp_trends.png')
# 
#   plot(detrended_dat$date, detrended_dat$air_trend, 
#        type = "l", ylab = "temp C", xlab = "date")
#   lines(detrended_dat$date, detrended_dat$wat_trend, col = "steelblue")
#   lines(decomp_noaa$date, decomp_noaa$air_trend, lty = 2)
#   legend("bottomright", 
#          c("50 year air trend", "4 year air trend", "4 year water trend"), 
#          lty = c(2, 1, 1), col = c("black","black","steelblue"), 
#          lwd = 1.2, bty = "n")
# 
# dev.off()

air_temp <- left_join(daily, decomp_noaa)
write_csv(air_temp, "data/noaa_air_temp.csv")

decomp_noaa$wat_trend <- decomp_noaa$air_trend * trendmod$coefficients[2] + 
  trendmod$coefficients[1]

