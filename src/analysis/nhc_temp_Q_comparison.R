#################################
# Temperature and discharge comparison 
# between Hall 1972 data and modern NHC data

# a carter
# 2020-10-14


library(tidyverse)
library(lubridate)
library(zoo)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl")

# load new NHC data - update if it's been a while since 10/14/20
# nhc_mega <- StreamPULSE::request_data("NC_NHC", variables=c("WaterTemp_C", "DO_mgL"))
# write_rds(nhc_mega, "data/NHC_watertemp.rds")

nhc_mega <- read_rds("data/NHC_watertemp.rds")
nhc_dat <- nhc_mega$data %>%
  mutate(datetime = with_tz(DateTime_UTC, tzone="EST"),
         date = as.Date(datetime)) %>%
  filter(! flagtype %in% c("Bad Data", "Questionable")) %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  select(date, datetime, water_temp_C = WaterTemp_C, DO_mgL) %>%
  as_tibble() 
  
nhc_19 <- nhc_dat %>%
  group_by(date) %>%
  summarize(water_temp_C = mean(water_temp_C, na.rm=T)) %>%
  ungroup() %>%
  mutate(doy = format(date, "%j")) 

# load Hall NHC data
# this data was digitized from Hall dissertation figure 27
# using WebPlotDigitizer from github

nhc_69 <- read_csv("data/hall/hall_figure27_digitized_mean_daily_temp.csv")
hall_dates <- data.frame(date = seq(nhc_69$date[1], nhc_69$date[nrow(nhc_69)], 
                                    by = 1))
hall_dates$doy <- format(hall_dates$date, "%j")

nhc_69 <- left_join(hall_dates, nhc_69)
nhc_69$water_temp_C <- na.approx(nhc_69$water_temp_C)

mean_temps <- nhc_19 %>% 
  group_by(doy) %>%
  summarize(mean_temp_new = mean(water_temp_C, na.rm=T),
            low_temp_new = min(water_temp_C, na.rm=T),
            high_temp_new = max(water_temp_C, na.rm=T)) %>%
  ungroup()

mean_temps <- nhc_69 %>% 
  group_by(doy) %>%
  summarize(mean_temp_hall = mean(water_temp_C, na.rm=T),
            low_temp_hall = min(water_temp_C, na.rm=T),
            high_temp_hall = max(water_temp_C, na.rm=T)) %>%
  ungroup() %>%
  full_join(mean_temps, by="doy") %>%
  mutate(doy = as.numeric(doy))

# png(width=7, height=4, units='in', type='cairo', res=300,
#     filename='figures/daily_water_temp_comp.png')
# 
#   plot(mean_temps$doy, mean_temps$mean_temp_hall, type = "l", 
#        lwd = 2, col = "darkred", ylim = c(0,29),
#        ylab = "Water Temp C", xlab = "Day of Year",
#        main = "Mean Daily Water Temperatures")
#   
#   lines(mean_temps$doy, mean_temps$mean_temp_new, lwd=2)
#   polygon(x = c(mean_temps$doy, rev(mean_temps$doy)),
#           y = c(mean_temps$low_temp_hall, rev(mean_temps$high_temp_hall)),
#           col = alpha("darkred", 0.3), border = NA)
#   polygon(x = c(mean_temps$doy, rev(mean_temps$doy)),
#           y = c(mean_temps$low_temp_new, rev(mean_temps$high_temp_new)),
#           col = alpha("black", 0.3), border = NA)
#   legend("topleft", c(paste0("2016-2020 temp (n = ", nrow(nhc_19), ")"),
#                       paste0("1968-1970 temp (n = ", nrow(nhc_69), ")")),
#          lty = 1, col = c("black","darkred"), bty = "n", lwd = 2, cex = .8)
# dev.off()



# load Hall rating curve, calculate equation

hall_rc <- read_csv("data/hall/hall_figure5_digitized_ratingcurve.csv")
# Calculate discharge from rating curves
# Q = a * L ^ b

m <- lm(log(hall_rc$discharge_m3s) ~ log(hall_rc$stage_cm))
a <- summary(m)$coefficients[1]
b <- summary(m)$coefficients[2]#Summary of the regression statistics
plot(hall_rc$stage_cm, hall_rc$discharge_m3s)#, log = "y")
lines(seq(10, 90, by = 1), exp(a + b * log(seq(10, 90, by = 1))))

# not good for stages above 80 cm
hall_Q <- read_csv("data/hall/hall_figure26_digitized_dailystage.csv") %>%
  mutate(discharge_m3s = exp(a + b * log(stage_cm))) 

hall_dates <- data.frame(date = seq(hall_Q$date[1], hall_Q$date[nrow(hall_Q)], 
                                    by = 1))
hall_dates$doy <- format(hall_dates$date, "%j")

hall_Q <- left_join(hall_dates, hall_Q)

hall_Q$discharge_m3s <- na.approx(hall_Q$discharge_m3s)
plot(hall_Q$date, hall_Q$discharge_m3s)
abline(h=3)

hall_Q <- hall_Q %>% 
  mutate(doy = format(date, "%j")) %>%
  group_by(doy) %>%
  summarize(discharge_m3s_hall = mean(discharge_m3s, na.rm=T)) %>%
  ungroup()

new_Q <- read_csv("../NHC_2019_metabolism/data/metabolism/processed/NHC.csv") %>%
  mutate(datetime = force_tz(DateTime_EST, tzone = "EST"),
         date = as.Date(datetime)) %>%
  select(datetime, date, level_m, discharge) %>%
  group_by(date) %>% 
  summarize(discharge_m3s = mean(discharge, na.rm=T)) %>%
  ungroup() %>%
  mutate(doy = format(date, "%j")) %>%
  group_by(doy) %>%
  summarize(discharge_m3s_new = mean(discharge_m3s, na.rm=T)) %>%
  ungroup() %>%
  full_join(hall_Q) %>%
  mutate(doy = as.numeric(doy))
  

# png(width=7, height=4, units='in', type='cairo', res=300,
#     filename='figures/daily_discharge_comp.png')
# 
#   plot(new_Q$doy, new_Q$discharge_m3s_hall, type = "l",
#        lwd = 2, col = "darkred", 
#        ylim = c(min(new_Q$discharge_m3s_new, na.rm=T), 
#                 max(new_Q$discharge_m3s_new, na.rm=T)),
#        ylab = "Discharge m3s", xlab = "Day of Year",
#        main = "Mean Daily Discharge", log = "y")
# 
#   lines(new_Q$doy, new_Q$discharge_m3s_new, lwd=2)
#   
#   # polygon(x = c(new_Q$doy, rev(new_Q$doy)),
#   #         y = c(new_Q$low_temp_hall, rev(new_Q$high_temp_hall)),
#   #         col = alpha("darkred", 0.3), border = NA)
#   # polygon(x = c(new_Q$doy, rev(new_Q$doy)),
#   #         y = c(new_Q$low_temp_new, rev(new_Q$high_temp_new)),
#   #         col = alpha("black", 0.3), border = NA)
#   legend("topright", c("2016-2020 Q", "1968-1970 Q"),
#          lty = 1, col = c("black","darkred"), bty = "n", lwd = 2)
# dev.off()


QT_comp <- full_join(new_Q, mean_temps, by = "doy") %>%
  mutate(log_Q_new = log10(discharge_m3s_new),
         log_Q_hall = log10(discharge_m3s_hall))

library(ks)

png(width=7, height=6, units='in', type='cairo', res=300,
    filename="figures/QT_contour_comparison.png")    

kernel <- kde(na.omit(QT_comp[,c("mean_temp_new","log_Q_new")]))
plot(kernel, xlab = "Temperature", ylab = "log(Q)", ylim = c(-3, 1), xlim = c(0,30),
     cont=c(30,60,90), col="grey35", lwd = 2,
     main = "temperature discharge regimes")

kernel_hall <- kde(na.omit(QT_comp[,c("mean_temp_hall","log_Q_hall")]))
par(new=T)
plot(kernel_hall, xlab = "", ylab = "",  ylim = c(-3, 1), xlim = c(0,30),
     cont = c(30,60,90), col = "darkred", lwd=2)
legend("topright", cex = 1.4,
       c("2019-2020","1968-1970"),
       col = c("grey35", "darkred"), lty = 1, lwd = 2, bty = "n")
dev.off()
