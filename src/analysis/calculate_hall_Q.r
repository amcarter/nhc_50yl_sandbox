# calculate daily discharge for the CBP from Hall 1970
# Use digitized plots for rating curve (figure 5) and daily stage (figure 27)
# Pair with temperature daily data and snap so both time series start on the same date

library(tidyverse)
stage <- read_csv("C:/Users/Alice Carter/git/nhc_50yl/hall_50yl/code/data/hall/hall_figure26_digitized_dailystage.csv") 
hall_rc <- read_csv("C:/Users/Alice Carter/git/nhc_50yl/hall_50yl/code/data/hall/hall_figure5_digitized_ratingcurve.csv")

# Calculate discharge from rating curves
# Q = a * L ^ b
 
m <- lm(log(hall_rc$discharge_m3s) ~ log(hall_rc$stage_cm))
a <- summary(m)$coefficients[1]
b <- summary(m)$coefficients[2]#Summary of the regression statistics
plot(hall_rc$stage_cm, hall_rc$discharge_m3s)#, log = "y")
lines(seq(10, 90, by = 1), exp(a + b * log(seq(10, 90, by = 1))))
rc_range <- range(hall_rc$stage_cm)/100


hall_Q <- stage %>%
  mutate(discharge_m3s = exp(a + b * log(stage_cm))) %>%
  arrange(date)

write_csv(hall_Q, "C:/Users/Alice Carter/git/nhc_50yl/hall_50yl/code/data/hall/hall_figure26_digitized_dailystage.csv")
# calculate stage at 1 cms discharge:
hall_1cms_stage <- exp((log(1) - a)/b)/100

# pair with temperature:
temp <- read_csv('C:/Users/Alice Carter/git/nhc_50yl/hall_50yl/code/data/hall/hall_figure27_digitized_mean_daily_temp.csv')
temp$date <- temp$date - 5
hall <- full_join(temp, hall_Q, by = 'date') %>%
  mutate(level_m = stage_cm/100) %>%
  select(-stage_cm) %>%
  arrange(date)

dates <- seq(hall$date[1], hall$date[nrow(hall)], by = 'day')
hall <- data.frame(date = dates) %>%
  left_join(hall, by = 'date') %>%
  mutate(across(c(-date, -notes), zoo::na.approx, na.rm = F)) %>%
  as_tibble()

write_csv(hall, "C:/Users/Alice Carter/git/nhc_50yl/hall_50yl/code/data/hall/hall_discharge_temp_daily.csv")
