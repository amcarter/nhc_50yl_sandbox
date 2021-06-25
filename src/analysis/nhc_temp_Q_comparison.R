# Temperature and discharge comparison 
# between Hall 1972 data and modern NHC data

# a carter
# 2020-10-14


library(tidyverse)
library(lubridate)
library(zoo)
library(ContDataQC)
library(RHydro)

# Water Temperature ####
# load new NHC data - update if it's been a while since 10/14/20
# nhc_mega <- StreamPULSE::request_data("NC_NHC", variables=c("WaterTemp_C", "DO_mgL"))
# write_rds(nhc_mega, "data/NHC_watertemp.rds")
nhc <- readRDS("NHC_2019_metabolism/data/metabolism/compiled/met_preds_direct_calc.rds")$filled %>%
  filter(site == "NHC") 
  

nhc_qt <- nhc %>% 
  group_by(doy) %>% 
  summarize(across(any_of(c("temp.water", "discharge")), 
                   .fns = list(max = ~max(., na.rm = T),
                               min = ~min(., na.rm = T)),
                   .names = '{col}_{fn}')) %>%
  ungroup() %>%
  right_join(nhc, by = "doy") %>%
  filter(year ==2019) %>%
  select(date, doy, year, starts_with(c("temp.water", "discharge"))) %>%
  mutate(discharge_min = ifelse(discharge_min < 0.00675, 0.00675, discharge_min))

hall <- read_csv("hall_50yl/code/data/hall/hall_discharge_temp_daily.csv") %>%
  mutate(doy = as.numeric(format(date, '%j')),
         year = case_when(date < as.Date("1969-04-08") ~ 1968, 
                          date < as.Date("1970-04-08") ~ 1969,
                          TRUE ~ 1970)) %>%
  rename(temp.water = water_temp_C, discharge = discharge_m3s)

hall_qt <- hall %>%
  group_by(doy) %>% 
  summarize(across(any_of(c("temp.water", "discharge")), 
                   .fns = list(max = ~max(., na.rm = T),
                               min = ~min(., na.rm = T)),
                   .names = '{col}_{fn}')) %>%
  ungroup() %>%
  right_join(hall, by = "doy") %>%
  filter(year == 1969) %>%
  select(date, doy, year, starts_with(c("temp.water", "discharge"))) %>%
  mutate(discharge_min = ifelse(discharge_min < 0.00675, 0.00675, discharge_min))

h <- hall %>% group_by(month(date)) %>% summarize_all(mean, na.rm = T)
nhc %>% group_by(month(date)) %>% summarize_all(mean, na.rm = T) %>%
  select('month(date)', temp_now = temp.water) %>%
  full_join(h, by = 'month(date)') %>%
  mutate(temp_diff = temp_now - temp.water)

png("figures/watertemp_comparison.png", width = 6, height = 4.5, 
    units = 'in', res = 300)
  par(mfrow = c(2,1), mar = c(0,4,.5,2), oma = c(4,0,1,0))
  plot(nhc_qt$doy, nhc_qt$temp.water, type = 'l', lwd = 2, ylim = c(0, 28),
       xaxt = 'n', xlab = 'date', ylab = '')
  mtext('Water temperature C', 2, 2.1)
  polygon(c(nhc_qt$doy, rev(nhc_qt$doy)), 
          c(nhc_qt$temp.water_min, rev(nhc_qt$temp.water_max)), 
          col = alpha(now_col, .5),  border = NA)
  polygon(c(hall_qt$doy, rev(hall_qt$doy)), 
          c(hall_qt$temp.water_min, rev(hall_qt$temp.water_max)),
          col = alpha(then_col, .3), border = NA)
  lines(nhc_qt$doy, nhc_qt$temp.water, lwd = 2 )
  lines(hall_qt$doy, hall_qt$temp.water, col = then_col,lwd = 2 )
  legend('topright', legend = c(2019, 1969), col = c(1, then_col), lwd = 2,
         lty = 1, ncol = 2, bty = 'n', cex = .8)
  plot(nhc_qt$doy, nhc_qt$discharge, log = 'y', xaxt = 'n', xlab = 'Date', 
       lwd = 2, type = 'l', yaxt = 'n', ylab = '')
  mtext('Discharge (m3/s)',2,2.1) 
  polygon(c(nhc_qt$doy, rev(nhc_qt$doy)),
          c(nhc_qt$discharge_min, rev(nhc_qt$discharge_max)), 
          col = alpha(now_col, .5), border = NA)
  polygon(c(hall_qt$doy, rev(hall_qt$doy)),
          c(hall_qt$discharge_min, rev(hall_qt$discharge_max)),
          col = alpha(then_col, .3), border = NA)
  lines(nhc_qt$doy, nhc_qt$discharge, lwd = 2)
  lines(hall_qt$doy, hall_qt$discharge, lwd = 2, col = then_col)
  axis(1, at = as.numeric(format(seq(as.Date('2020-01-01'), 
                                     by = 'month', length.out = 12),'%j')),
       labels = month.abb) 
  axis(2, at = c(.01, 1, 100), labels = c(0.01, 1, 100))
  
dev.off()

# RBI calculations ####

Q_stats <- nhc %>%
  arrange(date) %>%
  mutate(discharge = na.approx(discharge, na.rm = F),
         bf = RHydro::baseflow_sep(discharge),
         stormflow = discharge - bf,
         bf_frac = bf/discharge) %>%
  group_by(year) %>%
  summarize(RBI = RBIcalc(discharge),
            ar_1 = arima(discharge, order = c(1,0,0))$coef[1],
            q05 = quantile(discharge, .05, na.rm = T),
            mean = mean(discharge, na.rm = T),
            median = median(discharge, na.rm = T), 
            bf_frac = mean(bf_frac, na.rm = T))

write_csv(Q_stats, 'NHC_2019_metabolism/data/rating_curves/annual_Q_stats.csv')
bf_19 <- nhc %>%

  select(year, bf_frac) %>%
  group_by(year)
bf_69 <- Q %>%
  filter(year==1969,
         !is.na(discharge_m3s)) %>%
  mutate(discharge_m3s = na.approx(discharge_m3s, na.rm = F),
         bf = RHydro::baseflow_sep(discharge_m3s),
         stormflow = discharge_m3s - bf,
         bf_frac = bf/discharge_m3s) 


summary(bf_69)


Q_stats <- data.frame(bffrac = c(sum(bf_69$bf)/sum(bf_69$discharge_m3s),
                                 sum(bf_19$bf)/sum(bf_19$discharge_m3s))) %>%
  bind_cols(Q_stats)
write_csv(Q_stats, "../NHC_2019_metabolism/data/rating_curves/Q_stats_nowthen_hall.csv")

# old ####


# png("figures/discharge_comparison.png", width = 5, height = 4, 
#     units = 'in', res = 300, family = 'cairo')
# ggplot(nhc_qt, aes(doy, log(discharge))) +
#   geom_line(lwd = 1.2) +
#   geom_ribbon(aes(ymin = log(discharge_min), ymax = log(discharge_max)),
#               fill = alpha("grey40", .2)) +
#   geom_line(data = hall_qt, aes(doy, log(discharge)), col = "brown4", lwd = 1.2) +
#   geom_ribbon(data = hall_qt, aes(ymin = log(discharge_min),
#                                   ymax = log(discharge_max)),
#               fill = alpha("brown3", .2)) +
#   theme_bw() +
#   ylab("log(discharge)") +
#   xlab("Day of year") + 
#   ylim(-5,10) 
# dev.off()

nhc_19 <- nhc_mega %>%
  group_by(date = as.Date(with_tz(DateTime_UTC, tz = "EST"))) %>%
  summarize(water_temp_C = mean(temp.water, na.rm=T),
            discharge_m3s = mean(discharge, na.rm = T)) %>%
  ungroup() %>%
  mutate(doy = format(date, "%j")) 

# load Hall NHC data
# this data was digitized from Hall dissertation figure 27
# using WebPlotDigitizer from github

nhc_69 <- read_csv("hall_50yl/code/data/hall/hall_discharge_temp_daily.csv") %>%
  mutate(doy = format(date, "%j"),
         hall = "then") %>%
  select(-stage_cm, -notes)


mean_temps <- nhc_19 %>% 
  mutate(hall = "now") %>%
  bind_rows(nhc_69) %>%
  group_by(doy, hall) %>%
  summarize(across(-date, .fns = list(mean = ~mean(., na.rm = T),
                                      min = ~min(., na.rm = T),
                                      max = ~max(., na.rm = T),
                                      sd = ~sd(., na.rm = T)),
                   .names = "{col}_{fn}")) %>%
  ungroup() %>%
  mutate(doy = as.numeric(doy))

Q <- nhc_19 %>% 
  mutate(hall = "now") %>%
  bind_rows(nhc_69) %>%
  mutate(year = year(date), 
         doy = as.numeric(doy)) %>%
  filter(year %in% c(1969, 2019)) %>%
  left_join(mean_temps, by = c("hall", "doy")) %>%
  select(doy, year, discharge_m3s, discharge_m3s_sd)  


# Water Temp Plot ####

png(width=7, height=4, units='in', type='cairo', res=300,
    filename='figures/daily_water_temp_comp.png')

  plot(mean_temps$doy[mean_temps$hall == "then"], 
       na.approx(mean_temps$water_temp_C_mean[mean_temps$hall == "then"], 
                 na.rm = F),
       type = "l", lwd = 2, 
       col = "brown3", ylim = c(0,29),
       ylab = "Water Temp C", xlab = "Day of Year",
       main = "Mean Daily Water Temperatures")
  # points(mean_temps$doy, mean_temps$mean_temp_hall, pch = 20)
  polygon(x = c(mean_temps$doy[mean_temps$hall == "then"],
                rev(mean_temps$doy[mean_temps$hall == "then"])),
          y = na.approx(c(mean_temps$water_temp_C_min[mean_temps$hall == "then"],
                          rev(mean_temps$water_temp_C_max[mean_temps$hall == "then"])),
                          na.rm = F),
          col = alpha("brown3", 0.5), border = NA)
  lines(mean_temps$doy[mean_temps$hall == "now"], 
        mean_temps$water_temp_C_mean[mean_temps$hall == "now"], 
        lwd=2, col = "grey40")
  polygon(x = c(mean_temps$doy[mean_temps$hall == "now"],
                rev(mean_temps$doy[mean_temps$hall == "now"])),
          y = c(mean_temps$water_temp_C_min[mean_temps$hall == "now"],
                rev(mean_temps$water_temp_C_max[mean_temps$hall == "now"])),
          col = alpha("grey40", 0.5), border = NA)
  legend("topleft", c(paste0("2016-2020 temp (n = ", nrow(nhc_19), ")"),
                      paste0("1968-1970 temp (n = ", 
                             length(which(!is.na(nhc_69$water_temp_C))), ")")),
         lty = 1, col = c("grey40", "brown3"), bty = "n", lwd = 2)
dev.off()



# Monthly summaries ####
t19 <- nhc_19 %>%
  mutate(month = format(date, "%b")) %>%
  filter(date >= "2017-03-01",
         date <= "2020-03-01") %>%
  group_by(month) %>%
  summarize(mean_temp = mean(water_temp_C, na.rm = T),
            min_temp = min(water_temp_C, na.rm = T),
            max_temp = max(water_temp_C, na.rm = T))

t19 <- nhc_69 %>%
  mutate(month = format(date, "%b")) %>%
  group_by(month) %>%
  summarize(mean_temph = mean(water_temp_C, na.rm = T),
            min_temph = min(water_temp_C, na.rm = T),           
            max_temph = max(water_temp_C, na.rm = T))%>%
  full_join(t19, by = "month") %>%
  mutate(deltat = mean_temp - mean_temph)

maxh <- 25.42298
nhc_19 %>%
  mutate(month = format(date, "%b"),
         above = ifelse(water_temp_C >= maxh, T, F)) %>%
  filter(date >= "2017-03-01",
         date <= "2020-03-01") %>%
  group_by(month) %>%
  summarize(high = sum(above, na.rm = T),
            all = sum(!is.na(above)))
  
# Discharge plot ####

plot(Q$doy[Q$year == 2019], Q$discharge_m3s[Q$year == 2019],
     col = "grey40", lwd = 2, type = "l", ylab = "Discharge (m3s)",
     xlab = "Day of Year", log = "y") 
lines(Q$doy[Q$year == 1969], Q$discharge_m3s[Q$year == 1969],
     col = "brown3", lwd = 2) 

# RBI calculations ####

Q_stats <- Q %>%
  group_by(year) %>%
  summarize(RBI = RBIcalc(discharge_m3s),
            ar_1 = arima(discharge_m3s, order = c(1,0,0))$coef[1],
            q05 = quantile(discharge_m3s, .05, na.rm = T),
            mean = mean(discharge_m3s, na.rm = T),
            median = median(discharge_m3s, na.rm = T))


bf_19 <- Q %>%
  filter(year==2019,
         !is.na(discharge_m3s)) %>%
  mutate(discharge_m3s = na.approx(discharge_m3s, na.rm = F),
         bf = RHydro::baseflow_sep(discharge_m3s),
         stormflow = discharge_m3s - bf,
         bf_frac = bf/discharge_m3s)
bf_69 <- Q %>%
  filter(year==1969,
         !is.na(discharge_m3s)) %>%
  mutate(discharge_m3s = na.approx(discharge_m3s, na.rm = F),
         bf = RHydro::baseflow_sep(discharge_m3s),
         stormflow = discharge_m3s - bf,
         bf_frac = bf/discharge_m3s) 


summary(bf_69)


Q_stats <- data.frame(bffrac = c(sum(bf_69$bf)/sum(bf_69$discharge_m3s),
                      sum(bf_19$bf)/sum(bf_19$discharge_m3s))) %>%
  bind_cols(Q_stats)
write_csv(Q_stats, "../NHC_2019_metabolism/data/rating_curves/Q_stats_nowthen_hall.csv")
# Discharge and Level ####

# load Hall rating curve, calculate equation
# 
# hall_rc <- read_csv("data/hall/hall_figure5_digitized_ratingcurve.csv")
# # Calculate discharge from rating curves
# # Q = a * L ^ b
# 
# m <- lm(log(hall_rc$discharge_m3s) ~ log(hall_rc$stage_cm))
# a <- summary(m)$coefficients[1]
# b <- summary(m)$coefficients[2]#Summary of the regression statistics
# plot(hall_rc$stage_cm, hall_rc$discharge_m3s)#, log = "y")
# lines(seq(10, 90, by = 1), exp(a + b * log(seq(10, 90, by = 1))))
# min(hall_rc$stage_cm)
# 
# hall_1cms_stage <- exp((log(1) - a)/b)/100 
# 
# # not good for stages above 80 cm
# hall_Q <- read_csv("data/hall/hall_figure26_digitized_dailystage.csv") %>%
#   mutate(discharge_m3s = exp(a + b * log(stage_cm)),
#          date = as.Date(date, format = "%m/%d/%Y")) %>%
#   arrange(date)
# 
# hall_dates <- data.frame(date = seq(hall_Q$date[1], hall_Q$date[nrow(hall_Q)], 
#                                     by = 1))
# hall_dates$doy <- as.numeric(format(hall_dates$date, "%j"))
# 
# hall_Q <- left_join(hall_dates, hall_Q)
# 
# hall_Q$discharge_m3s <- na.approx(hall_Q$discharge_m3s)
# hall_Q$stage_m <- na.approx(hall_Q$stage_cm)/100
# # plot(hall_Q$date, hall_Q$discharge_m3s)
# # abline(h=3)
# 
# hall_Qlim <- hall_Q %>% 
#   group_by(doy) %>%
#   summarize(low_discharge_m3s = min(discharge_m3s, na.rm=T),
#             high_discharge_m3s = max(discharge_m3s, na.rm = T), 
#             low_stage_m = min(stage_m, na.rm = T),
#             high_stage_m = max(stage_m, na.rm = T))
# 
# hall_Q <- hall_Q %>%
# #  filter(year(date) == 1969) %>%
#   mutate(year = year(date)) %>%
#   left_join(hall_Qlim, by = "doy") %>%
#   select(-stage_cm) %>%
#   as_tibble()
# 
# new_Q <- read_csv("../NHC_2019_metabolism/data/metabolism/processed/NHC.csv") %>%
#   mutate(datetime = force_tz(DateTime_EST, tzone = "EST"),
#          date = as.Date(datetime)) %>%
#   select(datetime, date, level_m, discharge) %>%
#   group_by(date) %>% 
#   summarize(discharge_m3s = mean(discharge, na.rm=T),
#             stage_m = mean(level_m, na.rm = T)) %>%
#   mutate(doy = as.numeric(format(date, "%j")),
#          discharge_m3s = ifelse(discharge_m3s < 1e-4, 
#                             NA, discharge_m3s))
# 
#   
# m <- lm(log(new_Q$discharge_m3s) ~ log(new_Q$stage_m))
# a_n <- summary(m)$coefficients[1]
# b_n <- summary(m)$coefficients[2]#Summary of the regression statistics
# plot(new_Q$stage_m, new_Q$discharge_m3s, ylim = c(0,10))
# lines(seq(0.5, 1.50, by = .1), exp(a_n + b_n * log(seq(.5, 1.5, by = .1))))
# 
# # find what flow is at NHC when flow at CBP is 1 m3s
# # qq <- read_csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/siteData/interpolatedQ_allsites.csv") %>%
# #   filter(CBP.Q <=3)
# # plot(qq$CBP.Q, qq$NHC.Q)
# # abline(0,1, col = 2, lty = 2, lwd = 1)
# # qm <- lm(qq$NHC.Q ~ qq$CBP.Q)
# # qs = summary(qm)$coefficients[2]
# # abline(0,qs, col = 3, lty = 2, lwd = 2)
# 
# qs <- 1.06683
# 
# new_1cms_stage <- exp((log(qs) - a_n)/b_n) 
# offset <- new_1cms_stage - hall_1cms_stage
# 
# new_Q <- new_Q %>%
#   mutate(mod_stage = exp((log(discharge_m3s) - a_n)/ b_n),
#          stage_m = ifelse(is.na(stage_m), 
#                               mod_stage, 
#                               stage_m)) %>%
#   select(-mod_stage)
# # don't let offest be larger than min stage:
# offset <- min(c(offset, new_Q$stage_m), na.rm = T) 
# 
# new_Q$stage_m <- new_Q$stage_m - offset
# 
# qlim <- new_Q %>%
#   group_by(doy) %>%
#   summarize(low_discharge_m3s = min(discharge_m3s, na.rm = T), 
#             high_discharge_m3s = max(discharge_m3s, na.rm = T), 
#             low_stage_m = min(stage_m, na.rm = T),
#             high_stage_m = max(stage_m, na.rm = T))
# 
# Q <- new_Q %>%
# #  filter(year(date) == 2019) %>%
#   mutate(year = year(date)) %>%
# #  select(-date) %>%
#   left_join(qlim, by = "doy") %>%
#   bind_rows(hall_Q) %>%
#   arrange(date)

# plots of then and now Q  ####



png(width=7, height=4, units='in', type='cairo', res=300,
    filename='figures/daily_level_focalyears.png')

  plot(Q$doy, Q$stage_m,
       type = 'n',
       ylab = "level (m)", xlab = "Day of Year",
       main = "Mean Daily Level")#, log = "y")
  for(y in 2016:2020){
    lines(Q$doy[Q$year == y], Q$stage_m[Q$year == y],
          lwd = 1.5, col = alpha("black", 1))
    }
  lines(Q$doy[Q$year == 2019], Q$stage_m[Q$year == 2019], lwd = 2)
  for(y in 1968:1970){
    lines(Q$doy[Q$year == y], Q$stage_m[Q$year == y],
          lwd = 1.5, col = alpha("brown3", 1))
    }
  lines(Q$doy[Q$year == 1969], Q$stage_m[Q$year == 1969],
        lwd = 2, col = "brown3")
  legend("topright", 
         c("2019", "1969"),
         lty = 1, col = c("black", "brown3"),
         bty = "n", lwd = 2, ncol = 2)
dev.off()

Q$discharge_m3s[Q$discharge_m3s < 1e-3] <- 1e-3
png(width=7, height=4, units='in', type='cairo', res=300,
    filename='figures/daily_discharge_comp.png')
plot(Q$doy, Q$discharge_m3s,
     type = 'n', log = "y",
     ylab = "discharge (m3s)", xlab = "Day of Year",
     main = "Mean Daily Q")#, log = "y")
for(y in 2016:2020){
  lines(Q$doy[Q$year == y], Q$discharge_m3s[Q$year == y],
        lwd = 2, col = alpha("steelblue", .5))
}
lines(Q$doy[Q$year == 2019], Q$discharge_m3s[Q$year == 2019], lwd = 2, col = "steelblue")
for(y in 1968:1970){
  lines(Q$doy[Q$year == y], Q$discharge_m3s[Q$year == y],
        lwd = 2, col = alpha("black", .5))
}  
lines(Q$doy[Q$year == 1969], Q$discharge_m3s[Q$year == 1969],
      lwd = 2)#, col = "grey40")
  legend("topright", 
         c("2016-2019", "1968-1970"),
         lty = 1, col = c("steelblue", "black"),
         bty = "n", lwd = 2, ncol = 2)


dev.off()
plot(Q$doy, na.approx(Q$discharge_m3s, na.rm = F), type = "l",
       lwd = 2, col = "grey50",
       ylim = c(min(Q$low_discharge_new, na.rm=T),
                max(Q$high_discharge_new, na.rm=T)),
       ylab = "discharge m3s", xlab = "Day of Year",
       main = "Mean Daily Discharge", log = "y")
  polygon(c(Q$doy, rev(Q$doy)), 
          c(Q$low_discharge_m3s_hall, rev(Q$high_discharge_m3s_hall)),
          col = alpha("grey40", .5), border = NA)
  lines(Q$doy, Q$discharge_m3s_new, col = "steelblue", lwd=2)
  polygon(c(Q$doy, rev(Q$doy)), 
          c(Q$low_discharge_new, rev(Q$high_discharge_new)),
          col = alpha("steelblue", .5), border = NA)

  legend("topright", 
         c("2019 n = 365", 
           paste0("1969 n = ", length(which(!is.na(Q$discharge_m3s))))),
         lty = 1, pch = c(NA, 20),
         col = c("steelblue","grey50"),
         bty = "n", lwd = 2)

dev.off()

# QT kernel density plots ####

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


