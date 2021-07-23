# 1. Correct digitized temperature and discharge from Hall 1970 using 
#    peak flow dates to shift dates
# 2. Pair daily stage, discharge, and temperature with depths and K600 from
#    days with modeled metabolism in Hall 1970
# 3. plot distributions of discharge, depth, and K600 for CBP in 1969 and today


library(tidyverse)
library(lubridate)
library(zoo)
# load Hall data for CBP site:

hall <- read_csv("C:/Users/Alice Carter/git/nhc_50yl/hall_50yl/code/data/hall/hall_discharge_temp_daily.csv") %>%
  mutate(doy = as.numeric(format(date, '%j')),
         year = case_when(date < as.Date("1969-04-08") ~ 1968, 
                          date < as.Date("1970-04-08") ~ 1969,
                          TRUE ~ 1970)) %>%
  rename(temp.water = water_temp_C, discharge = discharge_m3s)
hall_met <- read_csv("C:/Users/Alice Carter/git/nhc_50yl/hall_50yl/code/data/hall/hall_table_15.csv") %>%
  filter(site == 'Concrete') %>%
  mutate(date = as.Date(date, format = '%m/%d/%Y'))
hall_dat <- hall_met %>% full_join(hall, by = 'date') %>%
  select(date, depth_m, discharge, level_m) %>%
  arrange(date)

# look at depth x stage and discharge relationships:
plot(hall_dat$depth_m, hall_dat$level_m)
abline(0,1)
plot(hall_dat$depth_m, log(hall_dat$discharge))
m <- lm(log(discharge) ~ depth_m, data = hall_dat)
lines(m)


# compare the depth and stage measurements:
hall_stages <- read_csv('C:/Users/Alice Carter/git/nhc_50yl/hall_50yl/code/data/hall/hall_table_13_p.csv') %>%
  mutate(level_meas = water_stage_cm/100) %>%
  filter(site == 'Concrete') %>%
  select(date, level_meas)
hall_stages <- 
  read_csv('C:/Users/Alice Carter/git/nhc_50yl/hall_50yl/code/data/hall/hall_table_5_flood_stages_cbp.csv') %>%
  mutate(level_flood = stage_cm/100) %>%
  select(-stage_cm) %>%
  bind_rows(hall_stages)
d <- hall_dat %>% 
  full_join(hall_stages) %>%
  arrange(date)
plot(d$date, d$level_m, type = 'l')
points(d$date, d$level_meas, col = 2)
points(d$date, d$level_flood, col = 3)

# move the stage up to match discharge, and shift to minimize gap?
# d <- hall_dat %>%
#   mutate(level_shifted = level_m + mean(depth_m - level_m, na.rm = T))
# plot(d$date, d$level_shifted, type = 'l')
# points(d$date, d$depth_m, col = 2)
# test different offsets for depth and level time series:
n <- nrow(d)
stage <- d$level_m[1:(n-20)]
offset <- data.frame(shift = seq(1:20),
                     diff = rep(NA_real_, 20))
for(i in 1:20){
  offset$diff[i] <- mean(abs(d$level_meas[i:(n - 21 + i)] - stage), na.rm = T) 
  plot(d$date[i:(n - 21+i)], d$level_m[1:(n - 20)], type = 'l', main = i)
  points(d$date[i:(n - 21+i)], d$level_meas[i:(n-21+i)], col = 2)
  points(d$date[i:(n - 21+i)], d$level_flood[i:(n - 21 + i)], col = 3)
  }

# a 13 day offset seems to do the best at lining up the storms where they should be:
hall$date <- hall$date + 13


hall_dat <- hall_met %>% 
  full_join(hall, by = 'date') %>%
  select(date, depth_m, discharge, level_m, temp.water) %>%
  arrange(date)
write_csv(hall_dat, 'C:/Users/Alice Carter/git/nhc_50yl/hall_50yl/code/data/hall_discharge_temp_daily_corrected_dates.csv')
# add K values ####
hall_k <- read_csv("C:/Users/Alice Carter/git/nhc_50yl/hall_50yl/code/data/hall/hall_tableA2_k_morphology.csv")

##Convert KO2 to K600
K600fromO2<-function(temp, KO2) {
  ((600/(1800.6-(120.1*temp)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5)*KO2
}

hall_kd <- hall_k %>%
  mutate(K_O2 = 2.3 * 9 * K2_perday_perarea/24,
         K600 = K600fromO2(20, K_O2*24)) %>%  
  select(depth_m, K600)
hall_kd <- bind_rows(hall_dat, hall_kd) %>%
  arrange(depth_m)
hall_kd$K600[1] <- 1
hall_kd$K600 <- zoo::na.approx(hall_kd$K600, x = hall_kd$depth_m, na.rm = F)
hall_dat <- hall_kd %>% 
  filter(!is.na(date)) %>%
  arrange(date)

plot(hall_dat$discharge, hall_dat$K600, log = 'x', xlim = c(1e-2, 10))


# load data for CBP site: ####
cbp <- read_csv("C:/Users/Alice Carter/git/nhc_50yl/NHC_2019_metabolism/data/metabolism/processed/CBP.csv")
met <- read_csv("C:/Users/Alice Carter/git/nhc_50yl/NHC_2019_metabolism/data/metabolism/compiled/daily_preds_stream_metabolizer.csv") %>%
  filter(site == "CBP") %>%
  select(date, K600)
today <- cbp %>%
  mutate(date = as.Date(DateTime_UTC, tz = 'EST')) %>%
  left_join(met, by = 'date') %>%
  group_by(date, site) %>%
  select(-DateTime_UTC) %>%
  summarize_all(mean, na.rm = T) %>%
  ungroup()



png("C:/Users/Alice Carter/git/nhc_50yl/figures/disch_depth_K600_distributions_CBP.png", 
    width = 6, height = 4, units = 'in', res = 300)

  par(mfrow = c(1,3), 
      mar = c(4, 1, 3, 1), 
      oma = c(0,3, 0, 1))
  
  plot(density(log(today$discharge), na.rm = T), ylim = c(0, 0.45), 
       xlim = c(-5, 5), main = "Discharge", xlab = "log discharge (m3s)")
  lines(density(log(hall_dat$discharge), na.rm = T), lty = 2)
  mtext("Density", 2, line = 2.3, cex = .8)
  
  plot(density(today$depth, na.rm = T), xlim = c(0, 1.5), main = "Depth", 
       xlab = "average depth (m)", yaxt = 'n')
  lines(density(hall_dat$depth_m, na.rm = T), lty = 2)
  legend('topright', legend = c('2019', '1969'), bty = 'n', lty = c(1,2),
         cex = 1.4)
  
  plot(density(today$K600, na.rm = T), main = "Gas Exchange", 
       xlab = "K600 (d-1)", yaxt = 'n')
  lines(density(hall_dat$K600, na.rm = T), lty = 2)

dev.off()



