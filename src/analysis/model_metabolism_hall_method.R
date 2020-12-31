# Calculate Metabolism using Hall method

library(tidyverse)
library(lubridate)
library(zoo)
library(pracma)
setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")
source("../src/streamMetabolizer/inspect_model_fits.r")


kq_hall <- read_csv("siteData/KQ_hall_prior_from_equation_daily.csv")
flow_dates <- read_csv("rating_curves/flow_dates_filter.csv")

# Metabolism functions ####
# calculate hall met given k2 and a day of data spaced at dt min intervals
fit_metab_hall <- function(day, k2, dt = 15, plot = TRUE){
  # calculate the slope (in mgo2/m3/hr) at each point
  d <- as.Date(day$DateTime_EST[1])
  day$dodt <- c(day$DO.obs[2] - day$DO.obs[1], diff(day$DO.obs)) * 60/dt 
  
  day <- day %>%
    mutate(k_gm3hr = (2.3*DO.sat*k2)/(24),          # k in g/m3/h/100% sat def
           dodt_diff = (DO.sat - DO.obs) * k_gm3hr, # diffusion due to sat def
           dodt_cor = dodt - dodt_diff,             # dodt corrected for diff
           hour = (seq(1:nrow(day))-1)/60 * dt)
  # if missing data, cannot calc met (these ts are already gap filled up to 3 hrs)
  if(sum(is.na(day$dodt_cor)) > 0) { 
    m <- data.frame(date = d)
    print(paste(d, "missing data"))
    return(m)
  }
  
  # find sunrise and sunset based on light
  sun <- which(day$light > 0)
  sunrise <- sun[1]-1
  sunset <- sun[length(sun)]+1
  # find minimum change in two hours post sunset
  ps_min <- sunset -1 + which.min(day$dodt_cor[sunset:(sunset+8)])
  
  # build polygon for gpp by connecting sunrise to post sunset minimum
  a_gpp <- day %>%
    slice((sunrise):(ps_min)) %>%
    select(hour, dodt_cor)
  a_gpp$line = NA_real_
  a_gpp$line[c(1, nrow(a_gpp))] <- a_gpp$dodt_cor[c(1, nrow(a_gpp))]
  a_gpp <- a_gpp %>%
    mutate(line = na.approx(line),
           dodt_cor = ifelse(dodt_cor < line, line, dodt_cor)) %>%
    select(-line)
  
  # build polygon for er that is all negative values of DO/dt
  a_er <- day %>%
    select(hour, dodt_cor) %>%
    mutate(dodt_cor = ifelse(dodt_cor > 0, 0, dodt_cor))
  start <- data.frame(hour = 0, dodt_cor = 0)
  end <- data.frame(hour = 24, dodt_cor = 0)
  a_er <- bind_rows(start, a_er, end)

  if(plot) {
  # par(mfrow = c(3,1))
  # plot(day$hour, day$DO.obs, type = "b")
  # plot(day$hour, day$DO.obs/day$DO.sat, type = "l")
  # abline(h = 1)
  plot(day$hour, day$dodt_cor, type = "l",lty = 2, main = d, lwd = 2, 
       ylim = range(c(day$dodt, day$dodt_cor), na.rm = T),
       ylab = "dDO/dt (gO2/m3/hr)", xlab = "hour")
  lines(day$hour, day$dodt, lwd = 2)
  abline(h = 0)
  abline(v = day$hour[c(sunrise, sunset)], lwd = 2, 
         lty = 3, col = "goldenrod")
  # lines(day$hour, day$dodt_diff, lty = 2)

  polygon(a_gpp$hour, a_gpp$dodt_cor, border = NA,
          col = alpha("forestgreen", 0.3))
  polygon(a_er$hour, a_er$dodt_cor, border = NA,
          col = alpha("sienna", 0.3))
  legend("topright",
         legend = c("dDO/dt uncorrected", "dDO/dt corrected",
                    "GPP", "ER", "sunrise/set"), 
         fill = c(NA, NA, alpha("forestgreen", 0.4), alpha("sienna", 0.4), NA),
         col = c(1, 1, NA, NA, "goldenrod"),
         border = NA, 
         lty = c(1, 2, NA, NA, 3),
         bty = "n")
  }
  
  m <- data.frame(date = d,
                  gpp_gO2m3d = abs(polyarea(a_gpp$hour, a_gpp$dodt_cor)),
                  er_gO2m3d = abs(polyarea(a_er$hour, a_er$dodt_cor)),
                  K2_day = k2)
  
  return(m)
}
fit_ts_metab_hall<- function(dat, kq_hall, s, plot = FALSE){
  d1 <- as.Date(dat$DateTime_EST[1]) + 1
  d2 <- as.Date(dat$DateTime_EST[nrow(dat)]) - 1
  
  dates <- seq(d1, d2, by = "day")
  met <- data.frame()
  for(i in 1:length(dates)){
    d <- dates[i]
    ds <- ymd_hms(paste(d, "00:00:00"), tz = "EST")
    day <- dat %>%
      filter(DateTime_EST >= ds, DateTime_EST <= ds + 60*60*24)
    
    kq <- kq_hall %>%
      filter(date == d,
             site == s)
    k2 <- kq$k2[1]
    
    m <- fit_metab_hall(day, k2, plot = plot)
    met <- bind_rows(met, m)
  }
  
  met <- dat %>%
    mutate(date = as.Date(DateTime_EST)) %>%
    group_by(date) %>%
    summarize(discharge = mean(discharge, na.rm = T),
              depth = mean(depth, na.rm = T),
              temp.water = mean(temp.water, na.rm = T)) %>%
    ungroup() %>%
    full_join(met) %>%
    mutate(site = s,
           gpp_gO2m2d = gpp_gO2m3d * depth, 
           er_gO2m2d = er_gO2m3d * depth) %>%
    arrange(date)
  
  return(met)
}
K600fromO2 <- function(temp, KO2) {
  ((600/(1800.6-(120.1*temp)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5)*KO2
}
filter_hall_met <- function(met, flow_dates){
  preds <- met %>%
    mutate(K600 = K600fromO2(temp.water, K2_day)) %>%
    select(date, discharge, depth, temp.water,
           GPP = gpp_gO2m2d, ER = er_gO2m2d, K600) %>%
    left_join(flow_dates[,c(1,4)]) %>%
    mutate(across(all_of(c("GPP", "ER", "K600")), 
                  ~ ifelse(good_flow, ., NA)))
  w <- range(which(!is.na(preds$GPP)), na.rm = T)
  preds <- preds[w[1]:w[2],]
  
  coverage <- data.frame(missing_data = sum(is.na(met$gpp_gO2m2d)),
                         bad_flow = sum(preds$good_flow == FALSE &
                                          !is.na(met$gpp_gO2m2d), na.rm = T),
                         total_days = nrow(preds))
  sump <- preds %>%
    summarize(gpp_mean = mean(GPP, na.rm = T),
              gpp_median = median(GPP, na.rm = T),
              gpp_max = max(GPP, na.rm = T),
              gpp_min = min(GPP, na.rm = T),
              er_mean = mean(ER, na.rm = T),
              er_median = median(ER, na.rm = T),
              er_max = max(ER, na.rm = T),
              er_min = min(ER, na.rm = T))
  
  cum <- data.frame(date = seq(preds$date[1], 
                               preds$date[nrow(preds)], 
                               by = "day")) %>%
        as_tibble() %>%
        left_join(preds) %>%
        select(-K600, -good_flow, -discharge, -depth) %>%
        mutate(across(-date, na.approx, na.rm = F)) %>%
        mutate(across(-date, cumsum, .names = "{col}_cum")) 

  n <- nrow(cum)
  l = (n-9)
  weekly <- tibble(date = cum$date[1:l],
                   GPP_week = rep(NA_real_, l),
                   ER_week = rep(NA_real_, l))
  for(i in 1:l){
    weekly$GPP_week[i] <- sum(cum$GPP[i:(i+9)]) 
    weekly$ER_week[i] <- sum(cum$ER[i:(i+9)]) 
  }
  sump$gpp_max10d <- weekly$date[which.max(weekly$GPP_week)]
  sump$er_max10d <- weekly$date[which.max(weekly$ER_week)]
  sump$gpp_cum <- cum$GPP_cum[n]*365/n
  sump$er_cum <- cum$ER_cum[n]*365/n
  sump$daterange <- as.character(paste(cum$date[1], "-", cum$date[nrow(cum)]))
  sump$pctcoverage <- sum(!is.na(preds$GPP))/n

  return(list(preds,
              sump,
              cum))
}

# Fit Metabolism ####
filestring <- "metabolism/processed/"
dat <- read_csv(paste0(filestring, "CBP.csv")) %>%
  mutate(DateTime_EST = with_tz(DateTime_UTC, tz = "EST"))
met <- fit_ts_metab_hall(dat, kq_hall, "CBP")

plot(met$date, met$gpp_gO2m2d, type = "l", col = "forestgreen",
     ylim = range(c(met$gpp_gO2m2d, -met$er_gO2m2d), na.rm = T))
lines(met$date, -met$er_gO2m2d, col = "sienna")

out <- filter_hall_met(met, flow_dates)

preds = out[[1]]
met_sum = out[[2]]
cum = out[[3]]

met$ER <- -met$ER
plot_hall_metab(met, error = F)

plot_KvER(preds)

plot_kde_hall_metab(met)
plot_k(preds, xlim = 1)
                