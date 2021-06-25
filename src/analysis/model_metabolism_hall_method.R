# Calculate Metabolism using Hall method

# library(tidyverse)
# library(lubridate)
# library(zoo)
library(pracma)
# setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")
source("NHC_2019_metabolism/src/streamMetabolizer/inspect_model_fits.r")

kq_hall <- read_csv("NHC_2019_metabolism/data/siteData/KQ_hall_prior_from_equation_daily.csv")
flow_dates <- read_csv("NHC_2019_metabolism/data/rating_curves/flow_dates_filter.csv")
met_68 = read_csv('C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/data/hall/hall_table_15.csv')

# test run with Hall data ####
# hall_test = read_csv('hall_50yl/code/data/hall/Hall_figure8_DO_CBP_1970_02_14.csv')
# hall_dodt = read_csv('hall_50yl/code/data/hall/Hall_figure8_DO_CBP_1970_02_14_dodt.csv')
# d <- as.Date("1970-02-14")
# # dodt <- diff(hall_test$DO.obs)/diff(hall_test$hour)
# # day <- data.frame(hour = c(0, hall_test$hour[1:(nrow(hall_test)-1)] + 
# #                              diff(hall_test$hour)/2, 24),
# #                   dodt = c(dodt[1], dodt, dodt[length(dodt)]))
# day <- hall_dodt
# day$k_gm3hr <- 0.58
# sat <- day %>%
#   select(hour) %>% full_join(hall_test) %>% 
#   arrange(hour) %>%
#   select(-DO.obs)
# sat$DO.sat[1] <- hall_test$DO.sat[1] 
# day <- sat %>%
#   mutate(DO.sat = na.approx(DO.sat, x = hour)) %>%
#   right_join(day, by = "hour")
# day <- day %>%
#   mutate(dodt_diff = (1-DO.sat)* k_gm3hr,
#          dodt_cor = dodt - dodt_diff)
# sun <- data.frame(hour = c(6,20))
# 
# day <- day %>%
#   full_join(sun, by = "hour") %>%
#   arrange(hour) %>%
#   mutate(across(-hour, ~na.approx(., x = hour))) 
# 
# sunrise = 6
# sunset = 20
# 
# ps <- day %>%
#   filter(hour >=18 & hour <=22)
# ps_min <- ps$hour[which.min(ps$dodt_cor)]
# a_gpp <- day %>% 
#   filter(hour >= sunrise & hour <= ps_min) %>%
#   select(hour, dodt_cor) %>%
#   mutate(line = NA_real_)
# a_gpp <- data.frame(hour = seq(a_gpp$hour[1], a_gpp$hour[nrow(a_gpp)], by = .25)) %>%
#   full_join(a_gpp, by = "hour") %>%
#   arrange(hour)
# 
# a_gpp$line[c(1, nrow(a_gpp))] <- a_gpp$dodt_cor[c(1, nrow(a_gpp))]
# a_gpp <- a_gpp %>%
#   mutate(line = na.approx(line, hour),
#          dodt_cor = na.approx(dodt_cor, x = hour),
#          dodt_cor = ifelse(dodt_cor < line, line, dodt_cor))
# 
# # build polygon for er that is all negative values of DO/dt
# a_er <- data.frame(hour = seq(0, 24, by = .25)) %>%
#   full_join(day, by = "hour") %>%
#   arrange(hour) %>%
#   select(hour, dodt_cor)
# a_er$dodt_cor[1] <- a_er$dodt_cor[which(!is.na(a_er$dodt_cor))[1]]  
# a_er <- a_er %>% 
#   mutate(dodt_cor = na.approx(dodt_cor, hour),
#          dodt_cor = ifelse(dodt_cor > 0, 0, dodt_cor)) %>%
#   left_join(a_gpp[,c(1,3)], by = "hour") %>%
#   mutate(dodt_cor = case_when(is.na(line) ~ dodt_cor,
#                               dodt_cor < line ~ dodt_cor,
#                               TRUE ~ line))
# start <- data.frame(hour = 0, dodt_cor = 0)
# end <- data.frame(hour = 24, dodt_cor = 0)
# a_er <- bind_rows(start, a_er, end)
# plot(day$hour, day$dodt_cor, type = "l",lty = 1, main = d, lwd = 2, 
#      ylim = range(c(day$dodt, day$dodt_cor), na.rm = T),
#      ylab = "dDO/dt (gO2/m3/hr)", xlab = "hour")
# lines(day$hour, day$dodt, lwd = 2,lty = 2,)
# abline(h = 0)
# abline(v = c(sunrise, sunset), lwd = 2, 
#        lty = 3, col = "goldenrod")
# # lines(day$hour, day$dodt_diff, lty = 2)
# 
# polygon(a_gpp$hour, a_gpp$dodt_cor, border = "black",
#         density = 20, angle = 45, 
#         col = alpha("forestgreen", 1))
# polygon(a_er$hour, a_er$dodt_cor, border = NA,
#         density = 20, angle = -45,
#         col = alpha("sienna", 1))
# # points(day$hour, day$dodt_pred, pch = 19, col = 5)
# legend("topright", cex = .5,
#        legend = c("dDO/dt uncorrected", "dDO/dt corrected", 
#                   "sunrise/set",
#                   "GPP", "ER"), 
#        fill = c(NA, NA, NA, "forestgreen","sienna"),
#        density = c(0, 0, 0, 45, 45),
#        angle = c(NA, NA, NA, 45, -45),
#        col = c(1, 1, "goldenrod", NA, NA),
#        border = c(NA, NA, NA, 1,1), 
#        lty = c(2, 1, 3, NA, NA),
#        x.intersp = c(2, 2, 2, 0.5, 0.5),
#        bty = "n")
# depth = .4
# m <- data.frame(date = d,
#                 gpp_gO2m3d = abs(polyarea(a_gpp$hour, a_gpp$dodt_cor))*depth,
#                 er_gO2m3d = abs(polyarea(a_er$hour, a_er$dodt_cor))*depth,
#                 k_gm3hr = day$k_gm3hr[1])

# Metabolism functions ####
# calculate hall met given k2 and a day of data spaced at dt min intervals

fit_metab_hall <- function(day, dt = 15, plot = TRUE){
  # calculate the slope (in mgo2/m3/hr) at each point
  d <- as.Date(day$DateTime_EST[1])
  s <- day$site[1]
  day$dodt <- c(day$DO.obs[2] - day$DO.obs[1], diff(day$DO.obs)) * 60/dt # mgO2/m3/hr
  
  day <- day %>%
    mutate(dodt_diff = (1 - DO.obs/DO.sat) * k_gm3hr, # diffusion due to sat def
           dodt_cor = dodt - dodt_diff,             # dodt corrected for diff
           hour = (seq(1:nrow(day))-1)/60 * dt)
  # if missing data, cannot calc met (these ts are already gap filled up to 3 hrs)
  if(sum(is.na(day$dodt_cor)) > 0) { 
    m <- data.frame(date = d)
    print(paste(d, "missing data"))
    return(m)
  }
  
  # build met polygons ####
  # find sunrise and sunset based on light
  sun <- which(day$light > 0)
  sunrise <- sun[1]-1
  sunset <- sun[length(sun)]+1
  # find minimum change in two hours post sunset
  ps_min <- sunset -2 + which.min(day$dodt_cor[(sunset-60/dt):(sunset+2*60/dt)])
  
  # build polygon for gpp by connecting sunrise to post sunset minimum
  a_gpp <- day %>%
    slice((sunrise):(ps_min)) %>%
    select(hour, dodt_cor) %>%
    mutate(line = NA_real_)
  a_gpp <- data.frame(hour = seq(a_gpp$hour[1], a_gpp$hour[nrow(a_gpp)], by = .25)) %>%
    left_join(a_gpp, by = "hour")
  
  if(a_gpp$dodt_cor[1] > 0){
    a_gpp <- bind_rows(data.frame(hour = a_gpp$hour[1], 
                      dodt_cor = 0), a_gpp)
  }
  if(a_gpp$dodt_cor[nrow(a_gpp)] > 0) {
    a_gpp <- bind_rows(a_gpp, data.frame(hour = a_gpp$hour[nrow(a_gpp)],
                                         dodt_cor = 0))
  }
  a_gpp$line[c(1, nrow(a_gpp))] <- a_gpp$dodt_cor[c(1, nrow(a_gpp))]
  a_gpp <- a_gpp %>%
    mutate(line = na.approx(line),
           dodt_cor = na.approx(dodt_cor),
           dodt_cor = ifelse(dodt_cor < line, line, dodt_cor))

  # build polygon for er that is all negative values of DO/dt
  a_er <- data.frame(hour = seq(0, 24, by = .25)) %>%
    left_join(day, by = "hour") %>%
    select(hour, dodt_cor) %>%
    mutate(dodt_cor = na.approx(dodt_cor),
           dodt_cor = ifelse(dodt_cor > 0, 0, dodt_cor)) %>%
    left_join(a_gpp[,c(1,3)], by = "hour") %>%
    mutate(dodt_cor = case_when(is.na(line) ~ dodt_cor,
                                dodt_cor < line ~ dodt_cor,
                                TRUE ~ line))
  start <- data.frame(hour = 0, dodt_cor = 0)
  end <- data.frame(hour = 24, dodt_cor = 0)
  a_er <- bind_rows(start, a_er, end)
  if(plot) {
  # par(mfrow = c(3,1))
  # plot(day$hour, day$DO.obs, type = "b")
  # plot(day$hour, day$DO.obs/day$DO.sat, type = "l")
  # abline(h = 1)
  plot(day$hour, day$dodt_cor, type = "l",lty = 1, 
       main = paste(d, s), lwd = 2, 
       ylim = range(c(day$dodt, day$dodt_cor), na.rm = T),
       ylab = "dDO/dt (gO2/m3/hr)", xlab = "hour")
  lines(day$hour, day$dodt, lwd = 2,lty = 2,)
  abline(h = 0)
  abline(v = day$hour[c(sunrise, sunset)], lwd = 2, 
         lty = 3, col = "goldenrod")
  # lines(day$hour, day$dodt_diff, lty = 2)

  polygon(a_gpp$hour, a_gpp$dodt_cor, border = "black",
          density = 20, angle = 45, 
           col = alpha("forestgreen", 1))
  polygon(a_er$hour, a_er$dodt_cor, border = NA,
          density = 20, angle = -45,
          col = alpha("sienna", 1))
  # points(day$hour, day$dodt_pred, pch = 19, col = 5)
  legend("topright", cex = .5,
         legend = c("dDO/dt uncorrected", "dDO/dt corrected", 
                    "sunrise/set",
                    "GPP", "ER"), 
          fill = c(NA, NA, NA, "forestgreen","sienna"),
         density = c(0, 0, 0, 45, 45),
         angle = c(NA, NA, NA, 45, -45),
         col = c(1, 1, "goldenrod", NA, NA),
         border = c(NA, NA, NA, 1,1), 
         lty = c(2, 1, 3, NA, NA),
         x.intersp = c(2, 2, 2, 0.5, 0.5),
         bty = "n")
  }
  
  K600 = day$K600[1]
  m <- data.frame(date = d,
                  gpp_gO2m3d = abs(polyarea(a_gpp$hour, a_gpp$dodt_cor)),
                  er_gO2m3d = abs(polyarea(a_er$hour, a_er$dodt_cor)),
                  k_gm3hr = day$k_gm3hr[1],
                  K600 = K600)
  
  
  return(m)
}
fit_ts_metab_hall<- function(dat, kq_hall, s, plot = FALSE, dt = 15){
  d1 <- as.Date(dat$DateTime_EST[1]) +1
  d2 <- as.Date(dat$DateTime_EST[nrow(dat)]) -1
  
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
    day$k_gm3hr <- kq$k_gm3hr[1]
    day$K600 <- kq$K600[1]
    
    # m <- fit_metab_hall(day, k2, plot = T, dt = 60)
    m <- fit_metab_hall(day, plot = plot, dt = dt)
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
filter_hall_met <- function(met, flow_dates){
  preds <- met %>%
    select(date, discharge, depth, temp.water,
           GPP = gpp_gO2m2d, ER = er_gO2m2d, K600) %>%
    left_join(flow_dates[,c(1,4)]) %>%
    mutate(across(all_of(c("GPP", "ER", "K600")), 
                  ~ ifelse(good_flow, ., NA)))
  coverage <- data.frame(missing_data = sum(is.na(met$gpp_gO2m2d)),
                         bad_flow = sum(preds$good_flow == FALSE &
                                          !is.na(met$gpp_gO2m2d), na.rm = T),
                         total_days = nrow(preds))
  w <- range(which(!is.na(preds$GPP)), na.rm = T)
  preds <- preds[w[1]:w[2],]
  
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
  sump <- bind_cols(coverage, sump)
  
  return(list(preds,
              sump,
              cum))
}
get_met_year<- function(dat, kq_hall, s, flow_dates, year, dt = 15, plot = FALSE){  
  met <- fit_ts_metab_hall(dat, kq_hall, s, dt = dt, plot = plot)
  # plot(met$date, met$gpp_gO2m2d, type = "l", col = "forestgreen",
  #      ylim = range(c(met$gpp_gO2m2d, -met$er_gO2m2d), na.rm = T))
  # lines(met$date, -met$er_gO2m2d, col = "sienna")
  out <- filter_hall_met(met, flow_dates)
  
  preds = out[[1]] %>%
    mutate(site = s,
           year = year, 
           ER = -ER)
  met_sum = out[[2]]%>%
    mutate(site = s,
           year = year)
  cum = out[[3]]%>%
    mutate(site = s,
           year = year, 
           ER = -ER)
  
  m <- matrix(c(1,2,1,3,1,4), nrow = 2)
  layout(m)
  plot_hall_metab(preds, error = F)
  mtext(s)
  plot_kde_hall_metab(preds, lim = 7)
  plot_KvER(preds)
  plot_k(preds)
     
  return(list(preds, met_sum, cum))
}           

# Fit Metabolism ####
filestring <- "NHC_2019_metabolism/data/metabolism/processed/"

preds_all <- data.frame()
met_summary <- data.frame()
preds_filled_all <- data.frame()

for(s in sites$sitecode){
 
  dat <- read_csv(paste0(filestring, s, ".csv"), guess_max = 10000) %>%
    mutate(DateTime_EST = with_tz(DateTime_UTC, tz = "EST"))
  
  # run on 60 min intervals
  dat <- dat %>%
    filter(substr(DateTime_EST, 15, 19) == "00:00")
  
  # run on smoothed data
  # dat <- dat %>%
  #   select(DateTime_EST, DO.obs, discharge, temp.water, DO.sat, depth, light) %>%
  #   mutate(across(c(-DateTime_EST, -light), ~ 
  #                   stats::filter(., rep(1/3, 3), sides = 2)))
  # seq <- seq(1, nrow(dat), by = 4)
  # dat <- dat[seq,]
  if(s %in% c("NHC", "UNHC")){
    for(y in 2017:2019){
      yy <- ymd_hms(paste0(y, "-03-01 00:00:00"), tz = "EST")
      dat1 <- dat %>%
        filter(DateTime_EST >= yy, 
               DateTime_EST <= yy + 24*60*60*366)
      year = y
      out <- get_met_year(dat1, kq_hall, s, flow_dates, year, dt = 60) 
      preds_all <- bind_rows(preds_all, out[[1]])
      met_summary <- bind_rows(met_summary, out[[2]])
      preds_filled_all <- bind_rows(preds_filled_all, out[[3]])
    }
    next
  }
  if(s == "WBP"){
    dat <- dat %>% filter(DateTime_EST <= ymd_hms("2020-03-20 00:00:00", 
                          tz = "EST"))
  }
  
  year = 2019
  
  out <- get_met_year(dat, kq_hall, s, flow_dates, year, dt = 60, plot = FALSE) 
  preds_all <- bind_rows(preds_all, out[[1]])
  met_summary <- bind_rows(met_summary, out[[2]])
  preds_filled_all <- bind_rows(preds_filled_all, out[[3]])
}

saveRDS(list(preds = preds_all, 
             summary = met_summary, 
             cumulative = preds_filled_all),
        "NHC_2019_metabolism/data/metabolism/hall/hall_met_60min_2021_01.rds")

# # inspect Metabolism ####
# met_15 <- readRDS("metabolism/hall/hall_met_15min.rds")
# met_60 <- readRDS("NHC_2019_metabolism/data/metabolism/hall/hall_met_60min_2021_01.rds")
# met_sm <- readRDS("metabolism/hall/hall_met_15min_smoothed.rds")
# met_ray <- readRDS("NHC_2019_metabolism/data/metabolism/compiled/raymond_met.rds")
# 
# write_csv(bind_rows(met_ray$summary, met_60$summary),
#           "metabolism/compiled/metabolism_summary_tables_v2.csv")
# 
# png("figures/tmp/met_estimate_method_comparison_kernels.png",
#   width = 10, height = 5, units = "in", res = 300)
#   par(mfrow = c(1,2))
#   plot_kde_metab(met_ray$preds, col = "steelblue", lim = 6)
#   par(new = T)
#   plot_kde_metab(met_60$preds, lim = 6)
#   legend("topright",
#          c("Hall method", "Stream Met"),
#          fill = c("grey", "steelblue"),
#          border = NA, bty = 'n', cex = 1.2)
#   plot_kde_hall_metab(met_60$preds, lim = 6)
# dev.off()
# 
# met <- met_15$preds %>%
#   as_tibble() %>%
#   select(date, site, discharge, depth, temp = temp.water,
#          gpp15 = GPP, er15 = ER) %>%
#   left_join(met_60$preds[,c(1,5,6,9)], by = c("date", "site")) %>%
#   rename(gpp60 = GPP, er60 = ER) %>%
#   left_join(met_sm$preds[,c(1,5,6,9)], by = c("date", "site")) %>% 
#   rename(gppsm = GPP, ersm = ER)
# png("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/figures/hall_metab_sample_frequency_comparison.png",
#     width = 5, height = 4, units = "in", res = 300)
#   par(mfrow = c(1,2))
#   plot(-met$er60, -met$er15, pch = 20, cex = .7, 
#        xlab = "60 min ER", ylab = "15 min ER")
#   points(-met$er60, -met$ersm, pch = 20, col = 5, cex = .7)
#   abline(0,1, col = 2)
#   plot(met$gpp60, met$gpp15, pch = 20, cex = .7, 
#        xlab = "60 min GPP", ylab = "15 min GPP")
#   points(met$gpp60, met$gppsm, pch = 20, col = 5, cex = .7)
#   abline(0,1, col = 2)
#   legend("topleft", 
#          c("raw", "smoothed"),
#          pch = 20, col = c(1,5), bty = "n")
#   par(new = T, mfrow = c(1,1))
#   mtext("Metabolism calculated as in Hall 1972", line = 1)
# dev.off()
# 
# plot_kde_hall_metab(met_15$preds, lim = 8)
# par(new = T)
# plot_kde_metab(met_sm$preds, col = "steelblue", lim = 8)
# par(new = T)
# plot_kde_hall_metab(met_60$preds, lim = 8, site = "CBP")
# 
# 
# tmp <- met_68 %>%
#   mutate(pr = GPP_gO2m2d/ER_gO2m2d) %>%
#   arrange(pr)
# 
# dat <- met_60$preds %>%
#   mutate(pr = -GPP/ER) %>%
#   arrange(pr)
# 
# summary(tmp)
# min(which(tmp$pr >= 1))/sum(!is.na(tmp$pr))
# median(tmp$GPP, na.rm = T)/median(tmp$ER, na.rm = T)
# median(tmp$GPP_gO2m2d, na.rm = T)/median(tmp$ER_gO2m2d, na.rm = T)
# median(met_60$preds$GPP, na.rm = T)/median(met_68$GPP_gO2m2d)
# median(met_60$preds$ER, na.rm = T)/median(met_68$ER_gO2m2d)
# 
# w(tmp)
# tmp$index = seq(1:nrow(tmp))
# plot(tmp$GPP, -tmp$ER)
# abline(0,1, col = 2)
# plot(log(tmp$pr), ylim = c(-5,5), pch = 20)
# abline(h=0)
# 
# ggplot(tmp, aes(index,  log(pr), color = site)) +
#   geom_point()
# 
# for(i in c(1:6)){
#   s = sites$sitecode[i]
#   tmp <- dat %>%
#     filter(site == s,
#            !is.na(pr),
#            year == 2019) %>%
#     arrange(pr) 
#   plot(log(tmp$pr), ylim = c(-3,2), type = "l", col = 1, lwd = 2, 
#        axes = F, ylab = "log(P/R)")
#   par(new = T)
# }
# abline(h = 0)
# axis(2)
# legend("topleft", 
#        legend = c("2019", "1969"), lty = 1, lwd = 2, cex = 1.2,
#        col = c(1, "brown3"), bty = "n")
# 
# legend("topleft", 
#        legend = sites$sitecode[c(1:6), drop = T], 
#        lty = 1, lwd = 2, cex = 1.2,col = c(1:6), bty = "n")
