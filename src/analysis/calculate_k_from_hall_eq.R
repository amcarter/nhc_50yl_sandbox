# Calculate k values for NHC based on the equation in Hall 1970

# Appendix
# Equation 4
# k2 = 5.026 * V ^ 0.969 * R(ft) ^ -1.673  at 20 C
# Hall approximates R, hydraulic radius as Depth
# k2(T) = K2(T=20) * 1.0241^(T-20)

# k = 2.3 * k2 * DO.sat / 24hr
library(tidyverse)
library(zoo)
library(streamMetabolizer)
library(lubridate)

# get Hall K values ####
setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/code")

##Convert KO2 to K600
K600fromO2<-function(temp, KO2) {
  ((600/(1800.6-(120.1*temp)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5)*KO2
}

hall_k <- read_csv("../data/hall/hall_tableA2_k_morphology_extended.csv")# %>%
# hall_k <- read_csv("../data/hall/hall_k_compiled.csv") %>%
#   filter(method == "empirical") %>%
#   mutate(K600 = K600fromO2(20, KO2.perday))

ar_k <- read_csv("data/estimated_k_values.csv")
  
# d = depth (m), v = velocity (m/s), DO.sat = DO at saturation (mg/L), t = temp C
setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")

# get measured k values and pair with discharge ####
# nhc <- read_csv("metabolism/processed/NHC.csv") %>%
#   mutate(DateTime_EST = force_tz(DateTime_EST, tz = "EST"),
#          date = as.Date(DateTime_EST)) %>%
#   filter(date == as.Date("2017-06-27")) %>%
#   select(DO.obs, temp.water, level_m, depth, discharge)
# unhc <- read_csv("metabolism/processed/NHC.csv") %>%
#   mutate(DateTime_EST = force_tz(DateTime_EST, tz = "EST"),
#          date = as.Date(DateTime_EST)) %>%
#   filter(date == as.Date("2017-07-11")) %>%
#   select(DO.obs, temp.water, level_m, depth, discharge)
# 
# summary(nhc)
# nhc_k <- ar_k %>%
#   filter(site =="NHC") %>%
#   mutate(depth = .23,
#          discharge = .14, 
#          watertemp = 22,
#          K600 = K600fromO2(watertemp, k_md/depth))
# 
# summary(unhc)
# unhc_k <- ar_k %>%
#   filter(site =="UNHC") %>%
#   mutate(depth = .14,
#          discharge = .02, 
#          watertemp = 25,
#          K600 = K600fromO2(watertemp, k_md/depth))
# 
# ar_k <- bind_rows(nhc_k, unhc_k)
# write_csv(ar_k, "gas_data/measured_k_Ar_releases.csv")


# calc k from Hall 1970 equation ####

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")
# This section is where I figured out how to do this and double checked it looked right
cbp <- read_csv("metabolism/processed/CBP_lvl.csv") %>%
  group_by(date = as.Date(force_tz(DateTime_EST, tz = "EST"))) %>%
  select(date, discharge, depth, level_m, DO.obs, DO.sat, temp_C = temp.water) %>%
  summarize_all(mean, na.rm = T) %>%
  mutate(v_ms = calc_velocity(discharge),
         v2_ms = discharge/depth/18,
         v_fs = v2_ms*3.28084,
         depth_f = depth*3.28084)

cbp_k <- cbp %>%
  mutate(k2_20 = 5.026 * (v_fs)^0.969 * ((depth_f))^(-1.673),
         k2 = k2_20 * 1.0241 ^(temp_C - 20),
         k = (2.3 * (DO.sat) * k2 )/ 24,
         K600 = K600fromO2(temp_C, k2_20),
         K600 = ifelse(K600 > 20, NA, K600))

plot(cbp_k$depth, cbp_k$v2_ms, log = "y")
points(hall_k$depth.m, hall_k$v_ms, col = 2, pch = 19)

plot(cbp_k$depth, cbp_k$k2_20, log = "y")
points(hall_k$depth.m, hall_k$k2_d, col = 2, pch = 19)

plot(cbp_k$depth, cbp_k$k, ylim = c(0,1))
points(hall_k$depth.m, hall_k$k_gm3hr, col = 2, pch = 19)

plot(cbp_k$depth, cbp_k$K600, ylim = c(0,2))
points(hall_k$depth.m, K600fromO2(28,hall_k$k2_d), col = 2, pch = 19)

# Calculate K600 for sites ####

# this width needs to be replaced with actual data from longitudinal surveys
widths <- data.frame(site = c("NHC", "PM", "CBP", "WB", "WBP", "UNHC"),
                    width_m = c(15, 17, 18, 15, 15, 13))

# This loop needs to be stepped through for each individual site. 
kk <- data.frame()
for(site in widths$site){  
    if(site %in% c("NHC", "UNHC")){
      siten = site
    } else {
      siten = paste0(site, "_lvl")
    }
    width <- widths$width_m[widths$site == site]
    
    dat <- read_csv(paste0("metabolism/processed/", siten, ".csv")) %>%
      group_by(date = as.Date(force_tz(DateTime_EST, tz = "EST"))) %>%
      select(date, discharge, depth, level_m,
             DO.obs, DO.sat, temp_C = temp.water) %>%
      summarize_all(mean, na.rm = T) %>%
      mutate(v_ms = calc_velocity(discharge),
             v2_ms = discharge/level_m/width,
             v_fs = v2_ms*3.28084,
             depth_f = depth*3.28084, 
             # calc k2 from Churchill 1962 equation in the Hall 1970 paper
             k2_20 = 5.026 * (v_fs)^0.969 * ((depth_f))^(-1.673),
             # temperature correct k2
             k2 = k2_20 * 1.0241 ^(temp_C - 20),
             # another equation from Hall 1970, appendix 
             k_gm3hr = (2.3 * (DO.sat) * k2 )/ 24,
             # convert to K600 for SM
             K600 = K600fromO2(temp_C, k2_20),
             K600 = ifelse(K600 > 20, NA, K600)) %>%
     # select(date, discharge, depth, temp_C, v_ms, K600, k2_20, k_gm3hr) %>%
      mutate(site = site)
    dat$depth <- dat$level_m
    # plot(dat$depth, dat$v_ms, xlim = c(0,1), ylim = c(0,.5))
    # points(hall_k$depth.m, hall_k$v_ms, col = 2, pch = 19)
    # plot(dat$depth, dat$k2_20, xlim = c(0,1), ylim = c(0,8))
    # points(hall_k$depth.m, hall_k$k2_d, col = 2, pch = 19)
    # plot(dat$depth, dat$K600, pch = 20, main = site, col = "grey50")
    # points(hall_k$depth.m, K600fromO2(20, hall_k$k2_d), col = 2, pch = 20)
    
    plot(dat$discharge, dat$K600, log = "xy", main = site)
    
    kk <- bind_rows(kk, dat)
}

write_csv(kk, "siteData/KQ_hall_prior_from_equation.csv")
# generate K/Q nodes for SM for each site year:
kq <- data.frame()
for(site in unique(kk$site)){
  dat <- kk %>% 
    filter(site == !!site)
  plot(dat$discharge, dat$K600, log = "xy", main = site)
  m <- lm(log(K600) ~ log(discharge), dat)
  mm <- summary(m)$coefficients[,1]
  Qrng <- range(log(dat$discharge), na.rm = T)
  delta = 2
  n = 6
  while(delta > 1){
    n = n + 1
    delta <- (Qrng[2]-Qrng[1])/n
  }  
  Qnodes <- seq(Qrng[1] + delta/2, by = delta, length.out = n)
  lnK600 <- mm[1] + mm[2] * Qnodes
  points(exp(Qnodes), exp(lnK600), col = 2, pch = 19)
  nodes <- data.frame(site = site, 
                      lnQ = Qnodes,
                      lnK600 = lnK600)
  kq <- bind_rows(nodes, kq)
}

write_csv(kq, "siteData/KQ_hall_prior_from_equation.csv")

 # Previous attempt ####
# Hall used the numbers from stream morphology for his calculations, so I will too
hallQ <- log(range(hall_k$discharge_m3s))
# #get the range of all Q's
# q <- read_csv("metabolism/processed/NHC.csv") %>%
#   select(discharge)
# Q <- log(range(q, na.rm = T))

Q <- seq(-9, 5)

comb <- hall_k %>% 
  mutate(logQ = log(discharge_m3s)) %>%
  select(K600, logQ) %>%
  full_join(data.frame(nodes = Q,
                       logQ = Q)) %>%
  arrange(logQ)

comb$K600[which(comb$logQ <= hallQ[1])] <- hall_k$K600[1]
comb$K600[which(comb$logQ >= hallQ[2])] <- 
  hall_k$K600[nrow(hall_k)]

comb_k <- transform(comb, K600 = na.approx(K600, logQ, na.rm = F)) 
 
plot(comb_k$logQ, comb_k$K600)
points(comb_k$nodes, comb_k$K600, col = "red", pch = 19)
write_csv(comb_k, "siteData/KQ_hall_prior.csv")

# 11/19/2020
# didn't finish here the plan was to use the hall data that already was 
# calculated for stream morphology, but that is problematic for a few reasons:
# 1. The range is smaller than that of our data
# 2. the hall rating curve doesn't necesarily apply to this depth
# 
# Another thing to try is calculating the k from my own data, this is 
# a rabbit hole though because the depth is all fucky from using calc_depth
# my plan is to use the level data to modify these depths, but I will need 
# to clean it up a lot. I could also use this to get something like velocity 
# as well.