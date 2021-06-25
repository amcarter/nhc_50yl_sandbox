# Calculate k values for NHC based on the equation in Hall 1970

# Appendix
# Equation 4
# k2 = 5.026 * V ^ 0.969 * R(ft) ^ -1.673  at 20 C (day-1, per volume)
# Hall approximates R, hydraulic radius as Depth
# k2(T) = k2(T=20) * 1.0241^(T-20)

# k (g/m3/d) = 2.3 * k2 * DO.sat / 24hr
# K (g/m2/d) = k * d
library(streamMetabolizer)

# get Hall K values ####
##Convert KO2 to K600
calc_k210_from_d_v <- function(d_f, v_fs){
  k2 = 5.026 * v_fs ^ 0.969 * d_f ^ -1.673
  
  return(k2)
}

calc_k_from_k210 <- function(k2, temp, DOsat = NULL, airP = 101.3){
  if(is.null(DOsat)){ DOsat <- calc_DO_sat(temp, airP * 10)}
  k = k2 * 1.024^(temp - 20) * DOsat * 2.3 / 24
  
  return(k)
}

calc_K2_from_K00 <- function(K600, temp){
  sa = 1568
  sb = -86.04
  sc = 2.142
  sd = -0.0216
  se = -0.5
  K2 = K600 * ((sa + sb * temp + sc * temp^2 + sd * temp^3)/600) ^ se
  
  return(K2)
}

calc_K600_from_k210 <- function(k210, temp){
  sa = 1568
  sb = -86.04
  sc = 2.142
  sd = -0.0216
  se = -0.5
  K2 = k210 * 2.3
  K600 = K2 / ((sa + sb * temp + sc * temp^2 + sd * temp^3)/600) ^ se
  
  return(K600)
}

# load hall data
hallk <- read_csv("hall_50yl/code/data/hall/hall_tableA2_k_morphology.csv")

hallk <- hallk %>%
  mutate(k2_day_vol = K2_perday_perarea / depth_m,
         K600 = calc_K600_from_k210(k2_day_vol, 20),
         v_ms = ((k2_day_vol / (5.026 * (depth_m * 3.28) ^ -1.673)) ^ 
                   (1/0.969))/3.28)

write_csv(hallk,"hall_50yl/code/data/hall/hall_tableA2_k_morphology_extended.csv")
# get measured K ####
# ar_k <- read_csv("data/estimated_k_values.csv")
  
# d = depth (m), v = velocity (m/s), DO.sat = DO at saturation (mg/L), t = temp C


# setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")
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
# load processed site data files
# Calculate K600 for sites ####
# 
# # this width needs to be replaced with actual data from longitudinal surveys
# widths <- read_csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/siteData/NHCsite_metadata.csv") %>%
#   select(site = sitecode, width_m = width_mar_m) %>%
#   slice(1:7)

# ZQ <- read_csv("rating_curves/depth_discharge_relationship_LM1953.csv") %>%
#   rename(site = sitename) %>%
#   left_join(widths) %>%
#   mutate(cv_ms = 1/c_m/width_m)
# ZQ$cv_ms[ZQ$site %in% c("NHC", "UNHC")] <- .194
# This loop needs to be stepped through for each individual site. 
kk <- data.frame()
par(mfrow = c(2,1))
for(site in sites$sitecode){  
  dat <- read_csv(paste0("NHC_2019_metabolism/data/metabolism/processed/", 
                         site, ".csv"), guess_max = 10000)
  
  dat <- dat %>%
    group_by(date = as.Date(with_tz(DateTime_UTC, tz = "EST"))) %>%
    select(date, discharge, depth, avg_velocity, DO.obs, 
           DO.sat, temp_C = temp.water) %>%
    summarize_all(mean, na.rm = T) %>%
    mutate(k210_vol = calc_k210_from_d_v(depth * 3.28, avg_velocity * 3.28),
           k_gm3hr = calc_k_from_k210(k210_vol, temp_C, DO.sat),
           K600 = calc_K600_from_k210(k210_vol, temp_C))

    # plot(dat$depth, dat$avg_velocity, pch = 20, main = site, 
    #      col = "grey50", ylim = c(0,1))
    # points(hallk$depth_m, hallk$v_ms, col = 2, pch = 20)
    # 
    # plot(dat$depth, dat$K600, pch = 20, main = site, col = "grey50",
    #      ylim = c(0,20))
    # points(hallk$depth_m, hallk$K600, col = 2, pch = 20)
  dat$site = site
    
    kk <- bind_rows(kk, dat)
}
# as.tibble(kk) %>%
#   ggplot(aes(depth, K600, color= site)) +
#   geom_point()

write_csv(kk, "NHC_2019_metabolism/data/siteData/KQ_hall_prior_from_equation_daily.csv")
# generate K/Q nodes for SM for each site year:
# kq <- data.frame()
# par(mfrow = c(1,1))
# for(site in unique(kk$site)){
#   dat <- kk %>% 
#     filter(site == !!site)
#   plot(dat$discharge, dat$K600, log = "xy", main = site)
#   m <- lm(log(K600) ~ log(discharge), dat)
#   mm <- summary(m)$coefficients[,1]
#   Qrng <- range(log(dat$discharge), na.rm = T)
#   delta = 2
#   n = 6
#   while(delta > 1){
#     n = n + 1
#     delta <- (Qrng[2]-Qrng[1])/n
#   }  
#   Qnodes <- seq(Qrng[1] + delta/2, by = delta, length.out = n)
#   lnK600 <- mm[1] + mm[2] * Qnodes
#   points(exp(Qnodes), exp(lnK600), col = 2, pch = 19)
#   nodes <- data.frame(site = site, 
#                       lnQ = Qnodes,
#                       lnK600 = lnK600)
#   kq <- bind_rows(nodes, kq)
# }
# 
# write_csv(kq, "siteData/KQ_hall_prior_from_equation.csv")
# 
#  # Previous attempt ####
# # Hall used the numbers from stream morphology for his calculations, so I will too
# hallQ <- log(range(hall_k$discharge_m3s))
# # #get the range of all Q's
# # q <- read_csv("metabolism/processed/NHC.csv") %>%
# #   select(discharge)
# # Q <- log(range(q, na.rm = T))
# 
# Q <- seq(-9, 5)
# 
# comb <- hall_k %>% 
#   mutate(logQ = log(discharge_m3s)) %>%
#   select(K600, logQ) %>%
#   full_join(data.frame(nodes = Q,
#                        logQ = Q)) %>%
#   arrange(logQ)
# 
# comb$K600[which(comb$logQ <= hallQ[1])] <- hall_k$K600[1]
# comb$K600[which(comb$logQ >= hallQ[2])] <- 
#   hall_k$K600[nrow(hall_k)]
# 
# comb_k <- transform(comb, K600 = na.approx(K600, logQ, na.rm = F)) 
#  
# plot(comb_k$logQ, comb_k$K600)
# points(comb_k$nodes, comb_k$K600, col = "red", pch = 19)
# write_csv(comb_k, "siteData/KQ_hall_prior.csv")
# 
# # 11/19/2020
# # didn't finish here the plan was to use the hall data that already was 
# # calculated for stream morphology, but that is problematic for a few reasons:
# # 1. The range is smaller than that of our data
# # 2. the hall rating curve doesn't necesarily apply to this depth
# # 
# # Another thing to try is calculating the k from my own data, this is 
# # a rabbit hole though because the depth is all fucky from using calc_depth
# # my plan is to use the level data to modify these depths, but I will need 
# # to clean it up a lot. I could also use this to get something like velocity 
# # as well.