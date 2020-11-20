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


setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/code")

#fix this version so that it is the one that has discharge
# you can find how to do that in analyze_hall_k_values
#hall_k <- read_csv("../data/hall/hall_tableA2_k_morphology_extended.csv") 
hall_k <- read_csv("../data/hall/hall_tableA2_k_morphology.csv") %>% 
filter(method == "empirical")

# d = depth (m), v = velocity (m/s), DO.sat = DO at saturation (mg/L), t = temp C
setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")
cbp <- read_csv("metabolism/processed/CBP.csv") %>%
  select(-SpecCond_uScm, -discharge, -WaterTemp_C)

# this should work for calculating things... but it doesnt for some reason
cbp$v_ms <- calc_velocity(cbp$discharge_m3s)
cbp$k2 <- 5.026 * (cbp$v_ms*3.28084)^0.969 * ((cbp$level_m-0.3) * 3.28084)^(-1.673)
cbp$kO2 <- (2.3 * cbp$DO.sat * cbp$k2 )/(cbp$level_m - 0.46)

#hall_k$v_ms <- calc_velocity(hall_k$discharge_m3s)
# hall_k$calc_v <- 
#   hall_k$k2 <- 5.026 * (hall_k$v_ms*3.28084)^0.969 * (hall_k$depth.m * 3.28084)^-1.673
# hall_k$kO2 <- (2.3 * 9 * hall_k$k2 * 1.0241^(20 - 20))/(hall_k$depth.m)

# Hall used the numbers from stream morphology for his calculations, so I will too
hallQ <- range(hall_k$discharge_m3s)
Q <- range(cbp$discharge_m3s, na.rm = T)
comb <- hall_k %>% 
  select(hall_depth = depth.m, hall_kO2 = KO2.perday, discharge_m3s) %>%
  full_join(cbp) %>%
  mutate(logQ = log(discharge_m3s))

comb$hall_kO2[which(comb$discharge_m3s <= hallQ[1])] <- hall_k$KO2.perday[1]
comb$hall_kO2[which(comb$discharge_m3s >= hallQ[2])] <- 
  hall_k$KO2.perday[nrow(hall_k)]

comb_k <- transform(comb, hall_kO2 = na.approx(hall_kO2, logQ)) %>%
  mutate(K600 = K600fromO2(temp.water, hall_kO2))
plot(comb_k$logQ, comb_k$K600)
ggplot(comb_k, aes(hall_kO2, K600))+
  geom_point(aes(col = temp.water))

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