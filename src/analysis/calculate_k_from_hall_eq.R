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
##Convert KO2 to K600
K600fromO2<-function(temp, KO2) {
  ((600/(1800.6-(120.1*temp)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5)*KO2
}

hall_k <- read_csv("../data/hall/hall_tableA2_k_morphology_extended.csv")# %>%
hall_k <- read_csv("../data/hall/hall_k_compiled.csv") %>%
  filter(method == "empirical") %>%
  mutate(K600 = K600fromO2(20, KO2.perday))

# d = depth (m), v = velocity (m/s), DO.sat = DO at saturation (mg/L), t = temp C
setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")
cbp <- read_csv("metabolism/processed/CBP_lvl.csv") 

# this should work for calculating things... but it doesnt for some reason
# cbp$v_ms <- calc_velocity(cbp$discharge)
# cbp$k2 <- 5.026 * (cbp$v_ms*3.28084)^0.969 * ((cbp$depth) * 3.28084)^(-1.673)
# cbp$kO2 <- (2.3 * cbp$DO.sat * cbp$k2 )/(cbp$depth)
# 
# plot(cbp$depth, cbp$k2)#, ylim = c(0,100))
# points(hall_k$depth.m, hall_k$K2, col = 2)

#hall_k$v_ms <- calc_velocity(hall_k$discharge_m3s)
# hall_k$calc_v <- 
#   hall_k$k2 <- 5.026 * (hall_k$v_ms*3.28084)^0.969 * (hall_k$depth.m * 3.28084)^-1.673
# hall_k$kO2 <- (2.3 * 9 * hall_k$k2 * 1.0241^(20 - 20))/(hall_k$depth.m)

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
                       logQ = Q)) 

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