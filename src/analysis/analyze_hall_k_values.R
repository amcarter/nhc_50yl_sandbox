# Analysis of Hall K values

# library(nleqslv)
library(tidyverse)
# 11/17/2020

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/code")

diurnalk <- read_csv("../data/hall/hall_tableA1_k_diurnalO2_unused.csv") %>%
  mutate(method = "nreg")
morphk <- read_csv("../data/hall/hall_tableA2_k_morphology.csv") %>%
  mutate(method = "empirical", 
         location = "concrete_bridge") 
domek <- read_csv("../data/hall/hall_tableA4_k_dome.csv") %>%
  mutate(method = "dome")

##Convert KO2 to K600
K600fromO2<-function(temp, KO2) {
  ((600/(1800.6-(120.1*temp)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5)*KO2
}

hall_K <- bind_rows(diurnalk, morphk, domek)


# Calculate Hall velocity, temp and DO sat ####
# to make sure you have the equations right

kk <- morphk %>%
  mutate(depth.f = depth.m*3.28084,
         v_fs = (k2_d / (5.026 * depth.f ^ (-1.673)))^(1/.969),
         v_ms = v_fs/3.28084,
         K600 = K600fromO2(20, k2_d)) %>%
  select(-v_fs)

T0 = 2.00856
T1 = 3.224
T2 = 3.99063
T3 = 4.80299
T4 = .978188
T5 = 1.71069
mgL.mlL = 1.42905

Ts <- function(tt) {
  Ts = log((298.15 - tt)/(273.15 + tt))
  return(Ts)
}

k2T <- function(k2, tt) {
  k2T = k2 * 1.0241^(tt - 20)
  return(k2T)
}

#solving this system of equations should give a number for DOsat and temp
# it doesn't though, because hall's equation is based on saturation deficit!
# that means that there are three unknowns and only 2 equations.
# I can't proceed this way, but I can get reasonable numbers for my K600's now.
csat_eqs <- function(tt) {
  csat <- numeric(2)
#  csat[1] = calc_DO_sat(tt, 1013)
  csat[1] <- mgL.mlL * exp(T0 +
                             T1*Ts(tt) +
                             T2*Ts(tt)^2 +
                             T3*Ts(tt)^3 +
                             T4*Ts(tt)^4 +
                             T5*Ts(tt)^5)
  csat[2] = DO + (k * 24) / (2.3 * k2T(tt, k2))
  return(csat)
}


k <- kk$k_gm3hr[1]
k2 <- kk$k2_d[1]
DO = 8
csat0 <- c(20, 9.09)
# this totally didn't work.

nleqslv(csat0, csat_eqs, control = list(allowSingular = TRUE))

csat = csat_eqs(20)

# calculate discharge ####

hall_K <- bind_rows(diurnalk, morphk, domek) %>%
  mutate(KO2.perday = k.mglperhr * 24,
         K600 = K600fromO2(temp = 20, KO2 = KO2.perday))
# hall_K <- bind_rows(diurnalk, morphk, domek) %>%
#   mutate(KO2.perday = k.mglperhr * 24/2.3/calc_DO_sat(20, 1013),
#          K600 = K600fromO2(temp = 20, KO2 = KO2.perday))


hall_rc <- read_csv("../data/hall/hall_figure5_digitized_ratingcurve.csv")
m <- lm(log(hall_rc$discharge_m3s) ~ log(hall_rc$stage_cm))
a <- summary(m)$coefficients[1]
b <- summary(m)$coefficients[2]#Summary of the regression statistics

hall_K$discharge_m3s = exp(a + b * log(hall_K$depth.m*100))
write_csv(hall_K, "../data/hall/hall_k_compiled.csv")

png(filename = "../figures/hall_K600.png",width = 6, height = 5, units = "in",
    res = 300)
ggplot(data = hall_K, aes(log(discharge_m3s), K600)) +
  geom_point(aes(col = method), size = 3) +
  ggtitle("Hall K values") +
  theme_bw()
dev.off()
