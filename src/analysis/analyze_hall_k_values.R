# Analysis of Hall K values

# 11/17/2020

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/code")

diurnalk <- read_csv("../data/hall/hall_tableA1_k_diurnalO2_unused.csv") %>%
  mutate(method = "nreg")
morphk <- read_csv("../data/hall/hall_tableA2_k_morphology.csv") %>%
  mutate(method = "empirical") 
domek <- read_csv("../data/hall/hall_tableA4_k_dome.csv") %>%
  mutate(method = "dome")

##Convert KO2 to K600
K600fromO2<-function(temp, KO2) {
  ((600/(1800.6-(120.1*temp)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5)*KO2
}


hall_K <- bind_rows(diurnalk, morphk, domek) %>%
  mutate(KO2.perday = k.mglperhr * 24,
         K600 = K600fromO2(temp = 20, KO2 = KO2.perday))


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
