# compare nutrients measured in Hall 1970 to today
library(tidyverse)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/")
chem <- read_csv("water_chemistry/water_chemistry_2019-2020_compiled.csv")
 
now <- chem %>% 
  mutate(across(5:15, as.character)) %>%
  pivot_longer(cols = 5:15, names_to = "dissolved", 
               values_to = "concentration") %>%
  mutate(concentration = case_when(concentration == "<0.03" ~ 0.015,
                                   concentration == "<0.01" ~ 0.005,
                                   TRUE ~ as.numeric(concentration))) %>%
  pivot_wider(names_from = "dissolved", values_from = "concentration") %>%
  select(-cl_mgl, -so4_mgl, -br_mgl, -na_mgl, -ca_mgl) %>%
  mutate(doy = format(date, "%j")) %>%
  filter(site != "mc751") %>%
  arrange(doy)

summary(now)
ggplot(now, aes(date, concentration, color = site)) +
  geom_point() +
  facet_wrap(.~dissloved )
setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/data")

hall_p <- read_csv("hall/hall_table_13_p.csv") 
hall_p$tdp_mgl <- NA_real_
for(i in 1:nrow(hall_p)){
  hall_p$tdp_mgl[i] <- mean(c(hall_p$TP_manual_ppm[i], hall_p$TP_mgl[i]), 
                            na.rm = T)
}
hall_n <- read_csv("hall/hall_table_14_nitrogen.csv") %>%
  mutate(no2_mgl = ifelse(no2_mgl == "nd", 0.01, as.numeric(no2_mgl)),
         no3_mgl = ifelse(no3_mgl == "nd", 0.01, as.numeric(no3_mgl)))

hall <- full_join(hall_p, hall_n, by = c("site", "date")) %>%
  mutate(organic_p_mgl_mod = organic_n_mgl * 30 / 14 / 16, # organic P based on redfield
         no3n_mgl = no3_mgl * 14 / (14 + 16 * 3),
         nh3n_mgl = nh3_mgl * 14 / (14 + 2 * 3),
         doy = format(date, "%j"),
         po4p_mgl_mod = (tdp_mgl - organic_p_mgl_mod) * 30 / (30 + 16 * 4)) %>%
  mutate(po4p_mgl_mod = ifelse(po4p_mgl_mod < 0, 0, po4p_mgl_mod)) %>%
  select(site, date, doy, 
         no3n_mgl, nh3n_mgl, don_mgl = organic_n_mgl, tdn_mgl = total_n_mgl,
         tdp_mgl, po4p_mgl_mod, dop_mgl_mod = organic_p_mgl_mod) %>%
  arrange(doy)

plot(now$doy, now$no3n_mgl, type = "l", lty = 2)
points(hall$doy, hall$no3n_mgl, pch = 19)
plot(now$doy, now$nh4n_mgl, type = "l", lty = 2, ylim = c(0,.17))
points(hall$doy, hall$nh3n_mgl, pch = 19)
plot(now$doy, now$po4p_mgl, type = "l", lty = 2, ylim = c(0,.04))
points(hall$doy, hall$po4p_mgl_mod, pch = 19)
points(hall$doy, hall$tdp_mgl)
plot(now$doy, now$tdn_mgl, type = "l", lty = 2)
lines(now$doy, now$no3n_mgl, type = "l", lty = 1)
points(hall$doy, hall$tdn_mgl, pch = 19)

hall_col = "brown3"
now_col = "black"

png("../figures/nutrient_comparison.png", width = 6, height = 5, 
    units = "in", res = 300)  
  par(mfrow=c(2,2), mar = c(4, 4, .5, 0.1))
  plot(density(hall$tdn_mgl, na.rm = T), type = "l", col = hall_col,
       main = "", xlab = "TDN (mg/L)", 
       xlim = c(-.2, 2), lwd = 2)
  
  legend("topright",
         legend = c("TDN 2019", "TDN 1969"),
         col = c(now_col, hall_col),
         lty = 1, bty = "n")
  par(new = T)
  plot(density(now$tdn_mgl, na.rm = T), type = "l", col = now_col,
       xlim = c(-.2, 2), main = "", xlab = "", ylab = "", 
       xaxt = "n", yaxt = "n", lwd = 2)
  plot(density(hall$tdp_mgl, na.rm = T), type = "l", col = hall_col,
       main = "", xlab = "TDP or PO4 (mg/L)", 
       xlim = c(-.02, .3), lwd = 2) 
  par(new = T)
  plot(density(now$po4p_mgl, na.rm = T), type = "l", col = now_col,
       xlim = c(-.02, .3), main = "", xlab = "", ylab = "", 
       xaxt = "n", yaxt = "n", lwd = 2)
  par(new = T)
  plot(density(hall$po4p_mgl_mod, na.rm = T), type = "l", col = hall_col,
       xlim = c(-.02, .3), main = "", xlab = "", ylab = "", lty = 2, 
       xaxt = "n", yaxt = "n", lwd = 2)
  legend("topright",
         legend = c("PO4-P 2019", "TDP 1969", "PO4-P 1969 est"),
         col = c(now_col, hall_col, hall_col),
         lty = c(1,1,2), bty = "n")
  plot(density(hall$no3n_mgl, na.rm = T), type = "l", col = hall_col,
       xlab = "NO3-N (mg/L)", xlim = c(-.1, .4), lwd = 2, main = "")
  legend("topright",
         legend = c("NO3-N 2019", "NO3-N 1969"),
         col = c(now_col, hall_col),
         lty = 1, bty = "n")
  par(new = T)
  plot(density(now$no3n_mgl, na.rm = T), type = "l", col = now_col,
       xlim = c(-.1, .4), main = "", xlab = "", ylab = "", 
       xaxt = "n", yaxt = "n", lwd = 2)
  plot(density(hall$nh3n_mgl, na.rm = T), type = "l", col = hall_col,
       main = "", xlab = "NH4-N (mg/L)", 
       xlim = c(-.05, .25), lwd = 2) 
  par(new = T)
  plot(density(now$nh4n_mgl, na.rm = T), type = "l", col = now_col,
       xlim = c(-.05, .25), main = "", xlab = "", ylab = "", 
       xaxt = "n", yaxt = "n", lwd = 2)
  legend("topright",
         legend = c("NH4-N 2019", "NH4-N 1969"),
         col = c(now_col, hall_col),
         lty = 1, bty = "n")

dev.off()
