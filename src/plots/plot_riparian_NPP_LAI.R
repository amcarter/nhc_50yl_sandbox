# plot Modis data over time
# Modis summarized for watershed area from GEE
# 2020-12-05
# 1. read data 
dat_list <- list.files("hall_50yl/code/data/gee_files/")


# Riparian Data ####

dat <- read_csv("hall_50yl/code/data/gee_files/gpp_riparian.csv") %>%
  mutate(year = as.numeric(substr(as.character(date), 1, 4)),
         doy = as.numeric(substr(as.character(date), 5, 8)),
         date = as.Date(paste(year, doy, sep = '-'), format = '%Y-%j')) %>%
  filter(!is.na(date))
GPP <- dat %>%
  group_by(year) %>%
  summarize(GPP_cum_kgC = sum(GPP_mean),
            GPP_std = sum(GPP_std)) %>%
  filter(year != 2020)


dat <- read_csv("hall_50yl/code/data/gee_files/npp_riparian.csv") %>%
  select(year = date, NPP_mean = annualNPP_mean, NPP_std = annualNPP_std) %>%
  filter(year != 2020)

# NPP <- full_join(GPP, dat)
# plot(GPP$year, GPP$GPP_cum_kgC)
# png("figures/modis_npp.png", 
#     width=6, height=5, units='in', res=300,)
#   plot(NPP$year, NPP$GPP_cum_kgC/1000, type = "l", lwd = 2,
#        ylim = c(min(NPP$NPP_mean - NPP$NPP_std, na.rm = T)/1000,
#                 max(NPP$GPP_cum_kgC + NPP$GPP_std, na.rm = T)/1000), 
#        ylab = "1000 kgC/m2/year", xlab = "year",
#        main = "Terrestrial Productivity")
#   polygon(c(NPP$year, rev(NPP$year)), 
#           c(NPP$GPP_cum_kgC + NPP$GPP_std, 
#             rev(NPP$GPP_cum_kgC - NPP$GPP_std))/1000,
#           col = alpha("forestgreen", .5), border = NA)
#   lines(NPP$year, NPP$NPP_mean/1000, lwd = 2)
#   polygon(c(NPP$year, rev(NPP$year)), 
#           c(NPP$NPP_mean + NPP$NPP_std, 
#             rev(NPP$NPP_mean - NPP$NPP_std))/1000,
#           col = alpha("grey50", .5), border = NA)
#   legend("topleft",
#          c("GPP", "NPP"),
#          col = c("forestgreen", "grey30"),
#          lty = 1, lwd = 2, ncol = 2, bty = "n")
# dev.off()
# 
# summary(lm(NPP_mean~year, data = NPP))



# 2. LAI ####

dat <- read_csv("hall_50yl/code/data/gee_files/lai_riparian.csv") %>%
  filter(!is.na(date)) %>%
  mutate(date = as.Date(date, format = "%Y_%m_%d"),
         year = year(date)) %>%
  select(date, year,lai = Lai_500m_mean) 
dates <- data.frame(date = seq(dat$date[1], dat$date[nrow(dat)], by = 'day'))
max <- dat %>%
  group_by(year) %>%
  summarize(lai_max = quantile(lai, 0.975, na.rm = T)) %>%
  ungroup() 
dat <- dates %>%
  left_join(dat, by = 'date') %>%
  mutate(year = year(date),
         doy = as.numeric(format(date, '%j')),
         lai = na.approx(lai)) %>%
  left_join(max, by = "year") %>%
  mutate(lai = lai/lai_max) %>%
  as.tibble()

lai <- data.frame()
for(y in unique(dat$year)){
  tmp <- dat %>%
    filter(year == y) %>%
    arrange(date)
  w <- which(tmp$lai >= .5)
  don <- tmp$doy[w[1]]
  doff <- tmp$doy[w[length(w)]]
  tmp <- data.frame(year = y, 
                    leaf_on = don,
                    leaf_off = doff) 
  lai <- bind_rows(lai, tmp)
}

lai  <- left_join(lai, max, by = "year") %>%
  mutate(leaf_on_days = leaf_off - leaf_on)
# plot(lai$year, lai$leaf_on_days, type = 'l', col = "sienna", lwd = 1.5)
# lines(lai$year, lai$leaf_on, col = "forestgreen", lwd = 1.5)
# legend('bottomleft',
#        legend = c("Leaf on", "Leaf off"),
#        lty = 1, lwd = 1.5, col = c("forestgreen", "sienna"),
#        bty = 'n')
# 

riparian <- left_join(NPP, lai, by = "year")
png("figures/riparian_zone_terrestrial_change.png", width = 6, height = 3, 
    units = 'in', type = 'cairo', res = 300)  
  par(mfrow = c(3,1), mar = c(0,0,1,0), oma = c(4, 4, 4, 0))
  plot(riparian$year, riparian$NPP_mean, xlim = c(1987, 2019),
       ylim = c(min(riparian$NPP_mean - riparian$NPP_std), 
                max(riparian$NPP_mean + riparian$NPP_std)),
       axes = FALSE, type = 'l', lwd = 2)
  polygon(c(NPP$year, rev(NPP$year)), 
          c(NPP$NPP_mean + NPP$NPP_std, 
            rev(NPP$NPP_mean - NPP$NPP_std)),
          col = alpha("forestgreen", .5), border = NA)
  mtext("Riparian Zone Terrestrial Ecosystem", 3, 1, cex = .8)
  mtext("NPP (kgC/m2/y)", 2, 2.4, cex = .5)
  axis(2, cex.axis = .7, at = seq(5000, 10000, by = 2500))
  plot(riparian$year, riparian$leaf_on_days, xlim = c(1987, 2019),
       type = "l", lwd = 2, axes = F)
  axis(2, cex.axis = .7,at = seq(170,190, by = 10))
  mtext("Leaf on days", 2, 2.4, cex = .5)
  par(mar = c(1,0,0,0))
  plot(riparian$year, riparian$lai_max, xlim = c(1987, 2019),
       type = "l", lwd = 2, axes = F)
  axis(2, cex.axis = .7, at = seq(48,52, by = 2))
  mtext("Max LAI", 2, 2.4, cex = .5)
  axis(1, cex.axis = .7,at = seq(1986, 2019, by = 3))
  mtext("No data", col = 'grey40', side = 3, line = 0, adj = .2, cex = .8)
dev.off()

# Land Cover ####
