# plot Modis data over time
# Modis summarized for watershed area from GEE
# 2020-12-05
# 1. read data 
dat_list <- list.files("hall_50yl/code/data/gee_files/")


# Full Watershed Data 
# Terrestrial Ecosystem ####

dat <- read_csv("hall_50yl/code/data/gee_files/gpp.csv") %>%
  mutate(year = as.numeric(substr(as.character(date), 1, 4)),
         doy = as.numeric(substr(as.character(date), 5, 8)),
         date = as.Date(paste(year, doy, sep = '-'), format = '%Y-%j')) %>%
  filter(!is.na(date))
GPP <- dat %>%
  group_by(year) %>%
  summarize(GPP_cum_kgC = sum(GPP_mean),
            GPP_std = sum(GPP_std)) %>%
  filter(year != 2020)


dat <- read_csv("hall_50yl/code/data/gee_files/npp.csv") %>%
  select(year = date, NPP_mean = annualNPP_mean, NPP_std = annualNPP_std) %>%
  filter(year != 2020)

NPP <- full_join(GPP, dat)

# LAI 

dat <- read_csv("hall_50yl/code/data/gee_files/lai.csv") %>%
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


riparian <- left_join(NPP, lai, by = "year")
png("figures/Full_watershed_terrestrial_change.png", width = 6, height = 3, 
    units = 'in', type = 'cairo', res = 300)  
riparian <- left_join(NPP, lai, by = "year")
  par(mfrow = c(3,1), mar = c(0,0,1,0), oma = c(4, 4, 4, 0))
  plot(riparian$year, riparian$NPP_mean, xlim = c(1987, 2019),
       ylim = c(min(riparian$NPP_mean - riparian$NPP_std), 
                max(riparian$NPP_mean + riparian$NPP_std)),
       axes = FALSE, type = 'l', lwd = 2)
  polygon(c(NPP$year, rev(NPP$year)), 
          c(NPP$NPP_mean + NPP$NPP_std, 
            rev(NPP$NPP_mean - NPP$NPP_std)),
          col = alpha("forestgreen", .5), border = NA)
  mtext("NHC Watershed Terrestrial Ecosystem", 3, 1, cex = .8)
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

dat <- read_csv("hall_50yl/code/data/nlcd_1992-2016_summary.csv") %>%
  filter(!is.na(category)) %>%
  mutate(category = case_when(id %in% c(22, 23, 24) ~ 'developed',
                           id %in% c(42, 43, 41, 90) ~ 'forested',
                           id %in% c(81, 82) ~ 'agriculture',
                           id %in% c(21, 52, 71) ~ 'grass_shrub',
                           # id %in% c(11) ~ 'water',
                           TRUE ~ 'other')) %>%
  select(-id) %>%
  group_by(category) %>%
  summarize_all(sum, na.rm = T) %>%
  ungroup()
  
ncells <- sum(dat$CellTally1992, na.rm = T)
nlcd <- dat %>%
  filter(category != 'other') %>%
  mutate(across(starts_with("Cell"), ~ . / ncells)) %>%
  pivot_longer(cols = starts_with('Cell'), names_to = 'year', 
               values_to = 'percent') %>%
  mutate(year = as.numeric(substr(year, 10, 13)),
         percent = ifelse(is.na(percent), 0, percent),
         category = factor(category, levels = c('agriculture', 'developed',
                                                'grass_shrub', 'forested')))
nn <- nlcd %>%
  pivot_wider(names_from = category, values_from = percent) %>%
  mutate(grass_shrub = grass_shrub+ forested,
         developed = developed + grass_shrub,
         agriculture = agriculture + developed)

png("figures/nlcd_landcover_change_NHC_watershed.png", width = 5.4, height = 2, 
    units = 'in', type = 'cairo', res = 300)  
  nlcd %>%
  ggplot( aes(year, percent*100, fill = category)) +
    geom_area() + 
    scale_fill_manual(values = c("gray", "black", alpha("forestgreen", .4), 
                                 alpha("forestgreen", .7))) +
    theme_bw() +
    scale_x_continuous(breaks = c(1992, 2000, 2008, 2016)) +
    ylab("Percent Cover")
dev.off()


png("figures/nlcd_landcover_change_NHC_watershed.png", width = 10, height = 9, 
    units = 'in', type = 'cairo', res = 300)  
par(mfrow = c(2,1), mar = c(1,3,1,0), oma = c(4, 1, 4, 0))
  plot(1, type = 'n', ylim = c(0,100), xlim = c(1987, 2019), axes = F)
  polygon(c(nn$year, rev(nn$year)), c(nn$forested, rep(0, nrow(nn)))*100,
          col = alpha('forestgreen', .7), border = NA)
  polygon(c(nn$year, rev(nn$year)), c(nn$grass_shrub, rev(nn$forested))*100,
          col = alpha('forestgreen', .4), border = NA)
  polygon(c(nn$year, rev(nn$year)), c(nn$grass_shrub, rev(nn$developed))*100,
          col = 'black', border = NA)
  polygon(c(nn$year, rev(nn$year)), c(nn$agriculture, rev(nn$developed))*100,
          col = 'grey', border = NA)
  abline(v = c(1992, 2001, 2004, 2006, 2008, 2011, 2013, 2016), 
         col = 'grey30', lty = 2)
  mtext('Percent Landcover', 2, -6.4)
  legend('left', inset = -0.05, xpd = T, 
         legend = c('agriculture', 'developed','shrub/grass','forested',
                    'sample year'),
         fill = c('grey', 'black', alpha('forestgreen', .4), 
                  alpha('forestgreen', .7), NA), 
         col = c(rep(NA, 4), 'grey30'), lty = c(rep(0, 4), 2),
         bty = 'n', border = NA, x.intersp = .3, seg.len = 1.3)
  axis(2, pos = 1991.7, cex.axis = .7,at = seq(0,100, by = 25))
  mtext('Full Watershed Landcover', 3, cex = 1.1, line = .8)
  plot(riparian$year, riparian$NPP_mean, xlim = c(1987, 2019),
       ylim = c(-5000, #min(riparian$NPP_mean - riparian$NPP_std), 
                max(riparian$NPP_mean + riparian$NPP_std)),
       axes = FALSE, type = 'l', lwd = 2)
  mtext('Riparian Zone Productivity', 3, cex = 1.1, line = -.5)
  polygon(c(NPP$year, rev(NPP$year)), 
          c(NPP$NPP_mean + NPP$NPP_std, 
            rev(NPP$NPP_mean - NPP$NPP_std)),
          col = alpha("forestgreen", .5), border = NA)
  mtext("NPP (kgC/m2/y)", 2, 2.4, adj = .88, cex = .9)
  axis(2, cex.axis = .7, at = seq(5000, 10000, by = 2500))
  par(new = T)
  plot(riparian$year, riparian$leaf_on_days, xlim = c(1987, 2019),
       type = "l", lwd = 2, axes = F, ylim = c(150, 227))
  axis(2, pos = 1998, cex.axis = .7,at = seq(170,190, by = 10))
  mtext("Leaf on days", 2, -14.5,adj = .36, cex = .9)
  par(new = T)
  plot(riparian$year, riparian$lai_max, xlim = c(1987, 2019),
       type = "l", lwd = 2, axes = F, ylim = c(47, 79))
  axis(2, pos = 1998, cex.axis = .7, at = seq(48,52, by = 2))
  mtext("Max LAI", 2, -14.5, cex = .9, adj = .03)
  axis(1, cex.axis = .9,at = seq(1986, 2019, by = 3))
  dev.off()
  # mtext("No data", col = 'grey40', side = 3, line = 0, adj = .2, cex = .8)
# climate ####

precip =
  
  
  pp_max10d, er_max10d, gpp_cum, er_cum, nep) %>%
  mutate(across(all_of(c("gpp_median", "er_median",
                         "gpp_max", "er_max",
                         "gpp_cum", "er_cum", "nep")), ~.*14/32))
  