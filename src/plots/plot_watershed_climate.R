# plot precip data for NHC watershed

pp_mon = read_csv('hall_50yl/code/data/prism/prism_raw.csv') %>%
  mutate(month = as.numeric(month(DateTime)),
         year = year(DateTime)) %>%
  dplyr::select(year, month, ppt_mm = '1') %>%
  filter(year < 1981)
  # pivot_longer(cols = precip_mm, names_to = "variable", values_to = "value") 
pp_day <- read_csv('hall_50yl/code/data/gee_files/ppt.csv') %>%
  filter(!is.na(date)) %>%
  mutate(year = as.numeric(substr(date, 1, 4)),
         month = as.numeric(substr(date, 5,6)),
         date = as.Date(paste(year, month, substr(date, 7, 8), sep = '-',
                              format = '%Y-%m-%d'))) %>%
  select(year, month, date, ppt_mm = ppt_mean)
pp_mon <- pp_day %>%
  group_by(year, month) %>%
  summarize(ppt_mm = sum(ppt_mm)) %>%
  bind_rows(pp_mon)
pp_year <- pp_mon %>%
  group_by(year) %>%
  summarize(ppt_cum = sum(ppt_mm)) 

pp_mon %>% filter(month == 9) %>% ggplot(aes(year, ppt_mm)) +geom_line()

p90 <- quantile(pp_day$ppt_mm, .9, na.rm = T)
pp <- pp_day %>% 
  group_by(year) %>%
  summarize(zero_days = length(which(ppt_mm == 0)),
            cum90 = sum(ppt_mm[which(ppt_mm >= p90)]),
            max_precip = max(ppt_mm, na.rm = T)) %>%
  ungroup() %>%
  full_join(pp_year, by = 'year') %>%
  mutate(percent_extreme = cum90/ppt_cum) %>%
  select(-cum90) %>%
  arrange(year)
  # pivot_longer(cols= -year, names_to = "variable", values_to = "value")
write_csv(pp, "hall_50yl/code/data/precip_annual_summary.csv")
write_csv(pp_mon, "hall_50yl/code/data/precip_monthly_summary.csv")


png("figures/watershed_climate.png", width = 10, height = 9, 
    units = 'in', type = 'cairo', res = 300)  


  par(mfrow = c(3,1), mar = c(0,2,0,2), oma = c(4, 2, 4, 0))
  # not changing
  plot(pp$year, pp$ppt_cum/1000, type = 'l', axes = F)
  axis(2, cex.axis = .7)
  mtext("annual precip (m)", 2, 2.4, cex = .9)
  mtext('Precipitation', 3, cex = 1.1, line = .8)
  plot(pp$year, pp$zero_days, type = 'l', axes = F)
  
  plot(pp$year, pp$percent_extreme, type = 'l', axes = F)
  axis(2)
  mtext("annual precip (m)", 2, 2.4, cex = .8)
  mtext("Not Changing", 3)
  value)) +
  geom_point() +
  facet_wrap(.~variable, scales = "free_y", dir = "v", switch = "y") +
  geom_smooth(method = lm, lwd = 1, col = "black") +
  theme_minimal()
dev.off()


# calculate precip during 2019 drought ####
p19 <- all %>%
  filter(as.numeric(substr(datetime, 1, 4)) == 2019, 
         as.numeric(substr(datetime, 6, 7)) %in% 9:10) %>%
  slice(c(6:37)) %>%
  select(datetime, precip_mmd = precipitation_amount) %>%
  mutate(pre_cum = cumsum(precip_mmd))
plot(p19$datetime, p19$pre_cum, type = "l")
p19
