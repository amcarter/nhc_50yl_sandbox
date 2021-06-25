# plot precip data for NHC watershed
nldas <- read_csv("hall_50yl/code/data/nldas.csv") %>%
  dplyr::select(datetime = DateTime, value = '1', variable) %>% 
  filter(!variable %in% c("wind_speed",
                          "max_relative_humidity",
                          "min_relative_humidity",
                          "surface_downwelling_shortwave_flux_in_air"))

ggplot(nldas, aes(x=datetime, y = value)) +
  geom_line() +
  facet_wrap(.~variable, scales = "free_y" )

  
precip = read_csv('hall_50yl/code/data/prism/prism_raw.csv') %>%
  as_tibble() %>%
  dplyr::select(datetime = DateTime, precip_mm = '1') %>% 
  pivot_longer(cols = precip_mm, names_to = "variable", values_to = "value")
  # group_by(year = as.numeric(substr(datetime, 1, 4))) %>%
  # filter(year < 2013) %>%
  # summarize(precip_mm = sum(precip_mm, na.rm = TRUE)) %>%
  # ungroup()
all <- bind_rows(nldas, precip) %>% 
  pivot_wider(names_from = "variable", values_from = "value")

all %>%
  rename(precip_daily_mm = precipitation_amount,
         precip_monthly_mm = precip_mm) %>%
  write_csv("hall_50yl/code/data/compiled_nldas_prisim_precip.csv")
            

plot(all$precip_mm, all$precipitation_amount)

p90 <- quantile(all$precipitation_amount, .9, na.rm = T)
pp <- all %>% 
  select(datetime, precip_mmd = precipitation_amount) %>%
  filter(!is.na(precip_mmd)) %>%
  group_by(year = year(datetime)) %>%
  summarize(zero_days = length(which(precip_mmd == 0)),
            cumulative_precip = sum(precip_mmd), 
            cum90 = sum(precip_mmd[which(precip_mmd >= p90)]),
            max_precip = max(precip_mmd, na.rm = T), 
            percent_extreme = cum90/cumulative_precip) %>%
  select(-cum90, -max_precip) %>%
  pivot_longer(cols= -year, names_to = "variable", values_to = "value")
png("figures/precip.png", width = 7.5, height = 6, units = "in", res = 300)
ggplot(pp, aes(x = year, y = value)) +
  geom_point() +
  facet_wrap(.~variable, scales = "free_y", dir = "v", switch = "y") +
  geom_smooth(method = lm, lwd = 1, col = "black") +
  theme_minimal()
dev.off()

# look at Florence
all %>% filter(datetime >= ymd_hms('2018-09-15 00:00:00'),
               datetime <= ymd_hms('2018-09-17 00:00:00')) %>%
  summarize(precip = sum(precipitation_amount))
  ggplot(aes(datetime, precipitation_amount))+
  geom_line()
# calculate precip during 2019 drought ####
p19 <- all %>%
  filter(as.numeric(substr(datetime, 1, 4)) == 2019, 
         as.numeric(substr(datetime, 6, 7)) %in% 9:10) %>%
  slice(c(6:37)) %>%
  select(datetime, precip_mmd = precipitation_amount) %>%
  mutate(pre_cum = cumsum(precip_mmd))
plot(p19$datetime, p19$pre_cum, type = "l")
p19
