## Download USGS gauge station data

# A Carter
# updated: 2020-10-07

# install.packages("dataRetrieval")
library(dataRetrieval) # interact with NWIS and WQP data
library(tidyverse)
library(nhdplusTools)
library(tmap)
library(lubridate)

# library(zoo)
# library(sf)
# library(rgdal)
# library(rgeos)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/code/")

# Get USGS site info from NHD ####
# tmap_mode("view")
# 
# nhc_gage_id <- "USGS-02097314"
# eno_gage_id <- "USGS-02085070"
# 
# #prepping below web quirie below, this just indicates the kind of data we have
# nhc_site <- list(featureSource = "nwissite",featureID = nhc_gage_id) 
# eno_site <- list(featureSource = "nwissite",featureID = eno_gage_id) 
# #returns all upstream flowlines
# nhc_line <- navigate_nldi(nhc_site, "UT", "") 
# eno_line <- navigate_nldi(eno_site, "UT", "") 
# #returns all upstream USGS gages
# nhc_sites_sf <- navigate_nldi(nhc_site, "UT", "nwissite") 
# eno_sites_sf <- navigate_nldi(eno_site, "UT", "nwissite") 
# 
# # now map flowlines and gages
# tm_shape(nhc_line) + tm_lines() +
#   tm_shape(nhc_sites_sf) +tm_dots(col="red") +
#   tm_shape(eno_line) + tm_lines() +
#   tm_shape(eno_sites_sf) +tm_dots(col="red") 
# 
# sites <- bind_rows(eno_sites_sf, nhc_sites_sf)
# site_dat <- sites %>%
#   mutate(lat = st_coordinates(.)[,2],
#          lon = st_coordinates(.)[,1],
#          site_no = substr(identifier, 6, nchar(identifier))) %>%
#   as_tibble() %>%
#   select(name, site_no, identifier, lat, lon)
# site_dat$site_id <- c("eno_hillsborough", "eno_colemill", "eno_roxboro", 
#                       "sandy_cornwallis", "thirdfork_mlk",
#                       "thirdfork_woodcroft", "nhc_blands")

site_dat <- read_csv("../data/USGS_gage_details.csv")

# Download discharge, gage height, and temperature where available ####
# This downloads only daily data, which for some reason is not available for 
#   all of the sites where there is instantaneous data? If I want to do this
#   I'll have to get the rest of the data.

# drop cole mill and sandy_mlk for now, they doesn't have daily discharge
# site_dat <- site_dat[c(1,7),]
# nwis_dat <- readNWISdv(siteNumbers = site_dat$site_no,
                      # parameterCd = c("00060", "00065"))
# nwis_dat <- dataRetrieval::addWaterYear(nwis_dat)
# write_csv(nwis_dat, "data/raw_usgs_download.csv")
# nwis_dat <- read_csv("../data/raw_usgs_download.csv")

nwis <- nwis_dat %>% 
  mutate(discharge_m3s = (0.3048)^3 * X_00060_00003, # discharge in cfs
         gage_height_m = 0.3048 * X_00065_00003      # gage height in ft
         ) %>%
  dplyr::select(site_no, Date, waterYear, 
        discharge_m3s, gage_height_m, cd = X_00060_00003_cd) %>%
  as_tibble()


nwis <- left_join(nwis, site_dat[,c(2,6)]) %>%
  filter(!is.na(site_id))
         site_id == 'nhc_blands') 
nwis %>%
  filter(waterYear>2017) %>%
  mutate(year= year(Date),
         doy = as.numeric(format(Date, '%j')),
         year = ifelse(doy<65, year-1, year)) %>%
  ggplot(aes(Date, discharge_m3s, col = site_id )) +
  geom_line()

ggplot(nwis, aes(Date, log(discharge_m3s))) +
  geom_line() +
  facet_wrap(~site_id, dir = 'v')
  



# RBI calculations ####
# library(remotes)
# remotes::install_github("leppott/ContDataQC")
library(ContDataQC)

nwis <- nwis %>%
  filter(site_id == 'eno_hillsborough')%>%
  mutate(year= year(Date),
         doy = as.numeric(format(Date, '%j')),
         year = ifelse(doy<65, year-1, year)) 
Q_stats <- nwis %>%
  group_by(year) %>%
  summarize(n = length(discharge_m3s), 
            RBI = RBIcalc(discharge_m3s),
            peak_Q = max(discharge_m3s, na.rm = T),
            ar_1 = arima(discharge_m3s, order = c(1,0,0))$coef[1],
            median = median(discharge_m3s, na.rm = T),
            mean = mean(discharge_m3s, na.rm = T),
            cumulative = sum(discharge_m3s, na.rm = T)*60*60*24,
            q05 = quantile(discharge_m3s, .05, na.rm = T)) %>%
  ungroup()%>%
  filter(n >= 365,
         year <2020) # Don't keep incomplete years


write_csv(Q_stats, 'hall_50yl/code/data/nhcblands_usgs_stats.csv')
png(width=7, height=6, units='in', type='cairo', res=300,
    filename='../figures/USGS_rbi_trends.png')
  
  ggplot(Q_stats, aes(waterYear, q05)) +
    geom_smooth() +
    geom_point() +
    facet_wrap(~site_id, dir = 'v')
dev.off()

# Baseflow separation ####
library(RHydro)
baseflow <- nwis %>% 
  dplyr::select(-site_no, -gage_height_m) %>%
  pivot_wider(names_from = site_id, 
              values_from = discharge_m3s, names_prefix = "q_") %>%
  arrange(Date)

eno_q <- baseflow %>%
  filter(!is.na(q_eno_hillsborough)) %>%
  dplyr::select(Date, waterYear, q_eno_hillsborough) %>%
  mutate(baseflow = baseflow_sep(q_eno_hillsborough),
         stormflow = q_eno_hillsborough - baseflow,
         bf_frac = baseflow/q_eno_hillsborough)

monthly_q <- eno_q %>%
  mutate(month = factor(format(Date, "%b"),
                        levels = month.abb)) %>%
  filter(waterYear > 1960, waterYear < 2020) %>%
  group_by(waterYear, month) %>%
  summarize(mean_bf = mean(baseflow, na.rm = T),
            med_bf = median(baseflow, na.rm = T),
            mean_sf = mean(stormflow, na.rm = T),
            med_sf = median(stormflow, na.rm = T),
            bf_frac = mean(bf_frac, na.rm = T),
            qu.05 = quantile(q_eno_hillsborough, .05, na.rm = T),
            qu.10 = quantile(q_eno_hillsborough, .1, na.rm = T),
            q.med = median(q_eno_hillsborough, .5, na.rm = T))

nhc_q <- baseflow %>%
  filter(!is.na(q_nhc_blands),
         waterYear <=2020) %>%
  dplyr::select(Date, waterYear, q_nhc_blands) %>%
  mutate(baseflow = baseflow_sep(q_nhc_blands),
         stormflow = q_nhc_blands - baseflow,
         bf_frac = baseflow/q_nhc_blands)
monthly_q <- nhc_q %>%
  mutate(month = factor(format(Date, "%b"), 
                        levels = month.abb)) %>%
  group_by(waterYear, month) %>%
  summarize(mean_bf = mean(baseflow, na.rm = T),
            med_bf = median(baseflow, na.rm = T),
            mean_sf = mean(stormflow, na.rm = T),
            med_sf = median(stormflow, na.rm = T),
            bf_frac = mean(bf_frac, na.rm = T),
            qu.05 = quantile(q_nhc_blands, .05, na.rm = T),
            qu.10 = quantile(q_nhc_blands, .1, na.rm = T),
            q.med = median(q_nhc_blands, .5, na.rm = T))


png(width=7, height=6, units='in', type='cairo', res=300,
    filename='../figures/nhc_baseflow_frac.png')

ggplot(monthly_q, aes(x = waterYear, y = bf_frac)) +
  #geom_line(aes(y = mean_sf), col = "sienna") +
  geom_point(col = "grey 20") +
  geom_smooth(col = "steelblue", lwd = 1.5)+
  facet_wrap(.~month, scales = "free_y") +
  ggtitle("NHC monthly baseflow fraction")+
  ylab("fraction of water in baseflow")
dev.off()

png(width=7, height=6, units='in', type='cairo', res=300,
    filename='../figures/NHC_median_flow.png')
ggplot(monthly_q, aes(x = waterYear)) +
  geom_point(aes(y = q.med), size = 1) +
  geom_smooth(aes(y = q.med), method = lm, col = "black", lwd = 1.5) +
  geom_point(aes(y = qu.05), col = "steelblue", size = 1) +
  geom_smooth(aes(y = qu.05), method = lm, col = "steelblue", lwd = 1.5)+
  facet_wrap(.~month, scales = "free_y") +
  ggtitle("NHC monthly discharge")+
  ylab("discharge (m3/s), median and 5th percentile")
  


med_bf <- nhc_q %>%
  group_by(waterYear) %>%
  summarize(nhc_mbf = median(baseflow, na.rm=T)) %>%
  ungroup()

#plot(med_bf$waterYear, med_bf$nhc_mbf, type = "l")

# Multipanel graph ####
library(ggpubr)

precip <- read_csv("data/nldas.csv") %>%
  dplyr::select(datetime = DateTime, value = '1', variable) %>% 
  filter(variable == "precipitation_amount") %>%
  mutate(Date = as.Date(datetime)) %>%
  select(Date, precip = value)

preq <- precip %>%
  full_join(nhc_q, by = "Date") %>%
  select(-Date, - bf_frac) %>%
  mutate(p_q = precip/q_nhc_blands) %>%
  group_by(waterYear) %>%
  summarize_all(median, na.rm = T)
cor(preq, use = "na.or.complete")


cum_dat <- precip %>%
  full_join(nhc_q, by = "Date") %>%
  mutate(year = year(Date),
         month = month(Date)) %>%
  group_by(waterYear) %>% #, month) %>%
  summarize(precip_0 = length(which(precip <= 0)),
            precip_mean = mean(precip, na.rm = T),
            precip_max = max(precip, na.rm = T),
            q_median = median(q_nhc_blands, na.rm = T),
            q_05 = quantile(q_nhc_blands, .05, na.rm = T),
            bf_median = median(baseflow, na.rm = T), 
            bf_cum = sum(baseflow),
            sf_mean = mean(stormflow, na.rm = T), 
            sf_cum = sum(stormflow)) %>%
  ungroup() %>%
 # mutate(date = as.Date(paste0(year, "-", month, "-01"), 
 #                       format = "%Y-%m-%d")) %>%
  rename(year = waterYear) %>%
  arrange(year) %>%
  mutate(precip_0 = ifelse(precip_0 > 365, NA, precip_0))


# precip <- read_csv('data/prism/prism_raw.csv') %>%
#   as_tibble() %>%
#   dplyr::select(datetime = DateTime, precip_mm = '1') %>%
#   group_by(year = as.numeric(substr(datetime, 1, 4))) %>%
#   filter(year < 2013) %>%
#   summarize(precip_mm = sum(precip_mm, na.rm = TRUE)) %>%
#   ungroup()



p1 <- ggplot(cum_dat, aes(year, precip_max)) +
        geom_point() +
        geom_smooth(method = lm, lwd = 1, col = "black") +
        theme_minimal() +
        xlab("") +
        ylab("max precip (mm/day)") +
        ggtitle("Precipitation in NHC watershed")
p2 <- ggplot(cum_dat, aes(year, precip_mean)) +
        geom_point() +
        geom_smooth(method = lm, lwd = 1, col = "black") +
        theme_minimal() +
        xlab("") 
p3 <- ggplot(cum_dat, aes(year, precip_0)) +
        geom_point() +
        geom_smooth(method = lm, lwd = 1, col = "black") +
        theme_minimal() +
        xlab("") +
        ylab("Days with 0 precip") 
  
#        xlim(as.Date("1970-01-01"), as.Date("2020-01-01"))
png("../figures/precipitation_trends.png", width = 6, height = 5, 
    units = "in", res = 300)
ggarrange(p1,p3, ncol = 1)
dev.off()
bf <- ggplot(cum_dat, aes(year, bf_mean)) +
        geom_point() +
        ylim(0, 10) + 
        ylab("baseflow (m3s)") +
        xlab("") +
        geom_smooth(method = lm, lwd = 1, col = "black") +
        theme_minimal() +
        xlim(1970, 2020)
#        xlim(as.Date("1970-01-01"), as.Date("2020-01-01")) 

qq <- ggplot(cum_dat, aes(year, q_median)) +
        geom_point() +
        geom_smooth(method = lm, lwd = 1, col = "black") +
        geom_point(aes(y = bf_median), col = "steelblue") +
        geom_smooth(aes(y = bf_median), method = lm, lwd = 1, col = "steelblue") +
        theme_minimal() +
#        ylim(0,6) +
        ylab("discharge (m3s)") +
        xlim(1970, 2020) +
        geom_text(x = 1981, y = 0.93, label = "median flow", 
                  family = "Arial", hjust = 1, size = 4) +
        geom_text(x = 1981, y = 0.715, label = "median baseflow",
                  family = "Arial", hjust = 1, size = 4, col = "steelblue") 
#        xlim(as.Date("1970-01-01"), as.Date("2020-01-01")) 
        
jj <- monthly_q %>%
  filter(month == 7)
jul_bf <- ggplot(jj, aes(waterYear, bf_frac)) +
        geom_point() +
        geom_smooth(method = lm, col = "steelblue", lwd = 1.5) +
        theme_minimal()+
        labs(x = "", y = "July baseflow fraction")
rbi <- ggplot(Q_stats[Q_stats$site_id=="nhc_blands",], aes(waterYear, RBI))+
        geom_point() +
        geom_smooth(method = lm, lwd = 1, col = "black") +
        theme_minimal() +
        xlab("") +
        xlim(1970, 2020)

air_t <- read_csv("../data/noaa_air_temp.csv") %>%
  rename(Date=date)
air_t <- left_join(baseflow_nhc, air_t)
art <- ggplot(air_t, aes(Date, air_trend))+
#        geom_point(col = "grey60") + 
        geom_line() +
        theme_minimal() +
        ylab("air temperature trend C")+
        xlab("")+
        xlim(as.Date("1970-01-01"), as.Date("2020-01-01"))  +
        ggtitle("Trends for NHC at Blands")
png(width=7, height=7, units='in', type='cairo', res=300,
    filename='../figures/nhc_blands_trends_3.png')

  ggarrange(art, qq, rbi, ncol = 1, nrow = 3)

dev.off()

#### contour plots T vs Q

library(ks)
air_t$logQ <- log10(air_t$q_nhc_blands)
wateryears <- c(1988, 1998, 2008, 2018)
cols <- c("grey75", "grey50","grey25","black")

plot(1,1, type = "n", xlab = "air temperature", ylab = "log discharge",
     xlim = range(air_t$air_temp, na.rm=T),
     ylim = range(air_t$logQ), 
     main = "T vs Q kernel density for NHC at Blands") 
for(i in 1:length(wateryears)){
  year <- wateryears[i]
  dat <- air_t[air_t$waterYear==year,]
  kernel <- kde(na.omit(dat[,c("air_temp","logQ")]))
  par(new = T)
  plot(kernel, xlab = "", ylab = "", axes = F, 
      col=cols[i])
}
  