###################################################################
# Download USGS gauge station data

# A Carter
# updated: 2020-10-07

# install.packages("dataRetrieval")
library(dataRetrieval) # interact with NWIS and WQP data
library(tidyverse)
library(nhdplusTools)
library(tmap)

# library(zoo)
# library(sf)
# library(rgdal)
# library(rgeos)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl")

#### Get USGS site info
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

site_dat <- read_csv("data/USGS_gage_details.csv")

#### Download discharge, gage height, and temperature where available
# This downloads only daily data, which for some reason is not available for 
#   all of the sites where there is instantaneous data? If I want to do this
#   I'll have to get the rest of the data.

# drop cole mill and sandy_mlk for now, they doesn't have daily discharge
site_dat <- site_dat[c(-2,-5),]
# nwis_dat <- readNWISdv(siteNumbers = site_dat$site_no,
#                       parameterCd = c("00060", "00065"))
# nwis_dat <- dataRetrieval::addWaterYear(nwis_dat)
# write_csv(nwis_dat, "data/raw_usgs_download.csv")
nwis_dat <- read_csv("data/raw_usgs_download.csv")

nwis <- nwis_dat %>% 
  mutate(discharge_m3s = (0.3048)^3 * X_00060_00003, # discharge in cfs
         gage_height_m = 0.3048 * X_00065_00003      # gage height in ft
         ) %>%
  select(site_no, Date, waterYear, 
        discharge_m3s, gage_height_m) %>%
  as_tibble()

nwis <- left_join(nwis, site_dat[,c(2,6)])

ggplot(nwis, aes(Date, log(discharge_m3s))) +
  geom_line() +
  facet_wrap(~site_id)
  



#### RBI calculations
# library(remotes)
# remotes::install_github("leppott/ContDataQC")
library(ContDataQC)

Q_stats <- nwis %>%
  group_by(site_id, waterYear) %>%
  summarize(n = length(discharge_m3s), 
            RBI = RBIcalc(discharge_m3s),
            peak_Q = max(discharge_m3s, na.rm = T),
            ar_1 = arima(discharge_m3s, order = c(1,0,0))$coef[1]) %>%
  ungroup() %>%
  filter(n >= 365) # Don't keep incomplete years
png(width=7, height=6, units='in', type='cairo', res=300,
    filename='figures/USGS_rbi_trends.png')
  
  ggplot(Q_stats, aes(waterYear, RBI)) +
    geom_smooth() +
    geom_point() +
    facet_wrap(~site_id)
dev.off()

#### Baseflow separation
library(RHydro)
baseflow <- nwis %>% 
  select(-site_no, -gage_height_m) %>%
  pivot_wider(names_from = site_id, 
              values_from = discharge_m3s, names_prefix = "q_")
baseflow <- baseflow[order(baseflow$Date),]
baseflow_nhc <- baseflow[!is.na(baseflow$q_nhc_blands),c(1,2,7)]

baseflow_nhc$baseflow <- baseflow_sep(baseflow_nhc$q_nhc_blands)

med_bf <- baseflow_nhc %>%
  group_by(waterYear) %>%
  summarize(nhc_mbf = median(baseflow, na.rm=T)) %>%
  ungroup()

#plot(med_bf$waterYear, med_bf$nhc_mbf, type = "l")

library(ggpubr)

bf <- ggplot(baseflow_nhc, aes(Date, baseflow))+
        geom_point(col = "grey60") +
        ylim(0,10) +
        ylab("baseflow  (m3s)")+
        geom_smooth(col = "steelblue", lwd = 1.5) +
        theme_minimal() 
rbi <- ggplot(Q_stats[Q_stats$site_id=="nhc_blands",], aes(waterYear, RBI))+
        geom_point(col = "grey25") +
        geom_smooth(col = "steelblue", lwd=1.5) +
        theme_minimal()+
        xlab("")

air_t <- read_csv("data/noaa_air_temp.csv") %>%
  rename(Date=date)
air_t <- left_join(baseflow_nhc, air_t)
art <- ggplot(air_t, aes(Date, air_trend))+
#        geom_point(col = "grey60") + 
        geom_line() +
        theme_minimal() +
        ylab("air temperature trend C")+
        xlab("")+
        ggtitle("Trends for NHC at Blands")
png(width=7, height=6, units='in', type='cairo', res=300,
    filename='figures/nhc_blands_trends.png')

  ggarrange(art, rbi, bf, ncol = 1, nrow = 3)

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
  