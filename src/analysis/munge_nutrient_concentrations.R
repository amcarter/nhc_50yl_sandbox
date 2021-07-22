# Compare nutrients measured in Hall to nutrients in NHC today

setwd("C:/Users/Alice Carter/git/nhc_50yl/hall_50yl/code/")
library(tidyverse)
library(lubridate)
library(zoo)

# 1. Load Hall Data ####
# table 26 has TP data

t26 <- read_csv("data/hall/hall_table_26.csv") %>%
  dplyr::rename(site = location)
t13 <- read_csv("data/hall/hall_table_13_p.csv")
t14 <- read_csv("data/hall/hall_table_14_nitrogen.csv") %>%
  mutate(no2_mgl = ifelse(no2_mgl == "nd", 0, as.numeric(no2_mgl)), 
         no3_mgl = ifelse(no3_mgl =="nd", 0, as.numeric(no3_mgl)))
hall <- full_join(t26, t13, by = c("date", "site", "TP_mgl")) %>%
  full_join(t14, by = c("date", "site")) %>%
  mutate(no3n_mgl = no3_mgl *14/(14+16*3),
         nh3n_mgl = nh3_mgl *14/17)

dat_mar <- read_csv("C:/Users/Alice Carter/git/nhc_50yl/NHC_2019_metabolism/data/longitudinal_sampling/NHCLongitudinalDO_20190308.csv") %>%
  mutate(depth_cm = as.numeric(depth_cm), 
         NO3.N_mgl = NO3.N_mgl/10)
dat_oct <- read_csv("C:/Users/Alice Carter/git/nhc_50yl/NHC_2019_metabolism/data/longitudinal_sampling/NHCLongitudinalDO_20191009.csv")

now <- bind_rows(dat_mar, dat_oct) %>%
  select(date, distance_m, depth_cm, temp_C, DO_mgL, SpC_uscm, 
         no3n_mgl = NO3.N_mgl, NPOC_mgl, tdn_mgl = TN_mgl, nh4n_mgl = NH4.N_mgl, 
         po4p_mgl = PO4.P_mgl)

summary(hall)
summary(now)


setwd("C:/Users/Alice Carter/git/ghg_patterns_nhc/")

nuts <- read_csv("data/water_chemistry/water_chemistry_2019-2020_compiled.csv") %>%
  select(-sample_name, -time) %>%
  mutate(site = toupper(site)) %>%
  filter(site != "MC751")

SP_wchem <- read_csv('data/water_chemistry/StreampulseWQDec2020.csv') %>%
  filter(site %in% c('NHC', 'UNHC')) %>%
  select(site, date, cl_mgl = Cl, so4_mgl = 'SO4 (mg/L)', br_mgl = Br, 
         no3n_mgl = 'NO3-N', na_mgl = 'Na (mg/L)', k_mgl = 'K (mg/L)', 
         mg_mgl = 'Mg (mg/L)', ca_mgl = 'Ca (mg/L)', doc_mgl = 'DOC (mg/L)', 
         tdn_mgl = 'TDN (mg/L)', nh4n_mgl = 'NH4-N (mg/L)', 
         po4p_mgl = 'PO4-P (mg/L)') %>%
  mutate(date = as.Date(date, format = '%m/%d/%Y'),           
         no3n_mgl = ifelse(no3n_mgl == "<0.001", 0.0005, 
                           as.numeric(no3n_mgl))) %>%
  bind_rows(nuts) %>%
  mutate(br_mgl = ifelse(br_mgl == "<0.03", 0.015, as.numeric(br_mgl)),
         nh4n_mgl = ifelse(nh4n_mgl == "<0.01", 0.005, as.numeric(nh4n_mgl)),
         po4p_mgl = ifelse(po4p_mgl == "<0.01", 0.005, as.numeric(po4p_mgl)),
         date = case_when(site == 'NHC' & date == as.Date('2020-01-29') ~
                            as.Date('2020-01-30'),
                          date == as.Date('2020-03-10') ~ 
                            as.Date('2020-03-11'),
                          TRUE ~ date)) %>%
  filter(date <as.Date('2020-04-01'))



spchem <- read_csv("data/water_chemistry/all_grab_data.csv") %>%
  filter(siteID %in% c("NHC", "UNHC")) %>%
  dplyr::select(-flagID, -flagComment, -methodDetail, -writeInMethod, -regionID, -method) %>% 
  group_by(siteID, dateTimeUTC, variable) %>%
  summarize_all(mean, na.rm = T) %>%
  ungroup() %>%#data.frame()
  filter(!is.na(siteID), !is.na(dateTimeUTC),
         !(variable %in% c('Br'))) %>%
  rename(site = siteID, DateTime_UTC = dateTimeUTC)%>%
  pivot_wider(values_from = value, names_from = variable) %>%
  mutate(date = as.Date(DateTime_UTC)) %>%
  mutate(nh4n_mgl = NH4 * 14 * 1000, 
         no3n_mgl = NO3 * 14 * 1000, 
         po4p_mgl = PO4 * 31 *1000, 
         doc_mgl = TOC, 
         tdn_mgl = TDN ) %>%
  select(site, date, nh4n_mgl, no3n_mgl, po4p_mgl, doc_mgl, tdn_mgl) %>%
  bind_rows(SP_wchem) %>%
  mutate(po4p_mgl = ifelse(po4p_mgl > .3, NA, po4p_mgl))

# pair with discharge and temperature
dat <- read_csv('C:/Users/Alice Carter/git/nhc_50yl/NHC_2019_metabolism/data/rating_curves/interpolatedQ_allsites_modified.csv') %>%
  rename(NHC = NHC.Q, UNHC = UNHC.Q) %>%
  dplyr::select(DateTime_UTC, NHC, UNHC) %>%
  pivot_longer(cols = -DateTime_UTC, names_to = "site", values_to = 'discharge') %>%
  mutate(datetime_est = with_tz(DateTime_UTC, tz = "EST")) %>%
  group_by(site, date = as.Date(datetime_est)) %>%
  summarize(discharge = mean(discharge, na.rm = T))

chem <- spchem %>% 
  left_join(dat, by = c('site', 'date')) %>%
  mutate(doy = as.numeric(format(date, '%j')))
chem %>%
  select( site, date, discharge, nh4n_mgl, no3n_mgl, po4p_mgl, doc_mgl, tdn_mgl)%>%
  pivot_longer(cols = ends_with('mgl'), names_to = 'variable', values_to = 'value') %>%
  ggplot(aes(date, value, col = site)) +
  geom_point()+
  # geom_smooth(method = lm) +
  facet_wrap(.~variable, scales = "free_y")

chem <- dat %>%
  ungroup()%>%
  filter(site == 'NHC') %>% 
  select(-site)%>%
  right_join(SP_wchem, by = c('date')) %>%
  mutate(doy = as.numeric(format(date, '%j')),
         po4p_mgl = ifelse(po4p_mgl > .3, NA, po4p_mgl))
ggplot(chem, aes(x = date, y = doc_mgl, col = site))+
  geom_line()+ geom_point()
now <- bind_rows(now, chem)
summary(now)

# compare nuts to metabolism ####
met <- readRDS("C:/Users/Alice Carter/git/nhc_50yl/NHC_2019_metabolism/data/metabolism/compiled/met_preds_stream_metabolizer.rds")$preds %>%
  filter(era == 'now', site !='PWC', year == 2019) %>%
  select(year, date, site, GPP, ER, discharge, DO_mgl = DO.obs, DO.sat) %>%
  mutate(DO_mgl = DO_mgl/DO.sat)
filled <- chem %>% 
  select(site, date, no3n_mgl, doc_mgl, tdn_mgl, nh4n_mgl, po4p_mgl) %>%
  full_join(met, by = c('site', 'date')) %>%
  arrange(date)%>%
  group_by(site) %>%
  mutate(across(-date, ~na.approx(., na.rm = F, maxgap = 3))) %>% ungroup() 
  
filled %>%
  # mutate(nh4n_mgl = ifelse(nh4n_mgl>0.027, NA, nh4n_mgl),
  #        no3n_mgl = ifelse(no3n_mgl>0.4, NA, no3n_mgl),
  #        tdn_mgl = ifelse(tdn_mgl>0.7, NA, tdn_mgl))%>%
  pivot_longer(cols = c(GPP, ER), values_to = 'gCm2d', names_to = 'met')%>%
  pivot_longer(cols = ends_with('mgl'), values_to = 'mgl', names_to = 'variable') %>%
ggplot(aes(x = mgl, y = gCm2d, col = site, group = interaction(site, met)))+
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(~variable, scales = 'free_x')

summary(lm(ER+GPP~no3n_mgl, data = filled))
summary(lm(GPP+ER~doc_mgl, data = filled))
summary(lm(ER~doc_mgl, data = filled))
summary(lm(GPP~DO_mgl, data = filled))
summary(lm(no3n_mgl~DO_mgl, data = filled))
summary(lm(doc_mgl~DO_mgl, data = filled))

met <- readRDS("C:/Users/Alice Carter/git/nhc_50yl/NHC_2019_metabolism/data/metabolism/compiled/met_preds_stream_metabolizer.rds")$preds %>%
  filter(era == 'now', site !='PWC', year != 2020) %>%
  select(year, date, site, GPP, ER, discharge, DO_mgl = DO.obs, DO.sat) %>%
  mutate(DO_mgl = DO_mgl/DO.sat)

DOfits <- data.frame()
for(s in unique(met$site)){
  dd <- met %>% filter(site == s)
  for(y in unique(dd$year)){
    yy <- dd %>% filter(year == y)
    m = summary(lm(-GPP-ER~DO_mgl, data = yy))
    rr <- data.frame(year = y, site = s,
      slope = m$coefficients[2,1],
      rsq = m$r.squared,
      p = m$coefficients[2,4])
    yy <- yy %>% filter(substr(date, 6,7) %in% c('10','11'))
    m = summary(lm(-GPP-ER~DO_mgl, data = yy))
    rr$slope_oct = m$coefficients[2,1]
    rr$rsq_oct = m$r.squared
    rr$p_oct = m$coefficients[2,4]
    
    DOfits <- bind_rows(DOfits, rr)
  }
}
summary(DOfits)
  
dd <- met %>% filter(site %in% c('NHC', 'UNHC')) 
s = 'both'  

