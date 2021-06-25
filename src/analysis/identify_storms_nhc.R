#########################################################################
# Identify storms in NHC through baseflow separation or through % change in flow

library(tidyverse)
library(lubridate)
setwd(hypox_projdir)

dat <- read_csv("data/raw/NHCdat.csv")
NHC <- dat %>%
  mutate(datetime = with_tz(DateTime_UTC, "EST"),
         year = year(datetime)) %>%
  filter(year == 2019) %>%
  select(datetime, discharge_cms)

stormdates <- as.Date(c("2019-01-13","2019-01-20","2019-01-24","2019-02-13",
                        "2019-02-16","2019-02-18","2019-02-21","2019-02-24",
                        "2019-03-21","2019-04-06","2019-04-09","2019-04-13",
                        "2019-04-20","2019-06-08","2019-06-19","2019-07-23",
                        "2019-08-05","2019-08-14","2019-10-20","2019-11-24",
                        "2019-12-01","2019-12-14", "2019-12-24"))
plot(NHC$datetime, NHC$discharge_cms, log = "y", type = "l")
abline(v = as.POSIXct(stormdates), col = "green")

NHC <- NHC%>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date) %>%
  summarize(minq = min(discharge_cms, na.rm = T),
            maxq = max(discharge_cms, na.rm = T)) %>%
  mutate(minq = ifelse(is.infinite(minq), NA_real_, minq),
         maxq = ifelse(is.infinite(maxq), NA_real_, maxq))
NHC$jump <- NHC$maxq/c(NA_real_, NHC$minq[1:(nrow(NHC)-1)]) 
plot(NHC$date, NHC$jump, log = "y")
abline(v = stormdates, col = "green")
