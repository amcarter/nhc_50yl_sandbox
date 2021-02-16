library(tidyverse)
library(lubridate)

setwd('~/git/papers/alice_nhc')
source('src/helpers.R')

d = read_csv('~/git/papers/alice_nhc/data/compiled_nhc_dat.csv')

d = d %>%
    select(DateTime_UTC, site, DO.obs, discharge) %>%
    pivot_wider(names_from = site,
                values_from = -all_of(c('DateTime_UTC', 'site'))) %>%
    select(DateTime_UTC, starts_with('DO.obs'), discharge_NHC)

d_plotinds = which(d$DateTime_UTC > as.POSIXct('2020-03-01', tz = 'UTC') &
                       d$DateTime_UTC < as.POSIXct('2020-04-01', tz = 'UTC'))
par(mfrow=c(4, 2), mar=c(0,0,0,0), oma=c(0,0,0,0))
for(dfc in colnames(d)[-1]){
    # plot(d$DateTime_UTC, d[, dfc, drop = TRUE],
    plot(d$DateTime_UTC[d_plotinds], d[d_plotinds, dfc, drop = TRUE],
         col='orange', type='p', lwd=2, pch = '.')
    # abline(v=d$DateTime_UTC[substr(d$DateTime_UTC, 12, 19) == '00:00:00'],
    #     lty=3, col='gray30')
    mtext(dfc, 3, line=-4)
}

#fill in any missing records
samp_int_m = determine_sample_interval(d)
d = populate_missing_rows(d, samp_int=samp_int_m)
# dtcol = d$DateTime_UTC
# d$DateTime_UTC = NULL

#keep track of NAs
na_inds = lapply(select(d, -DateTime_UTC),
                 function(x) which(is.na(x)))

# pldf = lin_interp_gaps(pldf, samp_int=samp_int_m, gap_thresh=180)
# na_inds = lapply(d, function(x) which(is.na(x)))

# d = bind_cols(tibble('DateTime_UTC'=dtcol), d)

#Nearest Days Interpolation
# ndiout = NDI(d, target_variable = 'DO.obs_UNHC', interv=samp_int_m)
ndiout = NDI(d, interv=samp_int_m)

#unhc has the march 7-19 gap
#all but nhc and unhc need to be extrapolated
ndi_plotinds = which(ndiout$DateTime_UTC > as.POSIXct('2020-03-01', tz = 'UTC') &
          ndiout$DateTime_UTC < as.POSIXct('2020-04-01', tz = 'UTC'))
d_plotinds = which(ndiout$DateTime_UTC > as.POSIXct('2020-03-01', tz = 'UTC') &
          ndiout$DateTime_UTC < as.POSIXct('2020-04-01', tz = 'UTC'))
plot(d$DateTime_UTC[ndi_plotinds], d$DO.obs_UNHC[ndi_plotinds], col='blue', type='b')
lines(ndiout$DateTime_UTC[d_plotinds], ndiout$DO.obs_UNHC[d_plotinds], col='red')
# plot(d$DateTime_UTC, d$DO.obs_UNHC, type='n')
# lines(d$DateTime_UTC, d$DO.obs_PM)

#assign qaqc code 2 to any gaps filled by NDI
dfcols = colnames(ndiout)

par(mfrow=c(4, 2), mar=c(0,0,0,0), oma=c(0,0,0,0))
for(dfc in dfcols[-1]){
    plot(ndiout$DateTime_UTC, ndiout[,dfc], col='orange', type='l', lwd=2)
    lines(d$DateTime_UTC, d[, dfc, drop=TRUE], lwd=2)
    # abline(v=d$DateTime_UTC[substr(d$DateTime_UTC, 12, 19) == '00:00:00'],
    #     lty=3, col='gray30')
    mtext(dfc, 3, line=-4)
}

# par(mfrow=c(1,1), mar=c(4,4,4,4))
# plot(ndiout$DateTime_UTC, ndiout[[c]], col='orange', type='l', lwd=2)
# lines(pldf$DateTime_UTC, pldf[[c]], lwd=2)
# # abline(v=pldf$DateTime_UTC[substr(pldf$DateTime_UTC, 12, 19) == '00:00:00'],
# #     lty=3, col='gray30')
# # for(i in 1:nrow(fullday_skips)){
# for(i in 1:2){
#     abline(v=pldf$DateTime_UTC[fullday_skips[i,2:3]], col='gray')
#     print(paste0(i, ': ', fullday_skips[i,2], '-', fullday_skips[i,3],
#         '; ', fullday_skips[i,3] - fullday_skips[i,2]))
#     # readline()
# }
# # abline(v=pldf$DateTime_UTC[fullday_skips[,2]], col='gray')
#
# xlims = as.numeric(as.POSIXct(c('2006-10-04', '2007-1-06')))
# plot(ndiout$DateTime_UTC, ndiout[[c]], col='orange', type='l',
#     lwd=2, xlim=xlims)
# lines(pldf$DateTime_UTC, pldf[[c]], lwd=2)
# abline(v=pldf$DateTime_UTC[substr(pldf$DateTime_UTC, 12, 19) == '00:00:00'],
#     lty=3, col='gray30')
# for(i in 1:2){
#     abline(v=pldf$DateTime_UTC[fullday_skips[i,2:3]], col='gray')
#     print(paste0(i, ': ', fullday_skips[i,2], '-', fullday_skips[i,3],
#         '; ', fullday_skips[i,3] - fullday_skips[i,2]))
# }

# plot(ndiout$DateTime_UTC, ndiout$AirPres_kPa, col='orange', type='l',
#     lwd=2, ylim=c(14.84, 14.9))
# lines(pldf$DateTime_UTC, pldf$AirPres_kPa, lwd=2)
ndiout = snap_days(dfcols, flagdf, ndiout, pldf, interv=samp_int_m)

# lines(ndiout$DateTime_UTC, ndiout[[c]], col='red')
# lines(pldf$DateTime_UTC, pldf[[c]])



#remove entirely empty rows?
