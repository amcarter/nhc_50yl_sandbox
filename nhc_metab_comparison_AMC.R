#################################
# Metabolism comparison NHC SP data and Hall data

#questions:
#combine data from historic sites? NHC UNHC %in% Concrete, Blackwood, Wood Bridge?
#combine years?

library(StreamPULSE)
library(viridis)
library(ggplot2)
library(beanplot)
library(dplyr)
library(scales)

setwd('C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl')


#setup ####

#switch this to TRUE if you want to use only modern metab estimates from days
#in which Q is in the same ballpark as it was for HALL's K estimates
filter_high_Q = FALSE

#read in historic data and average across multiple same-site, same-day estimates
nhc_68_70 = read.csv('code/mike_data/hall_table_15.csv', colClasses=c('date'='Date'))
nhc_68_70 = nhc_68_70 %>%
    group_by(date, site) %>%
    summarize_if(is.numeric, mean, na.rm=TRUE) %>%
    as.data.frame()

#subset historic data by site and year
gpp_concrete = nhc_68_70$GPP_gO2m2d[nhc_68_70$site == 'Concrete']
gpp_blackwood = nhc_68_70$GPP_gO2m2d[nhc_68_70$site == 'Blackwood']
gpp_wb = nhc_68_70$GPP_gO2m2d[nhc_68_70$site == 'Wood Bridge']

gpp_68_70 = nhc_68_70$GPP_gO2m2d
gpp_68 = nhc_68_70$GPP_gO2m2d[substr(nhc_68_70$date, 1, 4) == '1968']
gpp_69 = nhc_68_70$GPP_gO2m2d[substr(nhc_68_70$date, 1, 4) == '1969']
gpp_70 = nhc_68_70$GPP_gO2m2d[substr(nhc_68_70$date, 1, 4) == '1970']

er_68_70 = nhc_68_70$ER_gO2m2d
er_68 = nhc_68_70$ER_gO2m2d[substr(nhc_68_70$date, 1, 4) == '1968']
er_69 = nhc_68_70$ER_gO2m2d[substr(nhc_68_70$date, 1, 4) == '1969']
er_70 = nhc_68_70$ER_gO2m2d[substr(nhc_68_70$date, 1, 4) == '1970']

nep_68_70 = gpp_68_70 - er_68_70
nep_68 = gpp_68 - er_68
nep_69 = gpp_69 - er_69
nep_70 = gpp_70 - er_70

#retrieve contemporary data by year; get K and O2 data for later
# query_available_results('NC', 'NHC')
# nhc_17 = request_results(sitecode='NC_NHC', '2017')
# nhc_18 = request_results(sitecode='NC_NHC', '2018')
# nhc_new_K = nhc_17$model_results$data %>%
#     select(date, DO.mod, DO.sat, temp.water, depth) %>%
#     # bind_rows(nhc_18$model_results$data[,c('date','DO.obs','DO.sat')]) %>%
#     bind_rows(select(nhc_18$model_results$data, date, DO.mod, DO.sat,
#         temp.water, depth)) %>%
#     # group_by(date) %>%
#     # summarize_all(mean) %>%
#     # left_join(nhc_17$model_results$fit$daily[,c('date','K600_daily_mean')]) %>%
#     left_join(select(bind_rows(nhc_17$model_results$fit$daily,
#         nhc_18$model_results$fit$daily), date, K600_daily_mean)) %>%
#     as.data.frame()

nhc_new <- read_rds("data/NHC_metab_allsites.rds")

# if(filter_high_Q){
# 
#     #highest considered depth in Hall dissertation: 0.65m
#     #corresponding Q, based on modern Z-Q curve:
#     d17 = nhc_17$model_results$data
#     dd17 = nhc_17$model_results$data_daily
#     Q_cutoff = max(na.omit(d17$discharge[d17$depth > 0.64 & d17$depth < 0.65]))
#     Q_cutoff = max(na.omit(d17$discharge[d17$depth > 0.84 & d17$depth < 0.85]))
#     Qbool = dd17$discharge.daily < Q_cutoff
#     Qbool[is.na(Qbool)] = FALSE
#     dd17 = dd17[Qbool,]
#     pred17 = nhc_17$predictions[,c('date','GPP','ER')]
#     nhc_17 = inner_join(pred17, dd17[,c('date', 'discharge.daily')], by='date')
# 
#     d18 = nhc_18$model_results$data
#     dd18 = nhc_18$model_results$data_daily
#     # Q_cutoff = max(na.omit(d18$discharge[d18$depth > 0.64 & d18$depth < 0.65]))
#     Qbool = dd18$discharge.daily < Q_cutoff
#     Qbool[is.na(Qbool)] = FALSE
#     dd18 = dd18[Qbool,]
#     pred18 = nhc_18$predictions[,c('date','GPP','ER')]
#     nhc_18 = inner_join(pred18, dd18[,c('date', 'discharge.daily')], by='date')
# } else {
#     nhc_17 = as.data.frame(nhc_17$predictions[,c('date','GPP','ER')])
#     nhc_18 = as.data.frame(nhc_18$predictions[,c('date','GPP','ER')])
# }


gpp_new = nhc_new$GPP

er_new = nhc_new$ER

nep_new = gpp_new + er_new

dates_new = nhc_new$date

#time-series comparison of means then and now? ####
plot(gpp_new, type='l')
acf(gpp_new, na.action=na.pass)
pacf(gpp_new, na.action=na.pass)
#strong autocorrelation and partial autocorr;
#will have to model error as an autoregressive process
#not stationary; can't use pure AR
qqnorm(gpp_new); abline(1, 1, col='red', lty=2)
#normalish; no need to go bayesian

#nhc_17 is irregular, so can't use arima; only option would be GAM


#dist plots ####
png(width=9, height=6, units='in', type='cairo', res=300,
    filename='~/Dropbox/streampulse/figs/NHC_comparison/metab_distributions.png')

defpar = par(mfrow=c(2,3))

#plot GPP dists, then and now
plot(density(gpp_68_70, na.rm=TRUE), xlim=c(-3, 10), bty='l', col='sienna3',
    main='GPP 1968-70 vs. 2019', xlab='GPP', ylim=c(0,1.2))
lines(density(gpp_new, na.rm=TRUE), col='blue')
legend('topright', legend=c('68-70; n=79', paste0('19; n=', length(gpp_new))),
    col=c('sienna3','blue'), lty=1, bty='n', seg.len=1, cex=0.9, lwd=2)

#plot ER dists, then and now
plot(density(er_68_70 * -1, na.rm=TRUE), xlim=c(-15, 1), bty='l', col='sienna3',
    main='ER 1968-70 vs. 2019', xlab='ER', ylim=c(0,0.7))
lines(density(er_new, na.rm=TRUE), col='blue')
legend('topleft', legend=c('68-70; n=79', paste0('19; n=', length(gpp_new))),
    col=c('sienna3','blue'), lty=1, bty='n', seg.len=1, cex=0.9, lwd=2)

#plot NEP dists, then and now
plot(density(nep_68_70, na.rm=TRUE), xlim=c(-15, 2), bty='l', col='sienna3',
    main='NEP 1968-70 vs. 2019', xlab='NEP', ylim=c(0,1.0))
lines(density(nep_new, na.rm=TRUE), col='blue')
legend('topleft', legend=c('68-70; n=79', paste0('19-20; n=', length(nrow(gpp_new)))),
    col=c('sienna3','blue'), lty=1, bty='n', seg.len=1, cex=0.9, lwd=2)

#plot GPP dists by year
cols = viridis(4)
plot(density(gpp_68, na.rm=TRUE), xlim=c(-3, 10), bty='l', col=cols[1],
    main='GPP by year', xlab='GPP', ylim=c(0,1.3))
lines(density(gpp_69, na.rm=TRUE), col=cols[2])
lines(density(gpp_70, na.rm=TRUE), col=cols[3])
lines(density(gpp_new, na.rm=TRUE), col=cols[4])
legend('topright',
    legend=c('68; n=18', '69; n=49', '70; n=12', '19; n=2403'),
    col=cols, lty=1, bty='n', seg.len=1, cex=0.9, lwd=2)

#plot ER dists by year
cols = viridis(4)
plot(density(er_68 * -1, na.rm=TRUE), xlim=c(-18, 2), bty='l', col=cols[1],
    main='ER by year', xlab='ER', ylim=c(0,0.7))
lines(density(er_69 * -1, na.rm=TRUE), col=cols[2])
lines(density(er_70 * -1, na.rm=TRUE), col=cols[3])
lines(density(er_new, na.rm=TRUE), col=cols[4])
legend('topleft',
    legend=c('68; n=18', '69; n=49', '70; n=12', '19; n=2403'),
    col=cols, lty=1, bty='n', seg.len=1, cex=0.9, lwd=2)

#plot NEP dists by year
cols = viridis(4)
plot(density(nep_68, na.rm=TRUE), xlim=c(-18, 2), bty='l', col=cols[1],
    main='NEP by year', xlab='NEP', ylim=c(0,1.1))
lines(density(nep_69, na.rm=TRUE), col=cols[2])
lines(density(nep_70, na.rm=TRUE), col=cols[3])
lines(density(nep_new, na.rm=TRUE), col=cols[4])
legend('topleft',
    legend=c('68; n=18', '69; n=49', '70; n=12', '19; n=2403'),
    col=cols, lty=1, bty='n', seg.len=1, cex=0.9, lwd=2)

dev.off()

#plot temporal coverage for historic data ####

historic_dates = nhc_68_70$date[! is.na(nhc_68_70$GPP_gO2m2d)]
historic_year_agg = as.character(historic_dates)
substr(historic_year_agg, 1, 4) = '1970'
historic_year_agg = as.Date(historic_year_agg)
hy_num = as.numeric(historic_year_agg)

hd_concrete = nhc_68_70$date[! is.na(nhc_68_70$GPP_gO2m2d) &
        nhc_68_70$site == 'Concrete']
concrete_year_agg = as.character(hd_concrete)
substr(concrete_year_agg, 1, 4) = '1970'
concrete_year_agg = as.Date(concrete_year_agg)
concrete_num = as.numeric(concrete_year_agg)

hd_blackwood = nhc_68_70$date[! is.na(nhc_68_70$GPP_gO2m2d) &
        nhc_68_70$site == 'Blackwood']
blackwood_year_agg = as.character(hd_blackwood)
substr(blackwood_year_agg, 1, 4) = '1970'
blackwood_year_agg = as.Date(blackwood_year_agg)
blackwood_num = as.numeric(blackwood_year_agg)

hd_wb = nhc_68_70$date[! is.na(nhc_68_70$GPP_gO2m2d) &
        nhc_68_70$site == 'Wood Bridge']
wb_year_agg = as.character(hd_wb)
substr(wb_year_agg, 1, 4) = '1970'
wb_year_agg = as.Date(wb_year_agg)
wb_num = as.numeric(wb_year_agg)

# plot(historic_dates, rep(1, length(historic_dates), type='n', xlab='day'),
#     yaxt='n', ylab='', xlab='', main='Historic Coverage')
# abline(v=historic_dates, lty=2, col='gray')

png(width=7, height=6, units='in', type='cairo', res=300,
    filename='~/Dropbox/streampulse/figs/NHC_comparison/historic_coverage.png')

par(mfrow=c(4,1), mar=c(0,0,0,0), oma=c(3,4,3,0))

beanplot(hy_num, horizontal=TRUE, col='gray', xaxt='n',
    frame.plot=FALSE, ylim=lims)
    # main='Historic Annual Coverage Across Sites')
# axis(1, at=seq(as.Date('1970-01-01'), as.Date('1970-12-31'), length.out=13)[1:12],
#     labels=month.abb)
mtext('All sites', 2)
mtext('Historic Annual Coverage', 3)
legend('topright', legend=paste('n =', length(! is.na(hy_num))),
    bty='n', cex=1.3, text.font=2)

#plot temporal coverage by site
lims = c(min(hy_num), max(hy_num))
beanplot(concrete_num, horizontal=TRUE, col='yellow', xaxt='n',
    frame.plot=FALSE, ylim=lims)
mtext('Concrete', 2)
legend('topright', legend=paste('n =', length(! is.na(concrete_num))),
    bty='n', cex=1.3, text.font=2)

beanplot(blackwood_num, horizontal=TRUE, col='green', xaxt='n',
    frame.plot=FALSE, ylim=lims)
# legend('left', legend=c('Concrete', 'Blackwood', 'Wood Bridge'),
#     fill=c('yellow', 'green', 'orange'), cex=2, bty='n')
mtext('Blackwood', 2)
legend('topright', legend=paste('n =', length(! is.na(blackwood_num))),
    bty='n', cex=1.3, text.font=2)

beanplot(wb_num, horizontal=TRUE, col='orange', xaxt='n',
    frame.plot=FALSE, ylim=lims)
mtext('Wood Bridge', 2)
legend('topright', legend=paste('n =', length(! is.na(wb_num))),
    bty='n', cex=1.3, text.font=2)

# ax_dt = as.numeric(historic_dates)
# ax_seq = seq(ax_dt[1], ax_dt[which.max(ax_dt)], length.out=10)
# axis(1, at=ax_seq, labels=as.Date(ax_seq), las=1, cex.axis=1.7)
axis(1, at=seq(as.Date('1970-01-01'), as.Date('1970-12-31'), length.out=13)[1:12],
    labels=month.abb)
# mtext(text='Historic Annual Coverage by Site', side=3, outer=TRUE, cex=1.8)

dev.off()

#compare overall distributions then and now ####

#first assess normality
png(width=7, height=6, units='in', type='cairo', res=300,
    filename='~/Dropbox/streampulse/figs/NHC_comparison/normality_assessment.png')

par(mfrow=c(2,2), mar=c(0, 0, 0, 0), oma=rep(4, 4))
qqnorm(gpp_new, lty=2, xlab='', ylab='', main='', xaxt='n', yaxt='n', bty='o')
abline(0, 1, col='red', lty=2)
legend('bottom', legend='GPP 2017-18', bty='n', cex=1.3)
qqnorm(er_new, lty=2, xlab='', ylab='', main='', xaxt='n', yaxt='n', bty='o')
abline(0, 1, col='red', lty=2)
legend('bottom', legend='ER 2017-18', bty='n', cex=1.3)
qqnorm(gpp_68_70, lty=2, xlab='', ylab='', main='', xaxt='n', yaxt='n', bty='o')
abline(0, 1, col='red', lty=2)
legend('top', legend='GPP 1968-70', bty='n', cex=1.3)
qqnorm(er_68_70, lty=2, xlab='', ylab='', main='', xaxt='n', yaxt='n', bty='o')
abline(0, 1, col='red', lty=2)
legend('top', legend='ER 1968-70', bty='n', cex=1.3)
mtext('Theoretical Quantiles', 1, outer=TRUE, line=1.5)
mtext('Sample Quantiles', 2, outer=TRUE, line=1.5)
mtext('Normal Q-Q Plots (red line = 1:1)', 3, outer=TRUE, line=1.5)

dev.off()

#nonnormal, but CLT probably applies.
#let's assess equality of variance with an F-test
var.test(gpp_68_70, gpp_new) #not equal: p < 0.001
var.test(er_68_70, er_new) #not equal: p < 0.001

#unequal variance, so 2-sample t-test is out.
#not i.i.d., so welch's t-test is out (could transform)
#can't do Mann-Whitney-Wilcoxon Test either because of unequal var and autocorr

#bootstrap 2-samp t-test for GPP (and we'll go with Welch's) ####

#get observed t-statistic
t_obs_gpp = t.test(gpp_68_70, gpp_new, var.equal=FALSE)$statistic

#artificially make both sample means identical (satisfy the null)
gpp_68_70_mod = gpp_68_70 - mean(gpp_68_70, na.rm=TRUE) +
    mean(c(gpp_68_70, gpp_new), na.rm=TRUE)
gpp_new_mod = gpp_new - mean(gpp_new, na.rm=TRUE) +
    mean(c(gpp_68_70, gpp_new), na.rm=TRUE)

#verify
round(mean(gpp_68_70_mod, na.rm=TRUE), 7) ==
    round(mean(gpp_new_mod, na.rm=TRUE), 7)

#get historic monthly data coverage to use for sample weights
month_counts_68_70 = tapply(rep(1, nrow(nhc_68_70)),
    substr(nhc_68_70$date, 6, 7), sum)
month_proportions = month_counts_68_70 / sum(month_counts_68_70)

#split gpp vector by month for each dataset
#commented portions are for uniformly distributing monthly draw weights
gpp_68_70_bymo = split(gpp_68_70_mod, factor(substr(nhc_68_70$date, 6, 7)))
# gpp_68_70_bymo = split(gpp_68_70_mod,
#     factor(rep(c('01','02','03','04','05','06','07','08','10','11','12'),
#     length.out=length(gpp_68_70_mod))))
gpp_new_bymo = split(gpp_new_mod, factor(substr(dates_new, 6, 7)))
# gpp_new_bymo = split(gpp_new_mod,
#     factor(rep(c('01','02','03','04','05','06','07','08','10','11','12'),
#     length.out=length(gpp_new_mod))))
gpp_new_bymo = lapply(gpp_new_bymo, na.omit)
nsamp_new = sum(sapply(gpp_new_bymo, length))
nsamp_68_70 = length(gpp_68_70_mod)


#determine monthly sample sizes for modern dataset; deal with remainders
month_samp_new = month_proportions * nsamp_new
extra_sample_probs = month_samp_new %% 1
month_samp_new = floor(month_samp_new)

#get bootstrap estimate of sampling distribution of the t-stat if H0 is true;
#i.e. bootstrap the null distribution (weight draws by historic monthly coverage)
nsamp = 20000

t_vect_gpp = vector(length=nsamp)
for(i in 1:nsamp){
    samp_68_70_gpp = samp_new_gpp = c()
    remainder_months = sample(c(1:8, 10:12), size=sum(extra_sample_probs),
        prob=extra_sample_probs)
    for(j in c(1:8, 10:12)){
        extra_samp = ifelse(j %in% remainder_months, 1, 0)
        j = sprintf('%02d', j)
        samp_68_70_gpp = append(samp_68_70_gpp, sample(gpp_68_70_bymo[[j]],
            size=month_counts_68_70[j], replace=TRUE))
        samp_new_gpp = append(samp_new_gpp, sample(gpp_new_bymo[[j]],
            size=month_samp_new[j] + extra_samp, replace=TRUE))
    }
    # samp_68_70_gpp = sample(gpp_68_70_mod, size=nsamp_68_70, replace=TRUE)
    # samp_new_gpp = sample(gpp_new_mod, size=nsamp_new, replace=TRUE)
    t_vect_gpp[i] = t.test(samp_68_70_gpp, samp_new_gpp,
        var.equal=FALSE)$statistic
}

#p-val is proportion of times observed t-statistic >= bootstrap t-statistic
pval_gpp = (sum(t_vect_gpp <= t_obs_gpp) + 1) / (nsamp + 1)
if(pval_gpp == 1){
    pval_gpp = (sum(t_vect_gpp >= t_obs_gpp) + 1) / (nsamp + 1)
}

#bootstrap Welch's t-test for ER ####

#get observed t-statistic
t_obs_er = t.test(er_68_70, er_new, var.equal=FALSE)$statistic

#artificially make both sample means identical (satisfy the null)
er_68_70_mod = er_68_70 - mean(er_68_70, na.rm=TRUE) +
    mean(c(er_68_70, er_new), na.rm=TRUE)
er_new_mod = er_new - mean(er_new, na.rm=TRUE) +
    mean(c(er_68_70, er_new), na.rm=TRUE)

#verify
mean(er_68_70_mod, na.rm=TRUE) == mean(er_new_mod, na.rm=TRUE)

#split er vector by month for each dataset
er_68_70_bymo = split(er_68_70_mod, factor(substr(nhc_68_70$date, 6, 7)))
er_new_bymo = split(er_new_mod, factor(substr(dates_new, 6, 7)))
er_new_bymo = lapply(er_new_bymo, na.omit)
nsamp_new = sum(sapply(er_new_bymo, length))
nsamp_68_70 = length(er_68_70_mod)

#determine monthly sample sizes for modern dataset; deal with remainders
month_samp_new = month_proportions * nsamp_new
extra_sample_probs = month_samp_new %% 1
month_samp_new = floor(month_samp_new)

#get bootstrap estimate of sampling distribution of the t-stat if H0 is true;
#i.e. bootstrap the null distribution (weight draws by historic monthly coverage)
nsamp = 20000
t_vect_er = vector(length=nsamp)
for(i in 1:nsamp){
    samp_68_70_er = samp_new_er = c()
    remainder_months = sample(c(1:8, 10:12), size=sum(extra_sample_probs),
        prob=extra_sample_probs)
    for(j in c(1:8, 10:12)){
        extra_samp = ifelse(j %in% remainder_months, 1, 0)
        j = sprintf('%02d', j)
        samp_68_70_er = append(samp_68_70_er, sample(er_68_70_bymo[[j]],
            size=month_counts_68_70[j], replace=TRUE))
        samp_new_er = append(samp_new_er, sample(er_new_bymo[[j]],
            size=month_samp_new[j] + extra_samp, replace=TRUE))
    }
    # samp_68_70_er = sample(er_68_70_mod, size=nsamp_68_70, replace=TRUE)
    # samp_new_er = sample(er_new_mod, size=nsamp_new, replace=TRUE)
    t_vect_er[i] = t.test(samp_68_70_er, samp_new_er,
        var.equal=FALSE)$statistic
}

#p-val is proportion of times observed t-statistic >= bootstrap t-statistic
pval_er = (sum(t_vect_er >= t_obs_er) + 1) / (nsamp + 1)
if(pval_er == 1){
    pval_er = (sum(t_vect_er >= t_obs_er) + 1) / (nsamp + 1)
}
#visualize GPP hypothesis test ####

png(width=7, height=6, units='in', type='cairo', res=300,
    filename='~/Dropbox/streampulse/figs/NHC_comparison/bootstrap_welch_t_weighted_filtered.png')

par(mfrow=c(2,1), mar=c(4,4,1,2), oma=c(0,0,3,0))

plot(density(t_vect_gpp), xlab='t-value', main='', xlim = c(-8,4))
qs = quantile(t_vect_gpp, probs=c(0.025, 0.975))
dd = density(t_vect_gpp)
ddo = order(dd$x)
xdens = dd$x[ddo]
ydens = dd$y[ddo]
xdens_lt = xdens[xdens <= qs[1]]
ydens_lt = ydens[xdens <= qs[1]]
polygon(c(xdens_lt, rev(xdens_lt)), c(ydens_lt, rep(0, length(ydens_lt))),
    col='lightgreen', border='lightgreen')
xdens_ut = xdens[xdens >= qs[2]]
ydens_ut = ydens[xdens >= qs[2]]
polygon(c(xdens_ut, rev(xdens_ut)), c(ydens_ut, rep(0, length(ydens_ut))),
    col='lightgreen', border='lightgreen')
abline(v=t_obs_gpp, lty=2, col='red', lwd=2)
legend('topleft', legend='GPP', bty='n', text.font=2, cex=1)
legend('topleft', legend=paste('\np =', round(pval_gpp, 3)), bty='n',
    text.font=1, cex=1)

#why bimodality in the null dist?
#1. not from skewed historic draw weights; artificially uniformified them to test.
#2. is from skewed modern draw weights; artificially uniformified them to test.
#3. there is multimodality in some of the modern monthly GPP dists.
#if most draws come from one GPP peak or another,
#the t-val may land in one H0 peak or another

#visualize ER hypothesis test ####

plot(density(t_vect_er), xlim=c(-5, 40), xlab='t-value', main='')
qs = quantile(t_vect_er, probs=c(0.025, 0.975))
dd = density(t_vect_er)
ddo = order(dd$x)
xdens = dd$x[ddo]
ydens = dd$y[ddo]
xdens_lt = xdens[xdens <= qs[1]]
ydens_lt = ydens[xdens <= qs[1]]
polygon(c(xdens_lt, rev(xdens_lt)), c(ydens_lt, rep(0, length(ydens_lt))),
    col='lightgreen', border='lightgreen')
xdens_ut = xdens[xdens >= qs[2]]
ydens_ut = ydens[xdens >= qs[2]]
polygon(c(xdens_ut, rev(xdens_ut)), c(ydens_ut, rep(0, length(ydens_ut))),
    col='lightgreen', border='lightgreen')
abline(v=t_obs_er, lty=2, col='red', lwd=2)
legend('top', legend='ER', bty='n', text.font=2, cex=1)
legend('top', legend=paste('\np =', round(pval_er, 3)), bty='n',
    text.font=1, cex=1)

mtext("Observed Welch's t-values (red) relative to bootstrapped null dists", 3,
    outer=TRUE, line=1, font=2, cex=1.3)

dev.off()

#verify with Mann-Whitney-Wilcoxon Test? ####

#visualize dists again with non-bootstrapped means ####

png(width=7, height=6, units='in', type='cairo', res=300,
    filename='~/Dropbox/streampulse/figs/NHC_comparison/means_raw_boxplot.png')

gppHmean = paste('mean =', round(mean(gpp_68_70, na.rm=TRUE), 2))
gppCmean = paste('mean =', round(mean(gpp_new, na.rm=TRUE), 2))
erHmean = paste('mean =', round(mean(er_68_70, na.rm=TRUE), 2) * -1)
erCmean = paste('mean =', round(mean(er_new, na.rm=TRUE), 2))
boxplot(gpp_68_70, gpp_new * -1, er_68_70, er_new,
    ylab='', col='gray',
    names=c('GPP 1968-70', 'GPP 2019', 'ER 1968-70', 'ER 2019'))
axis(1, at=1:4, labels=c(gppHmean, gppCmean, erHmean, erCmean),
    line=1.5, col='transparent', tcl=0, font=2)
mtext(expression(paste("gm"^"-2" * " d"^"-1")), 2, line=2)
mtext('Another look at distributions, then and now (not bootstrapped)', 3,
    cex=1, font=2)

dev.off()

#bootstrap some confidence bounds (NEEDS TO BE WEIGHTED) ####
nsamp = 20000
mean_vect_er_68_70 = mean_vect_er_new = mean_vect_gpp_68_70 =
    mean_vect_gpp_new = vector(length=nsamp)
for(i in 1:nsamp){
    samp_68_70_er = sample(er_68_70, size=79, replace=TRUE)
    samp_new_er = sample(er_new, size=730, replace=TRUE)
    samp_68_70_gpp = sample(gpp_68_70, size=79, replace=TRUE)
    samp_new_gpp = sample(gpp_new, size=730, replace=TRUE)
    mean_vect_er_68_70[i] = mean(samp_68_70_er, na.rm=TRUE)
    mean_vect_er_new[i] = mean(samp_new_er, na.rm=TRUE)
    mean_vect_gpp_68_70[i] = mean(samp_68_70_gpp, na.rm=TRUE)
    mean_vect_gpp_new[i] = mean(samp_new_gpp, na.rm=TRUE)
}

# plot(density(mean_vect_er_68_70 * -1))
CI = data.frame('CI95_lower'=numeric(4), 'median'=numeric(4),
    'CI95_upper'=numeric(4),
    row.names=c('GPP_then', 'GPP_now', 'ER_then', 'ER_now'))
CI[1,] = quantile(sort(mean_vect_gpp_68_70), probs=c(0.025, 0.5, 0.975))
CI[2,] = quantile(sort(mean_vect_gpp_new), probs=c(0.025, 0.5, 0.975))
CI[3,] = quantile(sort(mean_vect_er_68_70) * -1, probs=c(0.025, 0.5, 0.975))
CI[4,] = quantile(sort(mean_vect_er_new), probs=c(0.025, 0.5, 0.975))
write.csv(CI, '~/Dropbox/streampulse/figs/NHC_comparison/metab_CIs.csv')

#evaluate k now vs. then ####

#read in hall's diffusion constant tables (attained via 3 different means)
k_diurnal = read.csv('code/mike_data/hall_k_diurnalO2_unused.csv', 
                     colClasses=c('date'='Date'))
k_diurnal$method = 'diurnal'
k_morphology = read.csv('code/mike_data/hall_k_morphology.csv')
k_morphology$method = 'morphology'
k_dome = read.csv('code/mike_data/hall_k_dome.csv', colClasses=c('date'='Date'))
k_dome$method = 'dome'
nhc_68_70_k = bind_rows(k_diurnal, k_morphology, k_dome)

nhc_68_70_k$k_daily = nhc_68_70_k$k * 24

#convert modern K600 to K2 (1/day)
SA = 1568
SB = -86.04
SC = 2.142
SD = -0.0216
SE = -0.5
if(filter_high_Q){ #prob should do this anyway
    nhc_new_K = nhc_new_K[nhc_new_K$date %in% dates_new,]
}
TT = nhc_new_K$temp.water
nhc_new_K$K2 = nhc_new_K$K600_daily_mean *
    ((SA + SB * TT + SC * TT ^ 2 + SD * TT ^ 3) / 600) ^ SE

# #convert K2 to D
# nhc_new_K$D = nhc_new_K$K2 * (nhc_new_K$DO.sat - nhc_new_K$DO.mod)
# D_daily = tapply(nhc_new_K$D, nhc_new_K$date, mean, na.rm=TRUE)
# # K2_daily = tapply(nhc_new_K$K2, nhc_new_K$date, mean, na.rm=TRUE)


#calculate 100% saturation deficit expressed as g m-3 (same as mg/L)
#S = saturation deficit = 100 - % saturation of stream
#can be expressed in terms of % or concentration
#but what is "100% saturation deficit"? is that when the stream has 0 O2?
#in that case, 100% saturation deficit = 100, or expressed as conc, = DO sat
#moving on...

#convert 15m K2 (sm parlance; would be k2 in hall parlance) to daily k (a la hall)
nhc_new_K$k = nhc_new_K$K2 * 2.3 * nhc_new_K$DO.sat #eq 9 in hall
k_daily = tapply(nhc_new_K$k, nhc_new_K$date, mean, na.rm=TRUE)

#plot distributions of historic and modern k
png(width=7, height=6, units='in', type='cairo', res=300,
    filename='~/Dropbox/streampulse/figs/NHC_comparison/k_dists_filtered.png')

xlims = range(c(k_daily, nhc_68_70_k$k_daily), na.rm=TRUE)
cur_dens = density(k_daily, na.rm=TRUE)
hist_dens_diurnal = density(nhc_68_70_k$k_daily[nhc_68_70_k$method == 'diurnal'],
    na.rm=TRUE)
hist_dens_morph = density(nhc_68_70_k$k_daily[nhc_68_70_k$method == 'morphology'],
    na.rm=TRUE)
hist_dens_dome = density(nhc_68_70_k$k_daily[nhc_68_70_k$method == 'dome'],
    na.rm=TRUE)
# ylims=range(c(max(hist_dens_diurnal$y), max(hist_dens_morph$y),
#     max(hist_dens_dome$y)), na.rm=TRUE)
plot(cur_dens, xlim=xlims, bty='l', col='sienna3',
    main='k distributions, then vs. now', yaxt='n', ylab='',
    xlab=expression(paste("gm"^"-3" * " d"^"-1")), lwd=2)
par(new=TRUE)
plot(hist_dens_diurnal, col='blue', main='', xlab='', ylab='', yaxt='n',
    xlim=xlims, bty='l', lwd=2)
par(new=TRUE)
plot(hist_dens_morph, col='darkgreen', main='', xlab='', ylab='', yaxt='n',
    xlim=xlims, bty='l', lwd=2)
par(new=TRUE)
plot(hist_dens_dome, col='purple', main='', xlab='', ylab='', yaxt='n',
    xlim=xlims, bty='l', lwd=2)
mtext('Density', 2, line=1)
legend('topright', legend=c('then (diurnal); n=9', 'then (morphology); n=14',
    'then (dome); n=5', 'now; n=475'), lty=1, bty='n', seg.len=1, cex=0.9, lwd=2,
    col=c('blue','darkgreen','purple','sienna3'))

dev.off()

#compare K by depth
k_by_depth = tapply(nhc_new_K$k, round(nhc_new_K$depth, 3), mean, na.rm=TRUE)
k_sub = k_by_depth[names(k_by_depth) %in% as.character(k_morphology$depth)]
still_need = k_morphology$depth[! as.character(k_morphology$depth) %in% names(k_sub)]
sort(unique(substr(names(k_by_depth), 1, 5)))
k_sub = c(k_by_depth[names(k_by_depth) == '0.413'], k_sub)

png(width=7, height=6, units='in', type='cairo', res=300,
    filename='~/Dropbox/streampulse/figs/NHC_comparison/k_by_depth.png')

#labels haven't been updated here, nor has anything been reconsidered since
#rethinking the conversions
cexes = scales::rescale(as.numeric(names(k_sub)), to=c(1,2))
plot(k_morphology$k[as.numeric(k_morphology$depth) >= 0.4], k_sub,
    xlab='K2 (1/day; historic, estimated analytically)',
    ylab='K2 (1/day; modern, modeled)',
    main='K2 paired by depth (seq 0.4 - 0.65m by 0.05)',
    cex=cexes)
abline(0, 1, lty=2, col='gray30')
legend('topright', legend=c('0.40 m', '0.65 m', ''), pt.cex=range(cexes), pch=1,
    col=c('black', 'black', 'transparent'))
legend('topright', legend=c('', '', '1:1'), lty=2, col=c('transparent', 'transparent',
    'gray30'), bg='transparent', bty='n')

dev.off()
