library(tidyr)
library(dplyr)

# 1. generating workable datasets for NHC
#(skip this section unless generating datasets for other sites) ####

#select site of interest
focus = 'NHC'

#read in full streampulse dataset (not including NEON or Powell)
sp = read.csv('~/Downloads/all_sp_data.csv', stringsAsFactors=FALSE)

#filter all but NC sites
sp = sp[sp$regionID == 'NC',]

#get rid of regionID since this is a huge dataset and it's just wasting memory
sp$regionID = NULL

#grab the NHC subset
subset = sp[sp$siteID == focus,]

#remove the main dataset from the global environment if youre low on RAM
rm(sp)

#remove duplicate rows
subset = subset[! duplicated(subset),]

#convert datetime from string to POSIXct
subset$dateTimeUTC = as.POSIXct(subset$dateTimeUTC, tz='UTC')

#remove bad data
subset = subset[subset$flagID != 'Bad Data',]

#get indices of all points that have been flagged with the word "storm".
#i just flagged all likely storms for Eno 2017 and NHC 2017
storm_related_rows = grep('storm', subset$flagComment, ignore.case=TRUE)

#filter all but storm-related rows (jk)
# subset = subset[storm_related_rows,]

#get set of all storm related flag comments
unique(subset$flagComment[storm_related_rows])

#get datetimes associated with those i just flagged (comment == 'storm')
#(you could include others, but look at them with the cleaning tool first)
storm_rows = grepl('storm', subset$flagComment)
storm_dt = subset$dateTimeUTC[storm_rows]

#filter all but storm rows (jk)
# subset = subset[subset$dateTimeUTC %in% storm_dt,]

#average rows with the same siteID, datetime, and variable (for shape format step below)
subset = aggregate(value ~ siteID + dateTimeUTC + variable, mean,
    data=subset, na.action=NULL)

# create boolean column for storm or not storm
subset$storm = as.numeric(subset$dateTimeUTC %in% storm_dt)

# convert from long to wide format
stormdata = tidyr::spread(subset, variable, value)

#bring in metabolism estimates; deal with \\N characters from MySQL
modres = read.csv('~/Downloads/all_daily_model_results.csv',
    stringsAsFactors=FALSE)
modres[modres == '\\N'] = NA
numind = grepl('(GPP|ER|K600)', colnames(modres))
modres[,numind] = apply(modres[,numind], 2, as.numeric)

#convert datetime to date (make new column)
modres$date = as.Date(as.POSIXct(modres$solar_date))

#clean
modres = filter(modres, site == focus & year != 1907) %>%
    select(-region, -site, -year, -solar_date)

#average storm data by day and merge with metab estimates
stormdaily = group_by(stormdata, 'date'=as.Date(substr(dateTimeUTC, 1, 10))) %>%
    summarize_if(is.numeric, mean, na.rm=TRUE) %>%
    left_join(modres, by='date') %>%
    as.data.frame()

#save csvs
write.csv(stormdata, '~/Desktop/NHC_stormdata.csv', row.names=FALSE)
write.csv(stormdaily, '~/Desktop/NHC_stormdaily.csv', row.names=FALSE)

# 2. some basic exploration ####
#read in daily storm data
stormdaily = read.csv('~/Downloads/NHC_stormdaily.csv', stringsAsFactors=FALSE,
    colClasses=c('date'='Date'))

#see what you're dealing with
head(stormdaily)

#plot storm occurance relative to metab series
ylims = range(stormdaily$GPP, stormdaily$ER, na.rm=TRUE)
plot(stormdaily$date, stormdaily$GPP, type='l', col='red', ylim=ylims,
    xlab='', ylab='g O2 m^-2 d^-1', main='NHC storms, 2017',
    xlim=c(as.Date('2017-01-01'), as.Date('2017-12-31')))
lines(stormdaily$date, stormdaily$ER, col='blue')
abline(v=stormdaily$date[stormdaily$storm > 0], col='gray', lty=2)

#get storm indices
stormbool = as.logical(stormdaily$storm)

stormrle = rle(stormbool)
stormends = cumsum(stormrle$lengths)[c(FALSE, TRUE)]
stormstarts = stormends - stormrle$lengths[c(FALSE, TRUE)] + 1

stormchunks = mapply(`:`, stormstarts, stormends)
stormfac = c()
for(i in 1:length(stormchunks)){
    chunk = rep(i, times=length(stormchunks[[i]]))
    stormfac = append(stormfac, chunk)
}

#get metab averages during storm days
duringGPP = tapply(stormdaily$GPP[stormbool], stormfac,
    mean, na.rm=TRUE)
duringER = tapply(stormdaily$ER[stormbool], stormfac,
    mean, na.rm=TRUE) * -1

#get pre- and post-storm index collections
during = which(stormbool)
day1pre = stormstarts - 1
day2pre = stormstarts - 2
day3pre = stormstarts - 3
day1post = stormends + 1
day2post = stormends + 2
day3post = stormends + 3

#get storm averages and pre - during, post - during differences
#pre - post would also be good to know
storm_summary = matrix(NA, 10, 7,
    dimnames=list(c('GPPmn', 'ERmn', 'GPP:ERmn', 'GPPse', 'ERse',
        'difGPPmn', 'difERmn', 'difGPP:ERmn', 'difGPPse', 'difERse'),
    c('day3pre', 'day2pre', 'day1pre', 'during', 'day1post', 'day2post', 'day3post')))
duringratio = mean(duringGPP, na.rm=TRUE) / mean(duringER, na.rm=TRUE)
for(i in 1:ncol(storm_summary)){
    day = colnames(storm_summary)[i]
    gppset = stormdaily$GPP[get(day)]
    erset = stormdaily$ER[get(day)] * -1

    storm_summary[1, i] = mean(gppset, na.rm=TRUE)
    storm_summary[2, i] = mean(erset, na.rm=TRUE)
    storm_summary[3, i] = storm_summary[1, i] / storm_summary[2, i]
    storm_summary[4, i] = sd(gppset, na.rm=TRUE) / sqrt(length(stormstarts))
    storm_summary[5, i] = sd(erset, na.rm=TRUE) / sqrt(length(stormstarts))

    if(day != 'during'){
        # gppdif = gppset - duringGPP
        # erdif = erset - duringER

        # storm_summary[6, i] = mean(gppdif, na.rm=TRUE)
        # storm_summary[7, i] = mean(erdif, na.rm=TRUE)
        storm_summary[6, i] = storm_summary[1,i] - mean(duringGPP, na.rm=TRUE)
        storm_summary[7, i] = storm_summary[2,i] - mean(duringER, na.rm=TRUE)
        storm_summary[8, i] = storm_summary[3, i] - duringratio
        storm_summary[9, i] = storm_summary[1,i] - sd(duringGPP, na.rm=TRUE) #wrong
        storm_summary[10, i] = storm_summary[2,i] - sd(duringER, na.rm=TRUE) #wrong
        # storm_summary[9, i] = sd(gppdif, na.rm=TRUE)
        # storm_summary[10, i] = sd(erdif, na.rm=TRUE)
    } else {
        storm_summary[6:10, i] = NA
    }
}

#plot
error.bars <- function(x, y, upper, lower=upper, cap.length=0.1, horiz=F,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("One or more vectors is not the same length")

    if(horiz==F) {
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=cap.length, ...)
    } else if (horiz==T) {
        arrows(x+upper,y, x-lower, y, angle=90, code=3, length=cap.length, ...)
    }
}

png(width=14, height=10, units='in', res=300,
    filename='~/Desktop/NHC_stormMetabDif.png', type='cairo')

par(mfrow=c(2, 1), mar=c(3, 5, 3, 1))
gppymax = colSums(storm_summary[c(1,4),])
erymax = colSums(storm_summary[c(2,5),])
barplot(storm_summary[1:3,], beside=TRUE, names.arg=rep('', 7), main='Day mean',
    ylim=c(0, max(c(gppymax, erymax))), col=c('red3', 'blue2', 'purple3'),
    border=c('red3', 'blue2', 'purple3'), ylab='g O2 m^-2 d^-1',
    cex.names=1.5, cex.lab=1.5, cex.main=1.5)
legend('topright', legend=c('GPP', 'ER', 'GPP:ER'),
    fill=c('red3', 'blue2', 'purple3'), bty='n')
barx = seq(1.5, 26.5, 4)
error.bars(barx, storm_summary[1,], upper=storm_summary[4,],
    lower=storm_summary[4,], cap.length=0.01)
error.bars(barx + 1, storm_summary[2,], upper=storm_summary[5,],
    lower=storm_summary[5,], cap.length=0.01)
barplot(storm_summary[6:8,], beside=TRUE, main='Day mean minus storm mean',
    col=c('red3', 'blue2', 'purple3'), ylab='g O2 m^-2 d^-1',
    border=c('red3', 'blue2', 'purple3'), cex.names=1.5, cex.lab=1.5, cex.main=1.5)
dev.off()
