
determine_sample_interval = function(d){

    intervals = int_diff(d$DateTime_UTC) %>%
        time_length(unit='minute') %>%
        table()
    primary_interval_m = as.numeric(names(intervals)[which.max(intervals)])

    return(primary_interval_m)
}

populate_missing_rows = function(d, samp_int){
    #samp_int is the sample interval in minutes

    daterange = range(d$DateTime_UTC, na.rm=TRUE)
    dt_filled = data.frame(DateTime_UTC=seq(daterange[1], daterange[2],
                                            by=paste(samp_int, 'min')))
    d = right_join(d, dt_filled, by='DateTime_UTC') %>%
        arrange(DateTime_UTC)

    if('upload_id' %in% colnames(d)){
        d = mutate(d, upload_id=imputeTS::na_locf(upload_id))
    }

    return(d)
}

nearest_k_days = function(dwg, dwog, cmat, k, minvars){
    # dl: a list of dataframes, named by datestring, each containing all samples
    #   for the specified day
    # k: number of similar days to find for each day with missing data
    # minvars: min number of variables that can be used to infer a
    #   missing value

    if(length(dwg) == 0){
        usr_msg_code <<- '1'
        return('nothing to do')
    }

    if(length(dwog) < 30){
        usr_msg_code <<- '2'
        return('nothing to do')
    }


    # dwg_dtcols = lapply(dwg, function(x) x[, 'DateTime_UTC', drop=FALSE])
    # dwog_dtcols = lapply(dwog, function(x) x[, 'DateTime_UTC', drop=FALSE])
    # nondt_cols = which(colnames(dwg[[1]]) != 'DateTime_UTC')
    # dwg = lapply(dwg, function(x) x[, nondt_cols])

    nearest_days = matrix(NA, length(dwg), k, dimnames=list(names(dwg), NULL))
    for(i in 1:length(dwg)){

        day_scores = vector('numeric', length(dwog))
        for(j in 1:length(dwog)){

            #RMSD: root mean squared difference
            #CSRMSD: correlation-scaled RMSD
            #MCSRMSD: mean CSRMSD

            #get MCSRSD between each day and the current dwg
            squarediffs = try( (dwg[[i]] - dwog[[j]])^2 )
            # (dwg[[i]][, nondt_cols] - dwog[[j]][, nondt_cols])^2
            if(class(squarediffs) == 'try-error'){
                day_scores[j] = Inf
            } else {
                RMSDs = sqrt(colMeans(squarediffs, na.rm=FALSE))
                gapvars = names(RMSDs[which(is.na(RMSDs))])
                inverse_mean_gapvar_abscorrs = 1 -
                    colMeans(abs(cmat[gapvars, , drop=FALSE]))
                CSRMSDs = RMSDs * inverse_mean_gapvar_abscorrs

                if(sum(! is.na(CSRMSDs)) < minvars){
                    day_scores[j] = NA
                } else {
                    MCSRMSD = mean(CSRMSDs, na.rm=TRUE)
                    day_scores[j] = MCSRMSD
                }
            }
        }

        nearest_days[i,] = names(dwog)[order(day_scores)[1:3]]
    }

    return(nearest_days)
}

standardize_df = function(x){

    x_std = x
    datecols = which(colnames(x) %in% c('date', 'DateTime_UTC'))
    x_means = colMeans(x[, -datecols], na.rm=TRUE)
    x_sds = apply(x[, -datecols], 2, sd, na.rm=TRUE)
    x_std[, -datecols] = scale(x[, -datecols])

    return(list(x_std=x_std, x_means=x_means, x_sds=x_sds))
}

unstandardize_df = function(x_std, x_means, x_sds){

    x = x_std
    datecol = which(colnames(x) %in% c('date', 'DateTime_UTC'))
    x[, -datecol] = mapply(function(w, y, z) w * y + z,
                           x_std[, -datecol], x_sds, x_means)

    return(x)
}

listify = function(dfx, incomplete_days){

    df_incompletes = dfx %>%
        filter(date %in% incomplete_days) %>%
        select(-date)

    dtcol_fulldays = dfx %>%
        filter(! date %in% incomplete_days) %>%
        select(DateTime_UTC)

    datelist_fulldays = dfx %>%
        select(-DateTime_UTC) %>%
        filter(! date %in% incomplete_days) %>%
        plyr::dlply(plyr::.(date))

    datelist_fulldays = Map(function(x) as.matrix(select(x, -date)),
                            datelist_fulldays)

    return(list(fulldays=datelist_fulldays, dtcol=dtcol_fulldays,
                incompletes=df_incompletes))
}

unlistify = function(dwg, dwog, orig_dtcol, incompletes){

    datelist_fulldays = append(dwg, dwog)
    listdates = as.Date(names(datelist_fulldays))
    datelist_fulldays = datelist_fulldays[order(listdates)]

    dfx = plyr::ldply(datelist_fulldays) %>%
        select(-.id) %>%
        bind_cols(orig_dtcol) %>%
        bind_rows(incompletes) %>%
        arrange(DateTime_UTC) %>%
        select(DateTime_UTC, everything())

    return(dfx)
}

gap_fill = function(dwg, dwog, nearest_days){

    for(i in 1:length(dwg)){
        ddate = names(dwg[i])
        # lapply(d, function(x) rle(is.na(x))) #start here if maxgap wanted
        # too_long_na_series =
        dnearest = nearest_days[ddate,]
        dnmean = Reduce(`+`, dwog[dnearest]) / length(dnearest)
        dyhat = dwg[[i]] - dnmean
        dyhat = apply(dyhat, 2,
                      function(x){
                          if(sum(! is.na(x)) < 4) return(x)
                          imputeTS::na_interpolation(x, option='linear')
                      })
        dwg[[i]] = dnmean + dyhat
    }

    return(dwg)
}

find_snappoints = function(x, ndi_sections_, original_indices, original_dates,
                           unpaired_ndi_tracker){

    ds = median(x)

    if(length(x) %% 2 == 0){
        ds = c(floor(ds), ceiling(ds))
        ndi_datetimes = ndi_sections_[ds, 'DateTime_UTC']
        ndi_difftimes = difftime(ndi_datetimes,
                                 as.Date(ndi_datetimes[2]), units='days')
        ds = ds[which.min(abs(ndi_difftimes))]
    }

    # plot(ndiout$DateTime_UTC, ndiout[[c]], col='orange', type='l', lwd=2)
    # lines(pldf$DateTime_UTC, pldf[[c]], lwd=2)
    gapsize = as.numeric(difftime(original_dates[ds],
                                  original_dates[ds - 1], units='days'))
    # dend = Position(function(x) ! is.na(x),
    #     ndiout_[lb:(ds - 1), c], right=TRUE) + lb - 1

    fullday_skip = unpaired = FALSE
    if(gapsize > 1) fullday_skip = TRUE

    ds = original_indices[ds]

    if(length(unpaired_ndi_tracker)){
        unpaired = any(sapply(unpaired_ndi_tracker,
                              function(x) ds >= x$start && ds <= x$stop))
        if(unpaired){
            fullday_skip = TRUE
            unpaired = TRUE
        }
    }

    # nsamps_to_midday = 24 * 60 / interv / 2
    # dbs = c(max(ds - nsamps_to_midday, original_indices[1]),
    #     #never going to reach midday, right?
    #     ds + nsamps_to_midday)
    # dbs = c()

    return(list('daystart'=ds, 'fullday_skip'=fullday_skip, 'unpaired'=unpaired))
}

NDI = function(df, knn=3, interv){

    samples_per_day = 24 * 60 / interv

    df_ = df %>%
        mutate(date=as.Date(lubridate::floor_date(DateTime_UTC, unit='day')))

    # date_counts = df_ %>%
    #     select(date, !!sym(target_variable)) %>%
    #     filter(! is.na(!!sym(target_variable))) %>%
    #     pull(date) %>%
    #     table()

    date_counts = table(df_$date)
    incomplete_days_ind = which(date_counts < samples_per_day)
    if(length(incomplete_days_ind) == 0){
        incomplete_days = NULL
    } else {
        incomplete_days = as.Date(names(incomplete_days_ind))
    }

    #make new df of standardized variables; grab moments by which to rescale
    std_out = standardize_df(df_)
    df_std = std_out$x_std

    listout = listify(df_std, incomplete_days)
    datelist_fulldays = listout$fulldays

    # find k nearest neighbors for each day with missing values
    days_with_gaps_bool = unlist(Map(function(x) any(is.na(x)), datelist_fulldays))
    days_with_gaps = datelist_fulldays[days_with_gaps_bool]
    days_without_gaps = datelist_fulldays[! days_with_gaps_bool]
    correlations = cor(select(df_std, -date, -DateTime_UTC),
                       use='na.or.complete', method='spearman')

    nearest_days = nearest_k_days(days_with_gaps, days_without_gaps,
                                  cmat=correlations, k=knn, minvars=3)

    if(length(nearest_days) == 1 && nearest_days == 'nothing to do'){
        return(df)
    }

    #perform NDI; return data to original shape and format
    gapdays_filled = gap_fill(days_with_gaps, days_without_gaps, nearest_days)

    df_std = unlistify(gapdays_filled, days_without_gaps, listout$dtcol,
                       listout$incompletes)
    df_ = unstandardize_df(df_std, std_out$x_means, std_out$x_sds)
    # df = as.data.frame(select(df_, -DateTime_UTC))
    df = as.data.frame(df_)

    return(df)
}

rle_custom = function(x){
    r = rle(x)
    ends = cumsum(r$lengths)
    r = cbind(values=r$values,
              starts=c(1, ends[-length(ends)] + 1),
              stops=ends, lengths=r$lengths, deparse.level=1)
    return(r)
}

snap_days = function(dfcols_, flagdf_, ndiout_, pldf_, interv){

    for(c in dfcols_[dfcols_ != 'DateTime_UTC']){
        print(c)

        twos = which(flagdf_[[c]] == 2) #code 2 marks NDI points
        # before_first_two = twos[1] - 1
        # if(before_first_two >= 1){
        #     actual_start_of_twos = Position(function(x) ! is.na(x),
        #         ndiout_[1:before_first_two, c], right=TRUE)
        #     if(actual_start_of_twos != before_first_two){
        #         twos = c(actual_start_of_twos: #ugh: doesn't just apply to the first index of twos
        # }
        # plot(ndiout_$DateTime_UTC, ndiout_[[c]], col='orange', type='l', lwd=2)
        # lines(pldf_$DateTime_UTC, pldf_[[c]], lwd=2)

        #some NDI sections won't be bounded by imputed points; extend them
        #to include adjacent NA sections
        na_runs = rle_custom(is.na(ndiout_[[c]]))
        na_runs = na_runs[na_runs[, 'values'] == 1, , drop=FALSE]
        # two_runs = rle_custom(diff(twos))
        # if(nrow(two_runs)){
        #     two_runs[, 2] = two_runs[, 2]
        two_runs = rle_custom(flagdf_[[c]] == 2)
        two_runs = two_runs[two_runs[, 'values'] == 1, , drop=FALSE]

        if(! nrow(two_runs)) next

        unpaired_ndi_tracker = list()
        for(i in 1:nrow(two_runs)){

            current_run = two_runs[i, ]

            abutting_na_run_bool = mapply(
                function(x, y, xx){
                    xx[2] == x - 1 || xx[1] == y + 1
                },
                na_runs[, 2], na_runs[, 3], MoreArgs=list('xx'=current_run[2:3]))

            if(any(abutting_na_run_bool)){

                abutting_na_run = na_runs[abutting_na_run_bool,]
                previous_value = flagdf_[abutting_na_run[2] - 1, c]
                ndi_unpaired_l = previous_value != 2 || is.na(previous_value)
                next_value = flagdf_[abutting_na_run[3] + 1, c]
                ndi_unpaired_r = next_value != 2 || is.na(next_value)

                ll = length(unpaired_ndi_tracker)
                if(ndi_unpaired_l){
                    current_run[2] = abutting_na_run[2]
                    unpaired_ndi_tracker[[ll + 1]] = list(side='l',
                                                          start=unname(current_run[2]), stop=unname(current_run[3]))
                } else if(ndi_unpaired_r){
                    current_run[3] = abutting_na_run[3]
                    unpaired_ndi_tracker[[ll + 1]] = list(side='r',
                                                          start=unname(current_run[2]), stop=unname(current_run[3]))
                }

                pre_splice = twos[twos < current_run['starts']]
                post_splice = twos[twos > current_run['stops']]
                twos = c(pre_splice, current_run['starts']:current_run['stops'],
                         post_splice)
            }
        }

        #carry on now that NDI sections include adjacent missing data
        if(! length(twos)) next

        twodates = ndiout_$DateTime_UTC[twos]
        ndi_sections = ndiout_[twos, c('DateTime_UTC', c)]

        #get new day indices within NDI sections
        ndi_rounded_dates = round(ndi_sections$DateTime_UTC, 'hour')
        daystarts = which(substr(ndi_rounded_dates, 12, 19) ==
                              '00:00:00')
        # twodates[daystarts]

        if(length(daystarts)){

            #this block FAILS if sample interval >= 30m; skipping this block
            #may result in failure
            if(interv < 30){

                #remove daystarts that don't correspond to NDI imputed values
                druns = rle_custom(diff(daystarts))
                druns = druns[druns[, 'lengths'] > 1, , drop=FALSE] #FAILS HERE
                druns[, 'stops'] = druns[, 'stops'] + 1 #account for diffing

                if(nrow(druns)){
                    forreals_daystart_inds = unlist(mapply(function(x, y) x:y,
                                                           druns[, 'starts'], druns[, 'stops'], SIMPLIFY=FALSE))

                    # ! daystarts %in% daystarts[forreals_daystart_inds]
                    daystarts = daystarts[forreals_daystart_inds]

                } else {
                    next
                }
            }

        } else {
            next
        }

        r = rle_custom(diff(daystarts))
        r = r[r[, 'values'] == 1, , drop=FALSE]

        if(nrow(r) > 1){
            r = r[-nrow(r), , drop=FALSE]
        }

        ds_fac = cut(1:length(daystarts), c(0, r[, 'stops'] + 1, Inf))
        ndi_ind_list = split(daystarts, ds_fac, drop=TRUE)
        #some mid-gap daystarts may have been included; remove them.
        ndi_ind_list = Filter(function(x) any(! is.na(ndiout_[twos[x], c])),
                              ndi_ind_list)
        ndi_list = unname(lapply(ndi_ind_list,
                                 find_snappoints, ndi_sections, twos, twodates, unpaired_ndi_tracker))
        # ndi_list = Filter(function(x) ! is.null(x), ndi_list)

        if(! length(ndi_list)){
            next
        }

        #identify NDI sections in which multiple days abut.
        #these need to be edge-snapped
        ndi_edges = rle_custom(flagdf_[[c]] == 2)
        ndi_edges = ndi_edges[ndi_edges[, 'values'] == 1, , drop=FALSE]
        na_edges = rle_custom(is.na(pldf_[[c]]))
        na_edges = na_edges[na_edges[, 'values'] == 1, , drop=FALSE]

        fullday_skips = matrix(ncol=4, nrow=2,
                               dimnames=list(NULL, c('values', 'starts', 'stops', 'lengths')))
        fullday_skips = unpaired_ndi = fullday_skips[-(1:2),] #create empty matrices

        for(j in nrow(ndi_edges):1){
            gapless_ndi = paste(ndi_edges[j, 'starts'],
                                ndi_edges[j, 'stops']) %in%
                paste(na_edges[, 'starts'], na_edges[, 'stops'])
            if(! gapless_ndi){
                fullday_skips = rbind(fullday_skips,
                                      ndi_edges[j, , drop=FALSE])
                ndi_edges = ndi_edges[-j, , drop=FALSE]
            }
        }
        if(nrow(fullday_skips)){
            fullday_skips = fullday_skips[nrow(fullday_skips):1, , drop=FALSE]
        }

        #sequester fullday skips that are also unpaired NDI sections
        if(length(unpaired_ndi_tracker)){

            unpaired_starts = sapply(unpaired_ndi_tracker,
                                     function(x) x$start)
            unpaired_stops = sapply(unpaired_ndi_tracker,
                                    function(x) x$stop)

            for(j in nrow(fullday_skips):1){

                fds = fullday_skips[j, ]
                if(fds[2] %in% unpaired_starts || fds[3] %in% unpaired_stops){
                    unpaired_ndi = rbind(unpaired_ndi,
                                         fullday_skips[j, , drop=FALSE])
                    fullday_skips = fullday_skips[-j, , drop=FALSE]
                }
            }
            unpaired_ndi = unpaired_ndi[nrow(unpaired_ndi):1, , drop=FALSE]
        }

        #ignore within-day NDI sections; these need no action
        if(nrow(ndi_edges)){

            all_daystarts = sapply(ndi_list, function(x) x$daystart)
            for(j in nrow(ndi_edges):1){

                sect = ndi_edges[j, ]
                contains_daystart = rep(NA, length(all_daystarts))
                for(k in 1:length(all_daystarts)){
                    contains_daystart[k] = sect[2] <= all_daystarts[k] &&
                        sect[3] >= all_daystarts[k]
                }

                if(! any(contains_daystart)){
                    ndi_edges = ndi_edges[-j, , drop=FALSE]
                }
            }
        }

        #make sure everything is chill.
        #associate NDI section boundaries with newday indices
        if(! nrow(fullday_skips) && ! length(unpaired_ndi_tracker) &&
           ! nrow(ndi_edges)) next

        if(nrow(fullday_skips)){

            fds_bool = sapply(ndi_list,
                              function(x) x$fullday_skip && ! x$unpaired)

            if(any(fds_bool)){
                for(k in 1:length(ndi_list[fds_bool])){
                    ds_k = ndi_list[fds_bool][[k]]$daystart
                    if(ds_k != fullday_skips[k * 2, 2]){
                        stop('not chill')
                    }
                }
            }
            # fullday_skips[c(1, 2, 9, 10, 13, 14, 19, 20),]

            #identify bounds for paired fullday skips
            for(k in 1:length(ndi_list[fds_bool])){
                fds_bounds = c(fullday_skips[k * 2 - 1, 2],
                               fullday_skips[k * 2, 3])
                ndi_list[fds_bool][[k]]$bounds = fds_bounds
            }

        }

        if(length(unpaired_ndi_tracker)){

            unpaired_ndi = alply(unpaired_ndi, 1, function(x){
                list(matrix(x, nrow=1))
            })

            unpaired_bool = sapply(ndi_list,
                                   function(x) x$fullday_skip && x$unpaired)

            #identify bounds for unpaired fullday skips
            ndi_list[unpaired_bool] = mapply(function(x, y){
                x$bounds = c(y$start, y$stop)
                return(x)
            }, ndi_list[unpaired_bool], unpaired_ndi_tracker, SIMPLIFY=FALSE)

        }

        if(nrow(ndi_edges)){

            ndi_edges = alply(ndi_edges, 1, function(x){
                list(matrix(x, nrow=1))
            })

            normday_bool = sapply(ndi_list,
                                  function(x) ! x$fullday_skip && ! x$unpaired)

            if(any(normday_bool)){
                ds_db_lined_up = all(mapply(function(x, y){
                    x$daystart >= y[[1]][2] && x$daystart <= y[[1]][3]
                }, ndi_list[normday_bool], ndi_edges))

                if(! ds_db_lined_up){
                    stop('not chill')
                }
            }

            #identify bounds for normal days
            ndi_list[normday_bool] = mapply(function(x, y){
                x$bounds = y[[1]][2:3]
                return(x)
            }, ndi_list[normday_bool], ndi_edges, SIMPLIFY=FALSE)

        }


        #snap NDI day edges:
        #adjust values linearly toward midpoint in proportion to the length
        #of each dangling side
        unpaired_counter = 0
        for(j in 1:length(ndi_list)){

            ndi_section = ndi_list[[j]]

            if(ndi_section$unpaired){
                unpaired_counter = unpaired_counter + 1
                unp_side = unpaired_ndi_tracker[[unpaired_counter]]$side
            }

            ds = ndi_section$daystart
            lb = ndi_section$bounds[1]
            rb = ndi_section$bounds[2]
            dend = Position(function(x) ! is.na(x),
                            ndiout_[lb:(ds - 1), c], right=TRUE) + lb - 1
            # points(ndiout_$DateTime_UTC[dend],
            #     ndiout_$AirPres_kPa[dend], col='blue')
            # points(ndiout_$DateTime_UTC[c(ds, dend)],
            #     ndiout_$AirPres_kPa[c(ds, dend)], col='blue')
            # points(ndiout_$DateTime_UTC[c(lb, rb)],
            #     ndiout_$AirPres_kPa[c(lb, rb)], col='darkgreen')

            ndilen = rb - lb + 1
            na_len = sum(is.na(ndiout_[lb:rb, c]))
            if(! ndi_section$unpaired){
                snapdist = diff(ndiout_[c(dend, ds), c])
            } else {
                if(unp_side == 'l'){
                    snapdist = diff(ndiout_[c(lb - 1, ds), c])
                } else {
                    snapdist = diff(ndiout_[c(dend, rb + 1), c])
                }
                if(is.na(snapdist)) next #unpaired side is the end of the series
            }

            snapdist = snapdist * ((ndilen - na_len) / ndilen)
            ndilen_l = dend - lb + 1
            ndilen_r = rb - ds + 1
            snapprop_l = ndilen_l / ndilen
            snapprop_r = ndilen_r / ndilen

            if(! is.na(snapprop_l)){
                snapshift_l = seq(0, snapdist * snapprop_l,
                                  length.out=ndilen_l)
                ndiout_[lb:dend, c] = ndiout_[lb:dend, c] + snapshift_l
            }

            if(! is.na(snapprop_r)){
                snapshift_r = seq(-1 * snapdist * snapprop_r, 0,
                                  length.out=ndilen_r)
                ndiout_[ds:rb, c] = ndiout_[ds:rb, c] + snapshift_r
            }

        }
    }

    return(ndiout_)
}
