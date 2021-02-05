library(geoknife)
library(prism)
library(tidyverse)
library(lubridate)
# library(plyr)
# library(sf)
# library(MODISTools)
library(nhdplusTools)
library(sf)
library(streamstats)
library(glue)
library(rgee)

#only for section 10
# rgee::ee_Initialize(email = 'email here',
#                     drive = TRUE)

rebuild_all = FALSE

# 1. setup and helper functions ####

setwd('~/git/papers/alice_nhc/')

sites = tibble(region = 'NC', sitecode = 'NHC', sitename = 'New Hope Creek',
               latitude = 35.9795, longitude = -79.0018)
EPSG = 4326 #EPSG code

#this should be updated with a "specs" option
delineate_watershed_from_point <- function(lat,
                                           long,
                                           crs,
                                           machine_status = 'n00b',
                                           write_dir,
                                           write_name,
                                           verbose = TRUE){

    #lat: numeric representing latitude of the pour point in decimal degrees
    #   (negative indicates southern hemisphere)
    #long: numeric representing longitude of the pour point in decimal degrees
    #   (negative indicates west of prime meridian)
    #crs: numeric representing the coordinate reference system (e.g. 4326 for WSG84)
    #machine_status: either '1337', indicating that your machine has >= 16 GB
    #   RAM, or 'n00b', indicating < 16 GB RAM. DEM resolution is chosen accordingly
    #write_dir: character. the directory in which to write output shapefile
    #write_name: character. the basename of the shapefile components to be written. i.e.
    #   <write_name>.shp, <write_name>.shx, <write_name>.prj, <write_name>.dbf
    #verbose: logical. determines the amount of informative messaging during run

    #details: Output will have CRS 4326 (WGS 84), though processing is done
    #   on projected data. the projection specifications are determined
    #   automatically, based on pour point location. Note that for mega-huge
    #   watersheds, this could result in an inaccurate watershed area calculation.

    #returns: a list containing the following components:
    #   watershed_area_ha: the area of the delineated watershed in hectares
    #       (meters squared divided by 10,000)
    #   buffer_radius_m: the width (meters) around the site location that was used when
    #       requesting a DEM (digital elevation model)
    #   snap_distance_m: the search radius (meters) around the pour point that was used
    #       to choose a stream to snap the pour point to.
    #   snap_method: either "standard", which snaps the pour point to the cell
    #       within snap_distance_m with highest flow accumulation, or "jenson",
    #       which snaps to the nearest flow line
    #   dem_resolution: passed to elevatr::get_elev_raster (z parameter).
    #       depends on supplied machine_status

    #NOTE: in addition to the packages loaded below, you'll need mapview to
    #   visualize delineated watershed boundaries

    require(dplyr)
    require(glue)
    require(sf)
    require(elevatr)
    require(raster)
    require(whitebox)

    sm <- suppressMessages
    sw <- suppressWarnings

    #the functions below are all helpers for subroutines of
    #   delineate_watershed_from_point. they should all stand alone pretty
    #   well, and some are generally useful. they are defined locally here,
    #   just for convenient distribution.

    #moving shapefiles can be annoying, since they're actually represented by
    #   3-4 files
    move_shapefiles <- function(shp_files, from_dir, to_dir,
                                new_name_vec = NULL){

        #shp_files is a character vector of filenames with .shp extension
        #   (.shx, .prj, .dbf are handled internally and don't need to be listed)
        #from_dir and to_dir are strings representing the source and destination
        #   directories, respectively
        #new_name_vec is an optional character vector of new names for each shape file.
        #   these can end in ".shp", but don't need to

        if(any(! grepl('\\.shp$', shp_files))){
            stop('All components of shp_files must end in ".shp"')
        }

        if(length(shp_files) != length(new_name_vec)){
            stop('new_name_vec must have the same length as shp_files')
        }

        dir.create(to_dir,
                   showWarnings = FALSE,
                   recursive = TRUE)

        for(i in 1:length(shp_files)){

            shapefile_base <- strsplit(shp_files[i], '\\.shp')[[1]]

            files_to_move <- list.files(path = from_dir,
                                        pattern = shapefile_base)

            extensions <- stringr::str_match(files_to_move,
                                             paste0(shapefile_base, '(\\.[a-z]{3})'))[, 2]

            if(is.null(new_name_vec)){
                new_name_base <- rep(shapefile_base, length(files_to_move))
            } else {
                new_name_base <- strsplit(new_name_vec[i], '\\.shp$')[[1]]
                new_name_base <- rep(new_name_base, length(files_to_move))
            }

            tryCatch({

                #try to move the files (may fail if they are on different partitions)
                mapply(function(x, nm, ext) file.rename(from = paste(from_dir,
                                                                     x,
                                                                     sep = '/'),
                                                        to = glue('{td}/{n}{ex}',
                                                                  td = to_dir,
                                                                  n = nm,
                                                                  ex = ext)),
                       x = files_to_move,
                       nm = new_name_base,
                       ext = extensions)

            }, warning = function(w){

                #if that fails, copy them and then delete them
                mapply(function(x, nm, ext) file.copy(from = paste(from_dir,
                                                                   x,
                                                                   sep = '/'),
                                                      to = glue('{td}/{n}{ex}',
                                                                td = to_dir,
                                                                n = nm,
                                                                ex = ext),
                                                      overwrite = TRUE),
                       x = files_to_move,
                       nm = new_name_base,
                       ext = extensions)

                lapply(paste(from_dir,
                             files_to_move,
                             sep = '/'),
                       unlink)
            })
        }

        return()
    }

    #prompt users for stuff, provide single-character responses,
    #   reprompt if they don't choose one of the expected responses
    get_response_1char <- function(msg, possible_chars,
                                   subsequent_prompt = FALSE){

        #msg: character. a message that will be used to prompt the user
        #possible_chars: character vector of acceptable single-character responses

        if(subsequent_prompt){
            cat(paste('Please choose one of:',
                      paste(possible_chars,
                            collapse = ', '),
                      '\n> '))
        } else {
            cat(msg)
        }

        ch <- as.character(readLines(con = stdin(), 1))

        if(length(ch) == 1 && ch %in% possible_chars){
            return(ch)
        } else {
            get_response_1char(msg, possible_chars, subsequent_prompt = TRUE)
        }
    }

    #chooose an appropriate projection, based on location
    choose_projection <- function(lat = NULL, long = NULL,
                                  unprojected = FALSE){

        if(unprojected){
            PROJ4 <- glue('+proj=longlat +datum=WGS84 +no_defs ',
                          '+ellps=WGS84 +towgs84=0,0,0')
            return(PROJ4)
        }

        if(is.null(lat) || is.null(long)){
            stop('If projecting, lat and long are required.')
        }

        if(lat <= 15 && lat >= -15){ #equatorial
        '+proj=laea +lat_0=35.9795 +lon_0=-79.0018'
            PROJ4 = glue('+proj=laea +lon_0=', long,
                         ' +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
        } else { #temperate or polar
            PROJ4 = glue('+proj=laea +lat_0=', lat, ' +lon_0=', long,
                         ' +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
        }

        return(PROJ4)
    }

    #for determining whether the DEM extent wasn't big enough to allow full
    #   delineation
    raster_intersection_summary <- function(wb, dem){

        #wb is a delineated watershed boundary as a rasterLayer
        #dem is a DEM rasterLayer

        summary_out <- list()

        #convert wb to sf object
        wb <- wb %>%
            raster::rasterToPolygons() %>%
            sf::st_as_sf()

        #get edge of DEM as sf object
        dem_edge <- raster::boundaries(dem) %>%
            raster::reclassify(matrix(c(0, NA),
                                      ncol = 2)) %>%
            raster::rasterToPolygons() %>%
            sf::st_as_sf()

        #tally raster cells
        summary_out$n_wb_cells <- length(wb$geometry)
        summary_out$n_dem_cells <- length(dem_edge$geometry)

        #tally intersections; calc percent of wb cells that overlap
        intersections <- sf::st_intersects(wb, dem_edge) %>%
            as.matrix() %>%
            apply(MARGIN = 2,
                  FUN = sum) %>%
            table()

        true_intersections <- sum(intersections[names(intersections) > 0])

        summary_out$n_intersections <- true_intersections
        summary_out$pct_wb_cells_intersect <- true_intersections /
            summary_out$n_wb_cells * 100

        return(summary_out)
    }

    #the workhorse
    delineate_watershed_apriori <- function(lat, long, crs,
                                            machine_status = 'n00b',
                                            verbose = FALSE){

        #lat: numeric representing latitude in decimal degrees
        #   (negative indicates southern hemisphere)
        #long: numeric representing longitude in decimal degrees
        #   (negative indicates west of prime meridian)
        #crs: numeric representing the coordinate reference system (e.g. WSG84)
        #machine_status: either '1337', indicating that your machine has >= 16 GB
        #   RAM, or 'n00b', indicating < 16 GB RAM. DEM resolution is chosen accordingly
        #verbose: logical. determines the amount of informative messaging during run

        #returns the location of candidate watershed boundary files

        tmp <- tempdir()
        inspection_dir <- glue(tmp, '/INSPECT_THESE')
        dem_f <- glue(tmp, '/dem.tif')
        point_f <- glue(tmp, '/point.shp')
        d8_f <- glue(tmp, '/d8_pntr.tif')
        flow_f <- glue(tmp, '/flow.tif')

        dir.create(path = inspection_dir,
                   showWarnings = FALSE)

        proj <- choose_projection(lat = lat,
                                  long = long)

        site <- tibble(x = lat,
                       y = long) %>%
            sf::st_as_sf(coords = c("y", "x"),
                         crs = crs) %>%
            sf::st_transform(proj)
        # sf::st_transform(4326) #WGS 84 (would be nice to do this unprojected)

        #prepare for delineation loops
        buffer_radius <- 100
        dem_coverage_insufficient <- FALSE
        while_loop_begin <- TRUE

        #snap site to flowlines 3 different ways. delineate watershed boundaries (wb)
        #for each unique snap. if the delineations get cut off, get more elevation data
        #and try again
        while(while_loop_begin || dem_coverage_insufficient){

            while_loop_begin <- FALSE

            if(machine_status == '1337'){
                dem_resolution <- case_when(
                    buffer_radius <= 1e4 ~ 12,
                    buffer_radius == 1e5 ~ 11,
                    buffer_radius == 1e6 ~ 10,
                    buffer_radius == 1e7 ~ 8,
                    buffer_radius == 1e8 ~ 6,
                    buffer_radius == 1e9 ~ 4,
                    buffer_radius >= 1e10 ~ 2)
            } else if(machine_status == 'n00b'){
                dem_resolution <- case_when(
                    buffer_radius <= 1e4 ~ 10,
                    buffer_radius == 1e5 ~ 8,
                    buffer_radius == 1e6 ~ 6,
                    buffer_radius == 1e7 ~ 4,
                    buffer_radius == 1e8 ~ 2,
                    buffer_radius >= 1e9 ~ 1)
            } else {
                stop('machine_status must be either "1337" or "n00b"')
            }

            site_buf <- sf::st_buffer(x = site,
                                      dist = buffer_radius)
            dem <- elevatr::get_elev_raster(locations = site_buf,
                                            z = dem_resolution,
                                            verbose = verbose)

            raster::writeRaster(x = dem,
                                filename = dem_f,
                                overwrite = TRUE)

            sf::st_write(obj = site,
                         dsn = point_f,
                         delete_layer = TRUE,
                         quiet = TRUE)

            whitebox::wbt_fill_single_cell_pits(dem = dem_f,
                                                output = dem_f)

            whitebox::wbt_breach_depressions(dem = dem_f,
                                             output = dem_f,
                                             flat_increment = 0.01)

            # whitebox::wbt_burn_streams_at_roads(dem = dem_f,
            #                                     streams = streams_f,
            #                                     roads = roads_f,
            #                                     output = dem_f)
            #                                     # width = 30)

            whitebox::wbt_d8_pointer(dem = dem_f,
                                     output = d8_f)

            whitebox::wbt_d8_flow_accumulation(input = dem_f,
                                               output = flow_f,
                                               out_type = 'catchment area')

            snap1_f <- glue(tmp, '/snap1_jenson_dist150.shp')
            whitebox::wbt_jenson_snap_pour_points(pour_pts = point_f,
                                                  streams = flow_f,
                                                  output = snap1_f,
                                                  snap_dist = 150)
            snap2_f <- glue(tmp, '/snap2_standard_dist50.shp')
            whitebox::wbt_snap_pour_points(pour_pts = point_f,
                                           flow_accum = flow_f,
                                           output = snap2_f,
                                           snap_dist = 50)
            snap3_f <- glue(tmp, '/snap3_standard_dist150.shp')
            whitebox::wbt_snap_pour_points(pour_pts = point_f,
                                           flow_accum = flow_f,
                                           output = snap3_f,
                                           snap_dist = 150)

            #the site has been snapped 3 different ways. identify unique snap locations.
            snap1 <- sf::st_read(snap1_f, quiet = TRUE)
            snap2 <- sf::st_read(snap2_f, quiet = TRUE)
            snap3 <- sf::st_read(snap3_f, quiet = TRUE)
            unique_snaps_f <- snap1_f
            if(! identical(snap1, snap2)) unique_snaps_f <- c(unique_snaps_f, snap2_f)
            if(! identical(snap1, snap3)) unique_snaps_f <- c(unique_snaps_f, snap3_f)

            #good for experimenting with snap specs:
            # delineate_watershed_test2(tmp, point_f, flow_f,
            #                           d8_f, 'standard', 1000)

            #delineate each unique location
            for(i in 1:length(unique_snaps_f)){

                rgx <- stringr::str_match(unique_snaps_f[i],
                                          '.*?_(standard|jenson)_dist([0-9]+)\\.shp$')
                snap_method <- rgx[, 2]
                snap_distance <- rgx[, 3]

                wb_f <- glue('{path}/wb{n}_buffer{b}_{typ}_dist{dst}.tif',
                             path = tmp,
                             n = i,
                             b = buffer_radius,
                             typ = snap_method,
                             dst = snap_distance)

                whitebox::wbt_watershed(d8_pntr = d8_f,
                                        pour_pts = unique_snaps_f[i],
                                        output = wb_f)

                wb <- raster::raster(wb_f)

                #check how many wb cells coincide with the edge of the DEM.
                #If > 0.1% or > 5, broader DEM needed
                smry <- raster_intersection_summary(wb = wb,
                                                    dem = dem)

                if(verbose){
                    print(glue('buffer radius: {br}; snap: {sn}/{tot}; ',
                               'n intersecting cells: {ni}; pct intersect: {pct}',
                               br = buffer_radius,
                               sn = i,
                               tot = length(unique_snaps_f),
                               ni = round(smry$n_intersections, 2),
                               pct = round(smry$pct_wb_cells_intersect, 2)))
                }

                if(smry$pct_wb_cells_intersect > 0.1 || smry$n_intersections > 5){
                    buffer_radius_new <- buffer_radius * 10
                    dem_coverage_insufficient <- TRUE
                } else {
                    buffer_radius_new <- buffer_radius

                    #write and record temp files for the technician to visually inspect
                    wb_sf <- wb %>%
                        raster::rasterToPolygons() %>%
                        sf::st_as_sf() %>%
                        sf::st_buffer(dist = 0.1) %>%
                        sf::st_union() %>%
                        sf::st_as_sf()#again? ugh.

                    wb_sf <- sf::st_transform(wb_sf, 4326) #EPSG for WGS84

                    wb_sf_f <- glue('{path}/wb{n}_BUF{b}{typ}DIST{dst}RES{res}.shp',
                                    path = inspection_dir,
                                    n = i,
                                    b = buffer_radius,
                                    typ = snap_method,
                                    dst = snap_distance,
                                    res = dem_resolution)

                    sf::st_write(obj = wb_sf,
                                 dsn = wb_sf_f,
                                 delete_dsn = TRUE,
                                 quiet = TRUE)

                    # zz = sf::st_read(wb_sf_f)
                    dem_conditioned = sf::st_read(dem_f)
                    mv(dem_conditioned) + mv(roads, color='gray') + mv(streams) + mv(wb_sf)
                }
            }

            buffer_radius <- buffer_radius_new
        } #end while loop

        if(verbose){
            message(glue('Candidate delineations are in: ', inspection_dir))
        }

        return(inspection_dir)
    }

    if(verbose){
        message('Beginning watershed delineation')
    }

    inspection_dir <- sw(delineate_watershed_apriori(
        lat = lat,
        long = long,
        crs = crs,
        machine_status = machine_status,
        verbose = verbose))

    files_to_inspect <- list.files(path = inspection_dir,
                                   pattern = '.shp')

    #if only one delineation, write it to write_dir
    if(length(files_to_inspect) == 1){

        selection <- files_to_inspect[1]

        move_shapefiles(shp_files = selection,
                        from_dir = inspection_dir,
                        to_dir = write_dir)

        if(verbose){
            message(glue('Delineation successful. Shapefile written to ',
                         write_dir))
        }

        #otherwise, inspect all delineations and choose one
    } else {

        nshapes <- length(files_to_inspect)

        wb_selections <- paste(paste0('[',
                                      c(1:nshapes, 'A'),
                                      ']'),
                               c(files_to_inspect, 'Abort delineation'),
                               sep = ': ',
                               collapse = '\n')

        helper_code <- glue('mapview::mapview(sf::st_read("{wd}/{f}"))',
                            wd = inspection_dir,
                            f = files_to_inspect) %>%
            paste(collapse = '\n\n')

        msg <- glue('Visually inspect the watershed boundary candidate shapefiles ',
                    'in {td}, then enter the number corresponding to the ',
                    'one that looks most legit. Here\'s some ',
                    'helper code you can paste into a separate R instance ',
                    ':\n\n{hc}\n\nChoices:\n{sel}\n\nEnter choice here > ',
                    hc = helper_code,
                    sel = wb_selections,
                    td = inspection_dir)

        resp <- get_response_1char(msg = msg,
                                   possible_chars = c(1:nshapes, 'A'))

        if(resp == 'A'){
            message('Delineation aborted')
            return()
        }

        selection <- files_to_inspect[as.numeric(resp)]

        move_shapefiles(shp_files = selection,
                        from_dir = inspection_dir,
                        to_dir = write_dir,
                        new_name_vec = write_name)

        message(glue('Selection {s}:\n\t{sel}\nwas written to:\n\t{sdr}',
                     s = resp,
                     sel = selection,
                     sdr = write_dir))
    }

    #calculate watershed area in hectares
    wb <- sf::st_read(glue('{d}/{s}.shp',
                           d = write_dir,
                           s = write_name),
                      quiet = TRUE)

    ws_area_ha <- as.numeric(sf::st_area(wb)) / 10000

    #return the specifications of the correctly delineated watershed, and some
    #   other goodies
    rgx <- stringr::str_match(selection,
                              paste0('^wb[0-9]+_BUF([0-9]+)(standard|jenson)',
                                     'DIST([0-9]+)RES([0-9]+)\\.shp$'))

    deets <- list(name = write_name,
                  watershed_area_ha = ws_area_ha,
                  buffer_radius_m = as.numeric(rgx[, 2]),
                  snap_distance_m = as.numeric(rgx[, 4]),
                  snap_method = rgx[, 3],
                  dem_resolution = as.numeric(rgx[, 5]))

    return(deets)
}

get_gee_imgcol = function(gee_id, band, start, end, subset){

    if(! missing(subset)){

        gee_imgcol = ee$ImageCollection(gee_id)$
            filterDate(start, end)$
            select(band)$
            # select(subset)$
            map(function(x){
                date <- ee$Date(x$get("system:time_start"))$format('YYYY_MM_dd')
                x$set("RGEE_NAME", date)
            })

    } else {

        gee_imgcol = ee$ImageCollection(gee_id)$
            filterDate(start, end)$
            select(band)$
            map(function(x){
                date <- ee$Date(x$get("system:time_start"))$format('YYYY_MM_dd')
                x$set("RGEE_NAME", date)
            })
    }

    return(gee_imgcol)
}

get_gee = function(gee_id,
                   band,
                   start,
                   end,
                   res,
                   wb){

    imgcol = get_gee_imgcol(gee_id = gee_id,
                            band = band,
                            start = start,
                            end = end)

    Gout = list()

    Gout$Gmed = ee_extract(x = imgcol,
                           y = wb,
                           scale = res,
                           fun = ee$Reducer$median(),
                           # sf = FALSE)
                           sf = T)

    Gout$Gmean = ee_extract(x = imgcol,
                            y = wb,
                            scale = res,
                            fun = ee$Reducer$mean(),
                            sf = FALSE)

    Gout$Gstd = ee_extract(x = imgcol,
                           y = wb,
                           scale = res,
                           fun = ee$Reducer$stdDev(),
                           sf = FALSE)

    Gout$Gcnt = ee_extract(x = imgcol,
                           y = wb,
                           scale = res,
                           fun = ee$Reducer$count(),
                           sf = FALSE)

    Gout$Gmax = ee_extract(x = imgcol,
                           y = wb,
                           scale = res,
                           fun = ee$Reducer$max(),
                           sf = FALSE)

    Gout$Gmin = ee_extract(x = imgcol,
                           y = wb,
                           scale = res,
                           fun = ee$Reducer$min(),
                           sf = FALSE)

    Gnames = names(Gout)

    for(gn in Gnames){

        if(band %in% c('Lai_500m', 'Fpar_500m', 'Percent_Tree_Cover',
                       'Percent_NonTree_Vegetation',
                       'Percent_NonVegetated')){
            rgx = 'X([0-9]{4}_[0-9]{2}_[0-9]{2})_(.+)'
            # } else if(band %in% c('Lai_500m', 'Fpar_500m')){
            #     rgx = 'X([0-9]{4}_[0-9]{2}_[0-9]{2})_(.+)'
        } else {
            rgx = 'X([0-9]+)_(.+)'
        }

        z = pivot_longer(Gout[[gn]],
                         cols = everything(),
                         names_pattern = rgx,
                         names_to = c('date', '.value'))

        colnames(z)[2] = paste(colnames(z)[2],
                               substr(gn, 2, nchar(gn)),
                               sep = '_')

        Gout[[gn]] = z
    }

    Gout = purrr::reduce(Gout,
                         full_join,
                         by = 'date')

    return(Gout)
}

# zz = get_gee_imgcol(gee_id = 'UMT/NTSG/v2/LANDSAT/NPP',
#         band = 'annualNPP',
#         start = '1968-01-01',
#         end = '2020-12-31')
# ee_s2 <- ee$ImageCollection('UMT/NTSG/v2/LANDSAT/NPP')$
#     filterDate("2016-01-01", "2016-01-31")$
#     filterBounds(sf_as_ee(wb))
# ee_s2$size()$getInfo() # 126
# # Get the first 5 elements
# ee_get(ee_s2, index = 0:4)$size()$getInfo() # 5

get_gee_chunk = function(gee_id,
                         band,
                         start,
                         end,
                         res,
                         chunk_size_yr,
                         wb){

    if(chunk_size_yr < 2) stop('chunk_size_yr must be at least 2')

    startyr = lubridate::year(as.Date(start))
    endyr = lubridate::year(as.Date(end))
    # dif = as.numeric(difftime(start, end, units = 'days') / 365)
    nyrs = endyr - startyr + 1

    startseq = seq(startyr, endyr - chunk_size_yr, chunk_size_yr)
    # endseq = seq(startyr + chunk_size_yr, endyr, chunk_size_yr)
    # nchunks = floor(nyrs / chunk_size_yr)
    nchunks = length(startseq)

    geeout = tibble()
    lastyr_adj = 1
    for(i in 1:nchunks){
        try_err <<- FALSE
        if(i == nchunks && nyrs %% 2 == 1) lastyr_adj = 0
        starti = paste(startseq[i], '01-01', sep='-')
        endi = paste(startseq[i] + chunk_size_yr - lastyr_adj, '12-31', sep='-')
        # geechunk = try(get_gee(gee_id, band, start = starti, end = endi, res))
        # if('try-error' %in% class(geechunk){
        geechunk = tryCatch({
            get_gee(gee_id, band, start = starti, end = endi, res, wb = wb)
        }, error = function(e){
            try_err <<- TRUE
            print(paste(starti, 'to', endi, 'failed'))
        })
        if(try_err) next

        geeout = bind_rows(geeout, geechunk)
    }

    return(geeout)
}

comid_from_point = function(lat, long, CRS){
    pt = sf::st_point(c(long, lat))
    ptc = sf::st_sfc(pt, crs=CRS)
    COMID = nhdplusTools::discover_nhdplus_id(ptc)
    if(! length(COMID)) COMID = NA
    return(COMID)
}

choose_projection <- function(lat = NULL,
                              long = NULL,
                              unprojected = FALSE){

    if(unprojected){
        PROJ4 <- glue('+proj=longlat +datum=WGS84 +no_defs ',
                      '+ellps=WGS84 +towgs84=0,0,0')
        return(PROJ4)
    }

    if(is.null(lat) || is.null(long)){
        stop('If projecting, lat and long are required.')
    }

    abslat <- abs(lat)

    if(abslat < 23){ #tropical
        PROJ4 = glue('+proj=laea +lon_0=', long)
    } else { #temperate or polar
        PROJ4 = glue('+proj=laea +lat_0=', lat, ' +lon_0=', long)
    }

    return(PROJ4)
}


# 2. [obsolete] delineate the NHC watershed; load as shapefile (see section 10) ####

if(rebuild_all){
    deets <- delineate_watershed_from_point(lat = sites$latitude,
                                            long = sites$longitude,
                                            crs = EPSG,
                                            write_dir = 'data/watershed_boundary',
                                            write_name = 'NHC')
}

wb = sf::st_read('data/watershed_boundary/NHC.shp')

# 3. [obsolete] get PRISM annual precip data 1970-2013 for NHC watershed (via geoknife) ####

if(rebuild_all){

    stencil = simplegeom(as(wb, 'Spatial'))
    fabric = webdata(list(times = as.POSIXct(c('1968-01-01', '2020-03-31')),
                          url = 'https://cida.usgs.gov/thredds/dodsC/prism',
                          variables = 'ppt'))

    job = geoknife(stencil, fabric, wait = TRUE)
    # successful(job)
    # job = cancel(job)

    precip = result(job)
    write_csv(precip, 'data/prism/prism_raw.csv') #NOTE: this is monthly
}

precip = read_csv('data/prism/prism_raw.csv') %>%
    as_tibble() %>%
    dplyr::select(datetime = DateTime, precip_mm = '1') %>%
    group_by(year = substr(datetime, 1, 4)) %>%
    summarize(precip_mm = sum(precip_mm, na.rm = TRUE)) %>%
    ungroup()

# 4. [obsolete?] get PRISM data and localize it to NHC watershed (old fashioned way) ####

if(rebuild_all){

    dir.create(paste0(getwd(), '/data/prism_raw'),
               showWarnings = FALSE)
    options(prism.path = paste0(getwd(), '/data/prism_raw'))

    get_prism_annual(type = 'ppt',
                     year = 1968:2020,
                     keepZip = FALSE)


    # stencil = simplegeom(as(wb, 'Spatial'))
    # fabric = webdata(list(times = as.POSIXct(c('1968-01-01', '2020-03-31')),
    #                       url = 'https://cida.usgs.gov/thredds/dodsC/prism',
    #                       variables = 'ppt'))
    #
    # job = geoknife(stencil, fabric, wait = TRUE)
    # # successful(job)
    # # job = cancel(job)
    #
    # precip = result(job)
    # write_csv(precip, 'data/prism/prism_raw.csv') #NOTE: this is monthly
}

precip = read_csv('data/prism/prism_raw.csv') %>%
    as_tibble() %>%
    dplyr::select(datetime = DateTime, precip_mm = '1') %>%
    group_by(year = substr(datetime, 1, 4)) %>%
    summarize(precip_mm = sum(precip_mm, na.rm = TRUE)) %>%
    ungroup()

# 5. get NLDAS data 1970-present for NHC watershed ####

if(rebuild_all){

    webdatasets = query('webdata')
    set_ind = grep('nldas', abstract(webdatasets), ignore.case = TRUE)
    endpoint = url(webdatasets[set_ind])
    wbd = webdata(webdatasets[set_ind])

    stencil = simplegeom(as(wb, 'Spatial'))
    fabric = webdata(list(times = as.POSIXct(c('1968-01-01', '2020-12-31')),
                          url = endpoint,
                          variables = query(wbd, 'variables')))

    job = geoknife(stencil, fabric, wait = TRUE)
    # successful(job)
    # job = cancel(job)

    nldas = result(job)

    nldas = nldas %>%
        as_tibble() %>%
        pivot_wider(names_from = 'variable',
                    values_from = '1') %>%
        dplyr::select(datetime = DateTime, -statistic)
        # group_by(year = substr(datetime, 1, 4)) %>%
        # summarize(precip_mm = sum(precip_mm, na.rm = TRUE)) %>%
        # ungroup()
    write_csv(nldas, 'data/nldas.csv')
}


# 6. [impossible in R] get NWALT land cover 1974-2012 (from Arc layers) ####

#the layer takes about 58GB RAM to load (srsly), and then fails with several errors.

library(rgdal)

nwalt_paths = c('https://water.usgs.gov/GIS/dsdl/ds948_landuse1974.zip')

dest = 'data/nwalt/raw'
dir.create(dest,
           showWarnings = FALSE,
           recursive = TRUE)

# nwalt_paths->p
for(p in nwalt_paths){
    path_components = strsplit(p, '/')[[1]]
    fn = path_components[length(path_components)]
    zipf = paste(dest, fn, sep = '/')
    download.file(url = p,
                  destfile = zipf,
                  cacheOK = FALSE,
                  method = 'curl')
    unzip(zipf, exdir = dest)
    unlink(zipf)

    #temporary dead end here. layer requires 16GB memory.
    #if picking up here, first must modify the below chunk to locate path to
    #layer programmatically. then must figure out how to subset various land
    #cover types from the layer programmatically
    x = new("GDALReadOnlyDataset", 'data/nwalt/raw/nwalt_landuse_1974/lu1974_050815/')
    getDriver(x)
    getDriverLongName(getDriver(x))
    xx = asSGDF_GROD(x)#, offset=c(10000, 10000))
    r = raster::raster(xx)
    plot(r)
}

# 7. [pending] evapotranspiration from geoknife ####

#adapt section 5 code

# x. [incomplete] geoknife experimentation (looking for other datasets) ####

if(rebuild_all){

    stencil = simplegeom(as(wb, 'Spatial'))
    # fabric = webdata(list(times = as.POSIXct(c('1970-01-01', '2019-12-31')),
    #                       url = 'https://cida.usgs.gov/thredds/dodsC/prism',
    #                       variables = 'ppt'))
    # fabric = webdata('nldas')#daymet”, “stageiv”, “topowx”, “solar”, “metobs
    fabric = webdata('prism') # good
    fabric = webdata('nldas') #gone?
    fabric = webdata('iclus')
    fabric = webdata('daymet') #?
    fabric = webdata('stageiv') #?
    fabric = webdata('solar') #just solar radiation
    fabric = webdata('metobs') #ocean stuff?
    query(fabric, 'variables')

    fabric = webdata('iclus',
                     times = as.POSIXct(c('1970-01-01', '2019-12-31')))

    job = geoknife(stencil, fabric, wait = FALSE)
    running(job)
    successful(job)
    error(job)
    job = cancel(job)

    x = result(job)
    # write_csv(precip, 'data/prism/prism_raw.csv')
}

# precip = read_csv('data/prism/prism_raw.csv') %>%
#     as_tibble() %>%
#     dplyr::select(datetime = DateTime, precip_mm = '1') %>%
#     group_by(month = substr(datetime, 1, 7)) %>%
#     summarize(precip_mm = sum(precip_mm, na.rm = TRUE)) %>%
#     ungroup()

# 8. [incomplete] delineate additional watersheds and calculate areas ####

#must copy function parts from macrosheds for this to work on all watersheds.
#but then, it will fail on bridges unless we work that out

more_sites = readr::read_csv('data/NHCsite_metadata_forWatershedDelin.csv')
more_sites = more_sites %>%
    filter(! grepl('^USGS', sitecode)) %>%
    mutate(ws_area_manual = NA_real_)

deetslist = list()

for(i in 1:nrow(more_sites)){

    current_site = more_sites$sitecode[i]

    deets <- delineate_watershed_from_point(lat = more_sites$latitude[i],
                                            long = more_sites$longitude[i],
                                            crs = EPSG,
                                            write_dir = 'data/watershed_boundary/more_boundaries',
                                            write_name = current_site)

    more_sites[i, 'ws_area_manual'] = deets$watershed_area_ha / 100
    deetslist[[i]] = deets
}

saveRDS('data/watershed_boundary/more_boundaries/deets.rds')

# 9. get full watershed boundary from NHD; isolate riparian zone of study area ####


if(rebuild_all){

    more_sites = readr::read_csv('data/NHCsite_metadata_forWatershedDelin.csv')

    nhc_mouth = more_sites %>%
        filter(sitecode == 'NHC') %>%
        st_as_sf(coords = c('longitude', 'latitude'),
                 remove = FALSE,
                 crs = 4326)

    comid = comid_from_point(lat = nhc_mouth$latitude,
                             long = nhc_mouth$longitude,
                             CRS = 4326)

    nhc_lines = nhdplusTools::navigate_nldi(
        nldi_feature = list(featureSource = 'comid',
                            featureID = as.character(comid)),
        mode = 'UT',
        data_source = '',
        distance_km = 100
    )

    nhc_wb = streamstats::delineateWatershed(
        xlocation = nhc_mouth$longitude,
        ylocation = nhc_mouth$latitude,
        crs = 4326,
        includeparameters = 'true'
        # includeflowtypes = TRUE
    )

    #still broken...
    # streamstats::writeShapefile(watershed = nhc_wb,
    #                             layer = 'nhc_wb_streamstats',
    #                             dir = 'data/watershed_boundary',
    #                             what = 'boundary')

    #save streamstats watershed boundary as shapefile
    #streamstats::toSp and streamstats::writeShapefile are broken;
    #the bodies of those functions are extracted and modified below:
    tpf = tempfile(fileext='.geojson')
    streamstats::writeGeoJSON(nhc_wb, file=tpf, what='boundary')
    spatialdf = rgdal::readOGR(tpf)
    unlink(tpf)
    rgdal::writeOGR(spatialdf,
                    dsn = 'data/watershed_boundary/',
                    layer = 'nhc_wb_streamstats',
                    driver = 'ESRI Shapefile')

    dir.create('data/other_watershed_stuff',
               showWarnings = FALSE)

    #read back in as sf
    nhc_wb = st_read('data/watershed_boundary/nhc_wb_streamstats.shp')

    #plot all just to make sure it's legit
    # nhc_points = st_as_sf(more_sites, coords = c(x = 'longitude', y = 'latitude'))
    # mv(nhc_wb) + mv(nhc_lines) + mv(nhc_points)

    #isolate sampled reach; combine and buffer to represent riparian area
    proj = choose_projection(lat = nhc_mouth$latitude,
                             long = nhc_mouth$longitude)

    nhc_ripar = filter(nhc_lines,
                       nhdplus_comid %in% c(8895490, 8895362, 8895420, 8895440)) %>%
        st_combine() %>%
        st_transform(crs = proj) %>%
        st_buffer(dist = 250) %>%
        st_transform(crs = 4326)

    st_write(nhc_ripar,
             dsn = 'data/other_watershed_stuff',
             layer = 'riparian',
             driver = 'ESRI shapefile',
             # driver = 'GeoJSON')
             delete_layer = TRUE)
}

nhc_wb = st_read('data/watershed_boundary/nhc_wb_streamstats.shp')
nhc_ripar = st_read('data/other_watershed_stuff/riparian.shp')

# 10. summarize GEE layers for the whole watershed ####

#annual NPP
Gnpp = get_gee(gee_id = 'UMT/NTSG/v2/LANDSAT/NPP',
               band = 'annualNPP',
               start = '1968-01-01',
               end = '2020-12-31',
               res = 30,
               wb = nhc_wb)
write_csv(Gnpp, '~/git/papers/alice_nhc/data/gee/npp.csv')

#GPP
Ggpp = get_gee(gee_id = 'UMT/NTSG/v2/LANDSAT/GPP',
               band = 'GPP',
               start = '1968-01-01',
               end = '2020-12-31',
               res = 30,
               wb = nhc_wb)
write_csv(Ggpp, '~/git/papers/alice_nhc/data/gee/gpp.csv')

#LAI
Glai = get_gee(gee_id = 'MODIS/006/MOD15A2H',
               band = 'Lai_500m',
               start = '1968-01-01',
               end = '2020-12-31',
               res = 500,
               wb = nhc_wb)
write_csv(Glai, '~/git/papers/alice_nhc/data/gee/lai.csv')

#FPAR
# Gfpar = get_gee(gee_id = 'MODIS/006/MOD15A2H',
#                 band = 'Fpar_500m',
#                 start = '1968-01-01',
#                 end = '2020-12-31',
#                 res = 500,
#                 wb = nhc_wb)
# write_csv(Gfpar, '~/git/papers/alice_nhc/data/gee/fpar.csv')

#landcover (MODIS) [NEEDS WORK]
# Gtyp1 = get_gee(gee_id = 'MODIS/006/MCD12Q1',
#                 band = 'LC_Type1',
#                 start = '1968-01-01',
#                 end = '2020-12-31',
#                 res = 500,
#                 wb = nhc_wb)

#tree cover (MODIS)
# Gtree = get_gee(gee_id = 'MODIS/006/MOD44B',
#                 band = 'Percent_Tree_Cover',
#                 start = '1968-01-01',
#                 end = '2020-12-31',
#                 res = 250,
#                 wb = nhc_wb)
# write_csv(Gtree, '~/git/papers/alice_nhc/data/gee/treed.csv')

# #non-tree veg cover (MODIS)
# Gveg = get_gee(gee_id = 'MODIS/006/MOD44B',
#                band = 'Percent_NonTree_Vegetation',
#                start = '1968-01-01',
#                end = '2020-12-31',
#                res = 250,
#                wb = nhc_wb)
# write_csv(Gveg, '~/git/papers/alice_nhc/data/gee/nontree_veg.csv')

# #vegless cover (MODIS)
# Gbare = get_gee(gee_id = 'MODIS/006/MOD44B',
#                 band = 'Percent_NonVegetated',
#                 start = '1968-01-01',
#                 end = '2020-12-31',
#                 res = 250,
#                 wb = nhc_wb)
# write_csv(Gbare, '~/git/papers/alice_nhc/data/gee/vegless.csv')

#precip (PRISM)
Gppt = get_gee_chunk(gee_id = 'OREGONSTATE/PRISM/AN81d',
                     band = 'ppt',
                     start = '1968-01-01',
                     end = '2020-12-31',
                     res = 4000,
                     chunk_size_yr = 2,
                     wb = nhc_wb)
write_csv(Gppt, '~/git/papers/alice_nhc/data/gee/ppt.csv')

#mean temp (PRISM)
Gtmean = get_gee_chunk(gee_id = 'OREGONSTATE/PRISM/AN81d',
                     band = 'tmean',
                     start = '1968-01-01',
                     end = '2020-12-31',
                     res = 4000,
                     chunk_size_yr = 2,
                     wb = nhc_wb)
write_csv(Gtmean, '~/git/papers/alice_nhc/data/gee/tmean.csv')

#max temp (PRISM)
Gtmax = get_gee_chunk(gee_id = 'OREGONSTATE/PRISM/AN81d',
                     band = 'tmax',
                     start = '1968-01-01',
                     end = '2020-12-31',
                     res = 4000,
                     chunk_size_yr = 2,
                     wb = nhc_wb)
write_csv(Gtmax, '~/git/papers/alice_nhc/data/gee/tmax.csv')

#min temp (PRISM)
Gtmin = get_gee_chunk(gee_id = 'OREGONSTATE/PRISM/AN81d',
                     band = 'tmin',
                     start = '1968-01-01',
                     end = '2020-12-31',
                     res = 4000,
                     chunk_size_yr = 2,
                     wb = nhc_wb)
write_csv(Gtmin, '~/git/papers/alice_nhc/data/gee/tmin.csv')

# 11. summarize GEE layers just for the riparian zone ####

#annual NPP
Gnpp_ripar = get_gee(gee_id = 'UMT/NTSG/v2/LANDSAT/NPP',
               band = 'annualNPP',
               start = '1968-01-01',
               end = '2020-12-31',
               res = 30,
               wb = nhc_ripar)
write_csv(Gnpp_ripar, '~/git/papers/alice_nhc/data/gee/npp_riparian.csv')

#GPP
Ggpp_ripar = get_gee(gee_id = 'UMT/NTSG/v2/LANDSAT/GPP',
               band = 'GPP',
               start = '1968-01-01',
               end = '2020-12-31',
               res = 30,
               wb = nhc_ripar)
write_csv(Ggpp_ripar, '~/git/papers/alice_nhc/data/gee/gpp_riparian.csv')

#LAI
Glai_ripar = get_gee(gee_id = 'MODIS/006/MOD15A2H',
               band = 'Lai_500m',
               start = '1968-01-01',
               end = '2020-12-31',
               res = 500,
               wb = nhc_ripar)
write_csv(Glai_ripar, '~/git/papers/alice_nhc/data/gee/lai_riparian.csv')

# x. [testing] messing with GEE. trying to get nlcd ####

zz = ee$FeatureCollection('USGS/NLCD')
subset = zz$filterBounds(wb)
rgee$
ee_as_sf(
    # filterDate(start, end)$

gee_imgcol = ee$ImageCollection(gee_id)$
    filterDate(start, end)$
    select(band)$
    map(function(x){
        date <- ee$Date(x$get("system:time_start"))$format('YYYY_MM_dd')
        x$set("RGEE_NAME", date)
    })

wb_ee = sf_as_ee(wb)
# wb_ee$geometry()

lcv <- ee$FeatureCollection('USGS/NLCD')$select('landcover')$
    filterBounds(wb_ee) %>%
    ee_as_sf()

lcv <- ee$ImageCollection('USGS/NLCD')$select('landcover')$
    filter(ee$Filter$eq('system:index', 'NLCD2011'))$first() %>%
    ee$Image$clip(wb_ee) %>%
    ee_as_raster(maxPixels = 20000000000)
# imp <- ee$ImageCollection('USGS/NLCD')$select('impervious') %>%
    ee_as_stars(maxPixels = 20000000000)
lcvv = lcv %>%
    sf::st_as_sf()
mapview::mapview(lcv)
# #impervious surface (NLCD)
# Gimperv = get_gee(gee_id = 'USGS/NLCD',
                  # band = 'impervious',
                  # start = '1968-01-01',
                  # end = '2020-12-31',
                  # res = 30)


#
# #tree cover (NLCD)
# Gtree2 = get_gee(gee_id = 'USGS/NLCD',
#                  band = '24',
#                  # band = 'landcover',
#                  # band = 'percent_tree_cover',
#                  start = '1968-01-01',
#                  end = '2020-12-31',
#                  res = 30)

# take 2 ####

# wdpa<-ee$FeatureCollection("WCMC/WDPA/current/points")
# columns<-wdpa$first()$getInfo()
# names(columns$properties)

# library(geojsonsf)

# nhc_wb_geojson = geojsonsf::sf_geojson(nhc_wb)
wb_ee = sf_as_ee(nhc_wb)

donkey = ee$ImageCollection('USGS/NLCD')
chili = donkey$first()$getInfo()
names(chili$properties)

gee_imgcol = ee$ImageCollection('USGS/NLCD')$
    # filterDate(start, end)$
    select('landcover')$
    filterBounds(wb_ee)$
    # filterDate('1968-01-01', '2020-12-31')$
    ee_print()
    # select(subset)$
    map(function(x){
        date <- ee$Date(x$get("system:time_start"))$format('YYYY_MM_dd')
        x$set("RGEE_NAME", date)
    })

point <- ee$Geometry$Point(-44.366,-18.145)
start <- ee$Date("2019-07-11")
end <- ee$Date("2019-07-20")
col.filt<-col$filterBounds(point)$filterDate(start,end)

col<-ee$ImageCollection('LANDSAT/LT05/C01/T2')$
    filterDate('1987-01-01','1990-05-01')$
    filterBounds(ee$Geometry$Point(25.8544,-18.08874))
filtered <- col$filterMetadata('IMAGE_QUALITY','equals',9)
notSoGood<-ee$Algorithms$Landsat$simpleComposite(col,75,3)
good<-ee$Algorithms$Landsat$simpleComposite(filtered,75,3)
Map$setCenter(25.8544,-18.08874,13)
Map$addLayer(notSoGood, list(bands=c('B3','B2','B1'), gain=3.5), 'Bad Composition')
Map$addLayer(good,list(bands=c('B3', 'B2', 'B1'), gain=3.5),'Good Composition')

Gout$Gmed = ee_extract(x = imgcol,
                       y = wb,
                       scale = res,
                       fun = ee$Reducer$median(),
                       # sf = FALSE)
                       sf = T)

#take 3####

wb_ee = sf_as_ee(nhc_wb)
img = ee$ImageCollection('USGS/NLCD')$select('landcover')$
    filter(ee$Filter$eq('system:index', 'NLCD2011'))$first()$clip(wb_ee)
dir.create('data/nlcd', showWarnings = FALSE)
rst = ee_as_raster(image = img,
                   region = wb_ee$geometry(),
                   dsn = 'data/nlcd/NLCD2011.tif')


ftr = ee$FeatureCollection('USGS/NLCD')$select('landcover')$
    filterBounds(wb_ee)
    filter(ee$Filter$eq('system:index', 'NLCD2011'))$first()

Map$centerObject(wb_ee)
Map$addLayer(eeObject = ee$FeatureCollection(wb_ee)) +
Map$addLayer(
    eeObject = wb_ee,
    visParams = list(palette = 'yellow'),
    name = 'ROI'
)

Map$centerObject(clipped)
Map$addLayer(
    eeObject = clipped,
    visParams = list(palette = 'red'),
    name = 'clipped') +
Map$addLayer(
    eeObject = ftr,
    name = 'Census roads'
)

    ee_as_raster()
    # ee_as_raster(maxPixels = 20000000000)
    # ee_as_sf()
    # ee_as_stars(maxPixels = 20000000000)
    # sf::st_as_sf()
clipped <- ftr$map(function(x) x$intersection(wb_ee))
mapview::mapview(rst)

# fc <- ee$FeatureCollection('TIGER/2016/Roads')$filterBounds(roi)
# clipped <- fc$map(function(x) x$intersection(roi))


# vPar <- list(bands = c('landcover'))#,min = 100,max = 8000)
# Map$setCenter(-44.366,-17.69, zoom = 10)
# Map$addLayer(img, vPar, "True Color Image")

# zz = img %>%
#     rgee::ee_extract(y = nhc_wb,
#                      scale = 30,
#                      fun = ,
#                      sf = TRUE)
# rgee::ee_as_sf()

Map$centerObject(wb_ee)
Map$addLayer(
    eeObject = img,
    visParams = list(bands = c('landcover'), min = 0, max = 95),
    name = 'nlcd'
)

Map$addLayer(
    eeObject = wb_ee,
    name = 'NHC'
)

# Get a dictionary of means in the region.  Keys are bandnames.
nlcd_sums <- img$reduceRegion(
    reducer = ee$Reducer$sum(),
    geometry = wb_ee,
    scale = 30
)

print(nlcd_sums$getInfo())


#for reals ####

dir.create('data/nlcd', showWarnings = FALSE)

wb_ee = sf_as_ee(nhc_wb)

img = ee$ImageCollection('USGS/NLCD')$select('landcover')$
    filter(ee$Filter$eq('system:index', 'NLCD2011'))$first()$clip(wb_ee)

rst = ee_as_raster(image = img,
                   region = wb_ee$geometry(),
                   dsn = 'data/nlcd/NLCD2011.tif')

values(rst)
nlcd2011 = raster::raster('data/nlcd/NLCD2011.tif')
