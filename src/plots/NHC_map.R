library(nhdplusTools)
library(tidyverse)
library(sf)
library(tmap)
library(rgdal)
library(maps)

setwd('~/git/papers/alice_nhc')

sites <- read.csv("data/map_files/NHCsite_50yl_coordinates.csv",
                  header=TRUE, stringsAsFactors=FALSE)

sites_sf <- st_as_sf(sites, coords=c("longitude", "latitude"), remove=FALSE, crs=4326)
carter_sites <- filter(sites_sf, type == 'now')
hall_sites <- filter(sites_sf, type == 'hall')
stream_line <- st_read("data/map_files/stream_lines.shp")
# riparian_boundary <- st_read("data/other_watershed_stuff/riparian.shp")
study_reaches_line <- st_read("data/map_files/study_reaches.shp")
duke_forest_boundary <- st_read("data/map_files/2019_boundary.shp")
korstian_div <- filter(duke_forest_boundary, DIVISION == 'Korstian')
watershed_boundary <- st_read("data/watershed_boundary/nhc_wb_streamstats.shp")

#remove the piece of the streamlines that's downstream of the lowest site
stream_line_shortened <- stream_line
ll = as(stream_line_shortened$geometry[[1]], 'list')
ll[1][[1]] = ll[1][[1]][1:21, ]
stream_line_shortened$geometry[[1]] = st_multilinestring(ll)

study_reaches_line_shortened <- study_reaches_line
ll = as(study_reaches_line_shortened$geometry[[1]], 'list')
ll[1][[1]] = ll[1][[1]][1:21, ]
study_reaches_line_shortened$geometry[[1]] = st_multilinestring(ll)

#make new riparian boundary
PROJ4 = paste0('+proj=laea +lat_0=', round(mean(sites$latitude), 4), ' +lon_0=',
               round(mean(sites$longitude), 4),
             ' +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
riparian_boundary = study_reaches_line_shortened %>%
    st_transform(crs = PROJ4) %>%
    st_buffer(dist = 250) %>%
    st_transform(crs = 4326)

# tmap_mode("view")
tmap_mode("plot")
# par(bg=NA)

map_with_colon = tm_shape(watershed_boundary) + tm_polygons(alpha=0, border.col="black", lwd=1) +
    tm_shape(korstian_div) + tm_polygons(alpha=0.3, col = 'springgreen3',
                                         border.col="transparent", lwd=.5) +
    tm_shape(riparian_boundary) + tm_polygons(alpha=0, col="black", lwd=1.5,
                                              border.col='steelblue3', border.alpha=0.8) +
    tm_shape(study_reaches_line_shortened) + tm_lines(col='steelblue3', lwd=2.5) +
    tm_shape(stream_line_shortened) + tm_lines(col='black', alpha=0.5, lwd=0.5) +
    tm_shape(carter_sites) + tm_symbols(shape=1, col="red2", size=0.6, border.lwd=2) +
    tm_shape(hall_sites) + tm_symbols(shape=3, col="black", size=0.6, border.lwd=2) +
    tm_scale_bar(text.size = 1, position="left") +
    tm_compass(type="arrow", position=c("right", "bottom", show.labels=3),
               size=5, text.size=1) +
    tm_style(style='white') +
    tm_layout(frame=TRUE, bg.color="white") +
    tm_add_legend(type='symbol', labels = '  Study sites', col = 'red2', size = 0.7,
                  shape=1) +
    tm_add_legend(type='symbol', labels = '  Hall 1972 sites', col = 'black',
                  size=0.5, shape=3, border.lwd=2) +
    tm_add_legend(type='line', labels = 'Study reach', col = 'steelblue3', lwd = 2.5) +
    tm_add_legend(type='line', labels = 'Riparian zone', col = 'steelblue3', lwd = 1) +
    tm_add_legend(type='fill', labels = 'Duke Forest', col = 'springgreen3', alpha=0.3,
                  border.col='transparent') +
    tm_legend(show=TRUE, position=c('left', 'top'), outside=FALSE, bg.color='gray97',
              frame=TRUE, text.size=1.1)

map_without_colon = tm_shape(watershed_boundary) + tm_polygons(alpha=0, border.col="black", lwd=1) +
    tm_shape(korstian_div) + tm_polygons(alpha=0.3, col = 'springgreen3',
                                         border.col="transparent", lwd=.5) +
    tm_shape(study_reaches_line_shortened) + tm_lines(col='steelblue3', lwd=2.5) +
    tm_shape(stream_line_shortened) + tm_lines(col='black', alpha=0.5, lwd=0.5) +
    tm_shape(carter_sites) + tm_symbols(shape=1, col="red2", size=0.6, border.lwd=2) +
    tm_scale_bar(text.size = 1, position="left") +
    tm_compass(type="arrow", position=c("right", "bottom", show.labels=3),
               size=5, text.size=1) +
    tm_style(style='white') +
    tm_layout(frame=TRUE, bg.color="white") +
    tm_add_legend(type='symbol', labels = '  Study sites', col = 'red2', size = 0.7,
                  shape=1) +
    tm_add_legend(type='line', labels = 'Study reach', col = 'steelblue3', lwd = 2.5) +
    tm_add_legend(type='fill', labels = 'Duke Forest', col = 'springgreen3', alpha=0.3,
                  border.col='transparent') +
    tm_legend(show=TRUE, position=c('left', 'top'), outside=FALSE, bg.color='gray97',
              frame=TRUE, text.size=1.1)

tmap_save(map_without_colon, filename="figs/map_without_colon.png",
          bg="white", dpi = 300)

tmap_save(map_with_colon, filename="figs/map_with_colon.png",
          bg="white", dpi = 300)
