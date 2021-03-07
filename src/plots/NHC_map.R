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
stream_line <- st_read("data/map_files/study_reaches.shp")
riparian_boundary <- st_read("data/other_watershed_stuff/riparian.shp")
study_reaches_line <- st_read("data/other_watershed_stuff/riparian.shp")
duke_forest_boundary <- st_read("data/map_files/2019_boundary.shp")
korstian_div <- filter(duke_forest_boundary, DIVISION == 'Korstian')
watershed_boundary <- st_read("data/watershed_boundary/nhc_wb_streamstats.shp")

tmap_mode("view")

tm_shape(watershed_boundary) + tm_polygons(alpha=0, border.col="black", lwd=.5) +
    tm_shape(riparian_boundary) + tm_polygons(alpha=0.1, col="black",
                                              border.col='black', border.alpha = 0.1) +
    tm_shape(korstian_div) + tm_polygons(alpha=0.2, col = 'darkgreen',
                                         border.col="transparent", lwd=.5) +
    tm_shape(study_reaches_line) + tm_lines(col='gray60', lwd=1.5) +
    tm_shape(stream_line) + tm_lines(col='gray40', lwd=1) +
    tm_shape(sites_sf) + tm_dots(col="brown3", size=.05) +
    tm_shape(wwtp_sf) + tm_markers(shape=3, col="lightblue", size=.05) +
    tm_scale_bar(text.size = 1, position = "left")

tmap_mode("plot")
par(bg=NA)
map<-tm_shape(nhc_ws)+tm_polygons(alpha=0, border.col="black", lwd=.5)+
#  tm_shape(mud_ws)+tm_polygons(alpha=0, border.col="black",lwd=.5)+
  tm_shape(cur_nhd)+tm_lines(col = "grey60") +
  tm_shape(longitudinal_transect) + tm_lines(lwd=2)+
  tm_shape(sites_sf)+tm_dots(col="brown3", size=.05)+
#  tm_shape(wwtp_sf)+tm_markers(shape=3, col="lightblue",size=.05)+
  tm_scale_bar(text.size = 1, position = "left") +
  tm_compass(type="arrow",position=c("right","bottom", show.labels=3))+
  tm_layout(frame=FALSE, bg.color="transparent")
tmap_save(map, filename="NHCmap_scalebar.eps", bg="transparent", dpi = 1200,
          )


# Plot of longitudinal transect
long_sites_sf <- sites_sf[sites_sf$site!="MC751",]
tmap_mode("view")

tm_shape(cur_nhd)+tm_lines(col = "grey80") +
  tm_shape(longitudinal_transect) + tm_lines(lwd=2)+
  tm_shape(long_sites_sf)+tm_dots(col="brown3", size=.05)+
  #  tm_shape(wwtp_sf)+tm_markers(shape=3, col="lightblue",size=.05)+
  tm_scale_bar(text.size = 1, position = "left")

# Plot North carolina with piedmont shape


pied <- readOGR(dsn="ncpiedmont_shape",
                layer="Piedmont_shape")

par(bg=NA)
png("NCmap.png",bg="transparent", type="windows")
map('state',region='North Carolina',fill=TRUE, col="white",bg="transparent",lwd=2)
  plot(pied,
     add=TRUE,
     col="grey90")
  points(wwtp_sf$Long, wwtp_sf$Lat, col="brown3", pch=22, cex=3)
dev.off()
