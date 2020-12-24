#===============================================================================
#Diagnostic plots for streams
#Created 4/7/2017
#Modified 11/13/2017
#===============================================================================
library(ks) #For making kernel density plots
library(tidyr)
library(dplyr)
library(zoo)
library(scales)
library(readr)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl")
ylims=c(-15,5)
GPP.col="#74B16D"
ER.col="#B99D82"
PAR.col = "#FFE083"
hall_met <- read_csv("data/hall/hall_table_15.csv")
hall_met$ER_gO2m2d <- -hall_met$ER_gO2m2d
# NHCdat <- readRDS("../NHC_2019_metabolism/data/metabolism/condensed/allNHCsites.rds")
# metab <- NHCdat$metab
# metab[!is.na(metab$GPP)&metab$GPP<0,c("GPP","GPP.upper","GPP.lower")]<-NA
# metab[!is.na(metab$ER)&metab$ER>0,c("ER","ER.upper","ER.lower")]<-NA
metab <- read_rds("data/NHC_metab_allsites.rds")

# data <- NHCdat$data
sites <- c("UNHC","WBP","WB","CBP","PM","NHC")
# sitematch <- data.frame(hall = c("Concrete", "Blackwood", "Wood Bridge"), 
#                         sitecode= c("CBP", "UNHC", "WB"))



# pdf("figures/NHC2019diagnostics3.pdf",width = 7, height = 5.83,onefile = TRUE)
png(width=7, height=6, units='in', type='cairo', res=300,
    filename="figures/metab_contour_comparison.png")    
    kplot = 5

    kernel <- kde(na.omit(metab[,c("GPP","ER")]))
    plot(kernel, xlab = "GPP", ylab = "ER", ylim = c(-kplot, 0), xlim = c(0, kplot),
         cont=c(30,60,90), col="grey35", lwd = 2,
         main = "All NHC sites metabolism")
    abline(0,-1)
    
    kernel_hall <- kde(na.omit(hall_met[,c("GPP_gO2m2d","ER_gO2m2d")]))
    par(new=T)
    plot(kernel_hall, xlab = "", ylab = "", ylim = c(-kplot, 0), xlim = c(0, kplot),
         cont = c(30,60,90), col = "darkred", lwd=2)
    legend("topright", cex = 1.4,
           c(paste0("2019-2020  n = ",nrow(na.omit(metab))),
             paste0("1968-1970  n = ", nrow(na.omit(hall_met)))),
           col = c("grey35", "darkred"), lty = 1, lwd = 2, bty = "n")
    #points(hall_met$GPP_gO2m2d, -hall_met$ER_gO2m2d, col="brown3", pch=20, cex=1.7)
dev.off()        




