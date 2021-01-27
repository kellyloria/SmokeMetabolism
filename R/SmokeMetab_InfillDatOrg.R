#' Data Organization for modeling metabolism for signals of smoke effect 
#' 
#' @description Infilling of NA's in data downloaded USGS and Mesowest 
#'@details Last update: 2020-20-21
#' See http://usgs-r.github.io/streamMetabolizer for vignettes on the web.

## ---------------------------
# File path setup:
if (dir.exists('/Users/kellyloria/Documents/UNR_2020/Fall2020Projects/')){
  inputDir<- '/Users/kellyloria/Documents/UNR_2020/Fall2020Projects/'
  outputDir<- '/Users/kellyloria/Documents/UNR_2020/Fall2020Projects/SmokeMetabolism/SmokeMetabData'
}
## ---------------------------

# Load packages 
library(StreamMetabolism)
library(streamMetabolizer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(rstan)
library(unitted)
library(zoo)
library(lubridate)

## ---------------------------
# SANTA MARGARITA R NR TEMECULA
## ---------------------------
# Load in data:
# https://www.dropbox.com/s/ujyad748ka7flli/nwis.11044000.csv?dl=0
# change dl=0 to dl=1
dat1 <- read.csv( "https://www.dropbox.com/s/ujyad748ka7flli/nwis.11044000.csv?dl=1" )
summary(dat1) # SANTA MARGARITA R NR TEMECULA CA https://waterdata.usgs.gov/nwis/inventory/?site_no=11044000&agency_cd=USGS
range(dat1$datetime)
dat1$timestampPDT <- as.POSIXct(dat1$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
dat1$timestamp <- as.POSIXct(round_date(dat1$timestampPDT, hour,unit="5 minutes"))

# left join in Mesowest data
clim2 < read.csv("https://www.dropbox.com/s/gs8w50soclliami/KPSP.2020-10-01.csv?dl=1")
clim2$timestampUTC <- as.POSIXct(clim2$date_time, format = "%Y-%m-%dT%H:%M", tz="UTC")
clim2$timestamp <- as.POSIXct(round_date(clim2$timestampUTC, hour,unit="5 minutes"))
summary(clim2)

# Merge data files
TemeculaDat <- left_join(dat1, clim2[c("timestamp", "pressure_set_1d", "air_temp_set_1")],
                         by = c("timestamp" = "timestamp"))
## 
## Infill missing data:
## 

summary(TemeculaDat)
# Infill some data with rolling averages:
# Water Temp
TemeculaDat$tempCMean <- rollapply(TemeculaDat$tempC, width=500,
                                   FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                   by.column=TRUE, partial=TRUE, fill=NA, align="center")

TemeculaDat$tempC<- as.numeric(ifelse(is.na(TemeculaDat$tempC),
                                      (TemeculaDat$tempCMean),
                                      (TemeculaDat$tempC)))

# Discharge:
TemeculaDat$dischargeMean <- rollapply(TemeculaDat$discharge, width=500,
                                       FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                       by.column=TRUE, partial=TRUE, fill=NA, align="center")

TemeculaDat$discharge <- as.numeric(ifelse(is.na(TemeculaDat$discharge), 
                                           (TemeculaDat$dischargeMean),
                                           (TemeculaDat$discharge)))
# Dissolved oxygen:
TemeculaDat$DOmgLMean <- rollapply(TemeculaDat$DOmgL, width=500,
                                   FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                   by.column=TRUE, partial=TRUE, fill=NA, align="center")

TemeculaDat$DOmgL<- as.numeric(ifelse(is.na(TemeculaDat$DOmgL), 
                                      (TemeculaDat$DOmgLMean),
                                      (TemeculaDat$DOmgL)))

# Pressure:
# Date range for both data restrict to 2016 for missing weather obs
TemeculaDat1 <- subset(TemeculaDat, timestamp >= as.POSIXct('2016-10-01 00:53:00') & 
                         timestamp <= as.POSIXct('2020-09-10 00:00:00'))

TemeculaDat1$pressureMean <- rollapply(TemeculaDat1$pressure_set_1d, width=500,
                                      FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                      by.column=TRUE, partial=TRUE, fill=NA, align="center")

TemeculaDat1$pressure_set_1d <- as.numeric(ifelse(is.na(TemeculaDat1$pressure_set_1d), 
                                                  (TemeculaDat1$pressureMean),
                                                  (TemeculaDat1$pressure_set_1d)))
summary(TemeculaDat1)

TemeculaDat1$pressure_millibar <- c(TemeculaDat1$pressure_set_1d * 0.01)


### Stream Metabolizer estimates ###
TemeculaDat1$depth <- calc_depth(Q=u(TemeculaDat1$discharge, "m^3 s^-1"), f=u(0.36))
latitude <- c(33.4739175)
longitude <- c(-117.1422536) 

TemeculaDat1$solar.time <- calc_solar_time(TemeculaDat1$timestamp, longitude)
TemeculaDat1$light <- calc_light(TemeculaDat1$solar.time,
                                latitude,
                                longitude,
                                max.PAR = u(2326, "umol m^-2 s^-1"),
                                attach.units = is.unitted(TemeculaDat1$solar.time))

TemeculaDat1$DO.sat <- calc_DO_sat(TemeculaDat1$tempC, 
                                  TemeculaDat1$pressure_millibar, sal=0) 

names(TemeculaDat1)
# colnames(TemeculaDat1)[4] <- "temp.water"
# colnames(TemeculaDat1)[8] <- "DO.obs"

TemeculaDatQ <- subset(TemeculaDat1, select= c(solar.time, DO.obs, DO.sat, depth, temp.water, light, discharge))

# write.csv(TemeculaDatQ, paste0(outputDir, "/11044000TemeculaDat_Infill.csv")) 



## ---------------------------
# SANTA YNEZ R A SOLVANG CA
## ---------------------------
dat2 <- read.csv("https://www.dropbox.com/s/eochnr2z8v3qsip/nwis.11128500.csv?dl=1")
dat2$timestampPDT <- as.POSIXct(dat2$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
dat2$timestamp <- as.POSIXct(round_date(dat2$timestampPDT, hour,unit="5 minutes"))

# Mesowest dtat near (34.58498706, -120.1445926)
KIZAdat <- read.csv("https://www.dropbox.com/s/6v1w2ubyx2hdpdo/KIZA.2020-10-01.csv?dl=1")
# Fix timestamp (round to nearest 5 minute):
KIZAdat$timestampUTC <- as.POSIXct(KIZAdat$Date_Time, format = "%Y-%m-%dT%H:%M", tz="UTC")
KIZAdat$timestamp <- as.POSIXct(round_date(KIZAdat$timestampUTC, hour,unit="5 minutes"))

# Merge data files 
SantaYnezDat <- left_join(dat2, KIZAdat[c("timestamp", "altimeter_set_1", "pressure_set_1d", "air_temp_set_1")],
                          by = c("timestamp" = "timestamp"))

SantaYnezDat <- subset(SantaYnezDat, timestamp >= as.POSIXct('2015-10-01 00:53:00') & 
                         timestamp <= as.POSIXct('2020-10-01 00:00:00'))
range(SantaYnezDat$timestamp)

## 
## Infill missing data: 
##    ** A lot of gaps in this data set **
##        - Especially temperature 

# Infill some data with rolling averages:
ggplot(SantaYnezDat, aes(x=timestamp, y=(DOmgL))) +
  geom_point(alpha = 0.7)  +
  theme_classic()

# Water Temp
SantaYnezDat$tempCMean <- rollapply(SantaYnezDat$tempC, width=500,
                                   FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                   by.column=TRUE, partial=TRUE, fill=NA, align="center")

SantaYnezDat$tempC<- as.numeric(ifelse(is.na(SantaYnezDat$tempC),
                                      (SantaYnezDat$tempCMean),
                                      (SantaYnezDat$tempC)))

# Discharge:
SantaYnezDat$dischargeMean <- rollapply(SantaYnezDat$discharge, width=500,
                                       FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                       by.column=TRUE, partial=TRUE, fill=NA, align="center")

SantaYnezDat$discharge <- as.numeric(ifelse(is.na(SantaYnezDat$discharge), 
                                           (SantaYnezDat$dischargeMean),
                                           (SantaYnezDat$discharge)))
# Dissolved oxygen:
SantaYnezDat$DOmgLMean <- rollapply(SantaYnezDat$DOmgL, width=500,
                                   FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                   by.column=TRUE, partial=TRUE, fill=NA, align="center")

SantaYnezDat$DOmgL<- as.numeric(ifelse(is.na(SantaYnezDat$DOmgL), 
                                      (SantaYnezDat$DOmgLMean),
                                      (SantaYnezDat$DOmgL)))

# Pressure:
SantaYnezDat$pressureMean <- rollapply(SantaYnezDat$altimeter_set_1_pascals, width=500,
                                       FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                       by.column=TRUE, partial=TRUE, fill=NA, align="center")

SantaYnezDat$altimeter_set_1_pascals <- as.numeric(ifelse(is.na(SantaYnezDat$altimeter_set_1_pascals), 
                                                  (SantaYnezDat$pressureMean),
                                                  (SantaYnezDat$altimeter_set_1_pascals)))
summary(SantaYnezDat)

SantaYnezDat$pressure_millibar <- c(SantaYnezDat$pressure_set_1d * 0.01)


### Stream Metabolizer estimates ###
SantaYnezDat$depth <- calc_depth(Q=u(SantaYnezDat$discharge, "m^3 s^-1"), f=u(0.36))

#(34.58498706, -120.1445926)
latitude <- c(34.58498706)
longitude <- c(-120.1445926) 

SantaYnezDat$solar.time <- calc_solar_time(SantaYnezDat$timestamp, longitude)
SantaYnezDat$light <- calc_light(SantaYnezDat$solar.time,
                                 latitude,
                                 longitude,
                                 max.PAR = u(2326, "umol m^-2 s^-1"),
                                 attach.units = is.unitted(SantaYnezDat$solar.time))

SantaYnezDat$DO.sat <- calc_DO_sat(SantaYnezDat$tempC, 
                                   SantaYnezDat$pressure_millibar, sal=0) 
names(SantaYnezDat)
# colnames(SantaYnezDat)[7] <- "temp.water"
# colnames(SantaYnezDat)[9] <- "DO.obs"

SantaYnezDatQ <- subset(SantaYnezDat, select= c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light", "discharge"))

# write.csv(SantaYnezDatQ, paste0(outputDir, "/11128500SantaYnezDat_Infill.csv")) 



## ---------------------------
# SAN JOAQUIN R AB MERCED R NR NEWMAN CA
## ---------------------------
dat3 <- read.csv("https://www.dropbox.com/s/d1taddci0q0lvqk/nwis.11273400.csv?dl=1")
summary(dat3)
dat3$timestampPDT <- as.POSIXct(dat3$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
dat3$timestamp <- as.POSIXct(round_date(dat3$timestampPDT, hour,unit="5 minutes"))
summary(dat3)

#Mesowest data near (37.3472151, -120.9761777)
KMOD <- read.csv("https://www.dropbox.com/s/icjidfpe5f44y5i/KMOD.2020-10-01.csv?dl=1")
summary(KMOD)

##
# Fix timestamp (round to nearest 5 minute):
KMOD$timestampUTC <- as.POSIXct(KMOD$Date_Time, format = "%Y-%m-%dT%H:%M", tz="UTC")
KMOD$timestamp <- as.POSIXct(round_date(KMOD$timestampUTC, hour,unit="5 minutes"))

# Merge data files 
SanJoaquin <- left_join(dat3, KMOD[c("timestamp", "air_temp_set_1", "pressure_set_1d")],
                        by = c("timestamp" = "timestamp"))

SanJoaquin <- subset(SanJoaquin, timestampPDT >= as.POSIXct('2015-10-01 00:53:00') & 
                       timestampPDT <= as.POSIXct('2020-10-01 00:00:00'))
range(SanJoaquin$timestampPDT)
summary(SanJoaquin)

## 
## Infill missing data: 
##   
# Water Temp
SanJoaquin$tempCMean <- rollapply(SanJoaquin$tempC, width=500,
                                    FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                    by.column=TRUE, partial=TRUE, fill=NA, align="center")

SanJoaquin$tempC<- as.numeric(ifelse(is.na(SanJoaquin$tempC),
                                       (SanJoaquin$tempCMean),
                                       (SanJoaquin$tempC)))

# Discharge:
SanJoaquin$dischargeMean <- rollapply(SanJoaquin$discharge, width=500,
                                        FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                        by.column=TRUE, partial=TRUE, fill=NA, align="center")

SanJoaquin$discharge <- as.numeric(ifelse(is.na(SanJoaquin$discharge), 
                                            (SanJoaquin$dischargeMean),
                                            (SanJoaquin$discharge)))
# Dissolved oxygen:
SanJoaquin$DOmgLMean <- rollapply(SanJoaquin$DOmgL, width=500,
                                    FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                    by.column=TRUE, partial=TRUE, fill=NA, align="center")

SanJoaquin$DOmgL<- as.numeric(ifelse(is.na(SanJoaquin$DOmgL), 
                                       (SanJoaquin$DOmgLMean),
                                       (SanJoaquin$DOmgL)))

# Pressure:
SanJoaquin$pressureMean <- rollapply(SanJoaquin$pressure_set_1d, width=500,
                                       FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                       by.column=TRUE, partial=TRUE, fill=NA, align="center")

SanJoaquin$pressure_set_1d <- as.numeric(ifelse(is.na(SanJoaquin$pressure_set_1d), 
                                                (SanJoaquin$pressureMean),
                                                (SanJoaquin$pressure_set_1d)))

SanJoaquin$pressure_millibar <- c(SanJoaquin$pressure_set_1d * 0.01)

### Stream Metabolizer estimates ###
SanJoaquin$depth <- calc_depth(Q=u(SanJoaquin$discharge, "m^3 s^-1"), f=u(0.36))

latitude <- c(37.3472151)
longitude <- c(-120.9761777) 

SanJoaquin$solar.time <- calc_solar_time(SanJoaquin$timestamp, longitude)
SanJoaquin$light <- calc_light(SanJoaquin$solar.time,
                               latitude,
                               longitude,
                               max.PAR = u(2326, "umol m^-2 s^-1"),
                               attach.units = is.unitted(SanJoaquin$solar.time))

SanJoaquin$Do.sat <- calc_DO_sat(SanJoaquin$tempC, 
                                 SanJoaquin$pressure_millibar, sal=0) 
names(SanJoaquin)
# colnames(SanJoaquin)[7] <- "temp.water"
# colnames(SanJoaquin)[9] <- "DO.obs"

SanJoaquinQ <- subset(SanJoaquin, select= c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light", "discharge"))

# write.csv(SanJoaquinQ, paste0(outputDir, "/11273400SanJoaquinDat_Infill.csv")) 



## ---------------------------
# SACRAMENTO R A FREEPORT CA
## ---------------------------
dat4 <- read.csv("https://www.dropbox.com/s/lke4l9r91ktculf/nwis.11447650.csv?dl=1")
dat4$timestampPDT <- as.POSIXct(dat4$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
dat4$timestamp <- as.POSIXct(round_date(dat4$timestampPDT, hour,unit="5 minutes"))
summary(dat4)

#Mesowest data near (38.4556, -121.50181)
KSAC <- read.csv("https://www.dropbox.com/s/pw3jvbvwvmac0sx/KSAC.2020-10-01.csv?dl=1")
summary(KSAC)

##
# Fix timestamp (round to nearest 5 minute):
KSAC$timestampUTC <- as.POSIXct(KSAC$Date_Time, format = "%Y-%m-%dT%H:%M", tz="UTC")
KSAC$timestamp <- as.POSIXct(round_date(KSAC$timestampUTC, hour,unit="5 minutes"))

names(KSAC)
# Merge data files 
SacramentoFP <- left_join(dat4, KSAC[c("timestamp", "air_temp_set_1", "pressure_set_1d")],
                          by = c("timestamp" = "timestamp"))

SacramentoFP <- subset(SacramentoFP, timestamp >= as.POSIXct('2015-10-01 00:53:00') & 
                         timestamp <= as.POSIXct('2020-10-01 00:00:00'))
range(SacramentoFP$timestamp)
summary(SacramentoFP)
## 
## Infill missing data: 
##   

# Water Temp
SacramentoFP$tempCMean <- rollapply(SacramentoFP$tempC, width=500,
                                  FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                  by.column=TRUE, partial=TRUE, fill=NA, align="center")

SacramentoFP$tempC<- as.numeric(ifelse(is.na(SacramentoFP$tempC),
                                     (SacramentoFP$tempCMean),
                                     (SacramentoFP$tempC)))

# Discharge:
SacramentoFP$dischargeMean <- rollapply(SacramentoFP$discharge, width=500,
                                      FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                      by.column=TRUE, partial=TRUE, fill=NA, align="center")

SacramentoFP$discharge <- as.numeric(ifelse(is.na(SacramentoFP$discharge), 
                                          (SacramentoFP$dischargeMean),
                                          (SacramentoFP$discharge)))
# Dissolved oxygen:
SacramentoFP$DOmgLMean <- rollapply(SacramentoFP$DOmgL, width=500,
                                  FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                  by.column=TRUE, partial=TRUE, fill=NA, align="center")

SacramentoFP$DOmgL<- as.numeric(ifelse(is.na(SacramentoFP$DOmgL), 
                                     (SacramentoFP$DOmgLMean),
                                     (SacramentoFP$DOmgL)))

# Pressure:
SacramentoFP$pressureMean <- rollapply(SacramentoFP$pressure_set_1d, width=500,
                                     FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                     by.column=TRUE, partial=TRUE, fill=NA, align="center")

SacramentoFP$pressure_set_1d <- as.numeric(ifelse(is.na(SacramentoFP$pressure_set_1d), 
                                                (SacramentoFP$pressureMean),
                                                (SacramentoFP$pressure_set_1d)))

SacramentoFP$pressure_millibar <- c(SacramentoFP$pressure_set_1d * 0.01)
summary(SacramentoFP)

### Stream Metabolizer estimates ###
SacramentoFP$depth <- calc_depth(Q=u(SacramentoFP$discharge, "m^3 s^-1"), f=u(0.36))

latitude <- c(38.4556)
longitude <- c(-121.50181) 

SacramentoFP$solar.time <- calc_solar_time(SacramentoFP$timestamp, longitude)
SacramentoFP$light <- calc_light(SacramentoFP$solar.time,
                                 latitude,
                                 longitude,
                                 max.PAR = u(2326, "umol m^-2 s^-1"),
                                 attach.units = is.unitted(SacramentoFP$solar.time))

SacramentoFP$DO.sat <- calc_DO_sat(SacramentoFP$tempC, 
                                   SacramentoFP$pressure_millibar, sal=0) 
names(SacramentoFP)
# colnames(SacramentoFP)[5] <- "temp.water"
# colnames(SacramentoFP)[11] <- "DO.obs"

SacramentoFPQ <- subset(SacramentoFP, select= c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light", "discharge"))

# write.csv(SacramentoFPQ, paste0(outputDir, "/11447650SacramentoFPDat_Infill.csv")) 



## ---------------------------
# SACRAMENTO R AB DELTA CROSS CHANNEL CA
## ---------------------------
dat5 <- read.csv("https://www.dropbox.com/s/us0cl9yygoo81ft/nwis.11447890.csv?dl=1")
dat5$timestampPDT <- as.POSIXct(dat5$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
dat5$timestamp <- as.POSIXct(round_date(dat5$timestampPDT, hour,unit="5 minutes"))
summary(dat5)


#Mesowest data near (38.25769218, -121.5182865)
KSAC <- read.csv("https://www.dropbox.com/s/pw3jvbvwvmac0sx/KSAC.2020-10-01.csv?dl=1")
summary(KSAC)

##
# Fix timestamp (round to nearest 5 minute):
KSAC$timestampUTC <- as.POSIXct(KSAC$Date_Time, format = "%Y-%m-%dT%H:%M", tz="UTC")
KSAC$timestamp <- as.POSIXct(round_date(KSAC$timestampUTC, hour,unit="5 minutes"))

names(KSAC)
# Merge data files 
SacramentoDCC <- left_join(dat5, KSAC[c("timestamp", "air_temp_set_1", "pressure_set_1d")],
                           by = c("timestamp" = "timestamp"))
summary(SacramentoDCC)

SacramentoDCC <- subset(SacramentoDCC, timestamp >= as.POSIXct('2015-10-01 00:53:00') & 
                          timestamp <= as.POSIXct('2020-10-01 00:00:00'))
range(SacramentoDCC$timestamp)

## 
## Infill missing data: 
##   

# Water Temp
SacramentoDCC$tempCMean <- rollapply(SacramentoDCC$tempC, width=500,
                                    FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                    by.column=TRUE, partial=TRUE, fill=NA, align="center")

SacramentoDCC$tempC<- as.numeric(ifelse(is.na(SacramentoDCC$tempC),
                                       (SacramentoDCC$tempCMean),
                                       (SacramentoDCC$tempC)))

# Discharge:
SacramentoDCC$dischargeMean <- rollapply(SacramentoDCC$discharge, width=500,
                                        FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                        by.column=TRUE, partial=TRUE, fill=NA, align="center")

SacramentoDCC$discharge <- as.numeric(ifelse(is.na(SacramentoDCC$discharge), 
                                            (SacramentoDCC$dischargeMean),
                                            (SacramentoDCC$discharge)))
# Dissolved oxygen:
SacramentoDCC$DOmgLMean <- rollapply(SacramentoDCC$DOmgL, width=500,
                                    FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                    by.column=TRUE, partial=TRUE, fill=NA, align="center")

SacramentoDCC$DOmgL<- as.numeric(ifelse(is.na(SacramentoDCC$DOmgL), 
                                       (SacramentoDCC$DOmgLMean),
                                       (SacramentoDCC$DOmgL)))

# Pressure:
SacramentoDCC$pressureMean <- rollapply(SacramentoDCC$pressure_set_1d, width=500,
                                       FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                       by.column=TRUE, partial=TRUE, fill=NA, align="center")

SacramentoDCC$pressure_set_1d <- as.numeric(ifelse(is.na(SacramentoDCC$pressure_set_1d), 
                                                  (SacramentoDCC$pressureMean),
                                                  (SacramentoDCC$pressure_set_1d)))

SacramentoDCC$pressure_millibar <- c(SacramentoDCC$pressure_set_1d * 0.01)
summary(SacramentoDCC)


### Stream Metabolizer estimates ###
SacramentoDCC$depth <- calc_depth(Q=u(SacramentoDCC$discharge, "m^3 s^-1"), f=u(0.36))

latitude <- c(38.25769218)
longitude <- c(-121.5182865) 

SacramentoDCC$solar.time <- calc_solar_time(SacramentoDCC$timestamp, longitude)
SacramentoDCC$light <- calc_light(SacramentoDCC$solar.time,
                                  latitude,
                                  longitude,
                                  max.PAR = u(2326, "umol m^-2 s^-1"),
                                  attach.units = is.unitted(SacramentoDCC$solar.time))

SacramentoDCC$DO.sat <- calc_DO_sat(SacramentoDCC$tempC, 
                                    SacramentoDCC$pressure_millibar, sal=0) 
names(SacramentoDCC)
# colnames(SacramentoDCC)[7] <- "temp.water"
# colnames(SacramentoDCC)[9] <- "DO.obs"


SacramentoDCCQ <- subset(SacramentoDCC, select= c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light", "discharge"))

# write.csv(SacramentoDCCQ, paste0(outputDir, "/11447890SacramentoDCCDat_Raw.csv")) 




## ---------------------------
# RUSSIAN R NR HOPLAND CA
## ---------------------------
dat6 <- read.csv("https://www.dropbox.com/s/6kq7ya5clbpzuum/nwis.11462500.csv?dl=1")
dat6$timestampPDT <- as.POSIXct(dat6$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
dat6$timestamp <- as.POSIXct(round_date(dat6$timestampPDT, hour,unit="5 minutes"))

#Mesowest data near (39.02656314, -123.1305588)
KUKI <- read.csv("https://www.dropbox.com/s/5oc6xcpzs75jy9i/KUKI.2020-10-01.csv?dl=1")
summary(KUKI)
##
# Fix timestamp (round to nearest 5 minute):
KUKI$timestampUTC <- as.POSIXct(KUKI$Date_Time, format = "%Y-%m-%dT%H:%M", tz="UTC")
KUKI$timestamp <- as.POSIXct(round_date(KUKI$timestampUTC, hour,unit="5 minutes"))

names(KUKI)
# Merge data files 
RussianH <- left_join(dat6, KUKI[c("timestamp", "air_temp_set_1", "pressure_set_1d")],
                      by = c("timestamp" = "timestamp"))

RussianH <- subset(RussianH, timestamp >= as.POSIXct('2016-10-01 00:53:00') & 
                     timestamp <= as.POSIXct('2020-10-01 00:00:00'))
range(RussianH$timestamp)
summary(RussianH)

## 
## Infill missing data: 
##   


# Infill some data with rolling averages:
ggplot(RussianH, aes(x=timestampPDT, y=(discharge))) +
  geom_point(alpha = 0.7)  +
  theme_classic()

# Water Temp
RussianH$tempCMean <- rollapply(RussianH$tempC, width=500,
                                     FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                     by.column=TRUE, partial=TRUE, fill=NA, align="center")

RussianH$tempC<- as.numeric(ifelse(is.na(RussianH$tempC),
                                        (RussianH$tempCMean),
                                        (RussianH$tempC)))

# Discharge:
RussianH$dischargeMean <- rollapply(RussianH$discharge, width=500,
                                         FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                         by.column=TRUE, partial=TRUE, fill=NA, align="center")

RussianH$discharge <- as.numeric(ifelse(is.na(RussianH$discharge), 
                                             (RussianH$dischargeMean),
                                             (RussianH$discharge)))
# Dissolved oxygen:
RussianH$DOmgLMean <- rollapply(RussianH$DOmgL, width=500,
                                     FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                     by.column=TRUE, partial=TRUE, fill=NA, align="center")

RussianH$DOmgL<- as.numeric(ifelse(is.na(RussianH$DOmgL), 
                                        (RussianH$DOmgLMean),
                                        (RussianH$DOmgL)))

# Pressure:
RussianH$pressureMean <- rollapply(RussianH$pressure_set_1d, width=500,
                                        FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                        by.column=TRUE, partial=TRUE, fill=NA, align="center")

RussianH$pressure_set_1d <- as.numeric(ifelse(is.na(RussianH$pressure_set_1d), 
                                                   (RussianH$pressureMean),
                                                   (RussianH$pressure_set_1d)))

RussianH$pressure_millibar <- c(RussianH$pressure_set_1d * 0.01)

summary(RussianH)


### Stream Metabolizer estimates ###
RussianH$depth <- calc_depth(Q=u(RussianH$discharge, "m^3 s^-1"), f=u(0.36))

latitude <- c(39.02656314)
longitude <- c(-123.1305588) 

RussianH$solar.time <- calc_solar_time(RussianH$timestamp, longitude)
RussianH$light <- calc_light(RussianH$solar.time,
                             latitude,
                             longitude,
                             max.PAR = u(2326, "umol m^-2 s^-1"),
                             attach.units = is.unitted(RussianH$solar.time))

RussianH$DO.sat <- calc_DO_sat(RussianH$tempC, 
                               RussianH$pressure_millibar, sal=0) 
names(RussianH)
colnames(RussianH)[5] <- "temp.water"
colnames(RussianH)[9] <- "DO.obs"

RussianHQ <- subset(RussianH, select= c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light", "discharge"))

# write.csv(RussianHQ, paste0(outputDir, "/11462500RussianHQDat_Infill.csv")) 



## ---------------------------
# RUSSIAN R A DIGGER BEND NR HEALDSBURG CA
## ---------------------------
dat7 <- read.csv("https://www.dropbox.com/s/gjr11v5t6pxjn2k/nwis.11463980.csv?dl=1")
dat7$timestampPDT <- as.POSIXct(dat7$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
dat7$timestamp <- as.POSIXct(round_date(dat7$timestampPDT, hour,unit="5 minutes"))
summary(dat7)

#Mesowest data near (38.6329653, -122.8555494)
KO69 <- read_csv("https://www.dropbox.com/s/zu85xhoqrdqd9rm/KO69.2020-10-01.csv?dl=1")
summary(KO69)

##
# Fix timestamp (round to nearest 5 minute):
KO69$timestampUTC <- as.POSIXct(KO69$Date_Time, format = "%Y-%m-%dT%H:%M", tz="UTC")
KO69$timestamp <- as.POSIXct(round_date(KO69$timestampUTC, hour,unit="5 minutes"))

names(KO69)
# Merge data files 
RussianDB <- left_join(dat7, KO69[c("timestamp", "air_temp_set_1", "pressure_set_1d")],
                       by = c("timestamp" = "timestamp"))

RussianDB <- subset(RussianDB, timestamp >= as.POSIXct('2015-10-01 00:53:00') & 
                      timestamp <= as.POSIXct('2020-10-01 00:00:00'))
range(RussianDB$timestamp)
summary(RussianDB)
## 
## Infill missing data: 
##   


# Infill some data with rolling averages:
ggplot(RussianDB, aes(x=timestampPDT, y=(discharge))) +
  geom_point(alpha = 0.7)  +
  theme_classic()

# Water Temp
RussianDB$tempCMean <- rollapply(RussianDB$tempC, width=500,
                                FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                by.column=TRUE, partial=TRUE, fill=NA, align="center")

RussianDB$tempC<- as.numeric(ifelse(is.na(RussianDB$tempC),
                                   (RussianDB$tempCMean),
                                   (RussianDB$tempC)))

# Discharge:
RussianDB$dischargeMean <- rollapply(RussianDB$discharge, width=500,
                                    FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                    by.column=TRUE, partial=TRUE, fill=NA, align="center")

RussianDB$discharge <- as.numeric(ifelse(is.na(RussianDB$discharge), 
                                        (RussianDB$dischargeMean),
                                        (RussianDB$discharge)))
# Dissolved oxygen:
RussianDB$DOmgLMean <- rollapply(RussianDB$DOmgL, width=500,
                                FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                by.column=TRUE, partial=TRUE, fill=NA, align="center")

RussianDB$DOmgL<- as.numeric(ifelse(is.na(RussianDB$DOmgL), 
                                   (RussianDB$DOmgLMean),
                                   (RussianDB$DOmgL)))

# Pressure:
RussianDB$pressureMean <- rollapply(RussianDB$pressure_set_1d, width=500,
                                   FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                   by.column=TRUE, partial=TRUE, fill=NA, align="center")

RussianDB$pressure_set_1d <- as.numeric(ifelse(is.na(RussianDB$pressure_set_1d), 
                                              (RussianDB$pressureMean),
                                              (RussianDB$pressure_set_1d)))

RussianDB$pressure_millibar <- c(RussianDB$pressure_set_1d * 0.01)

summary(RussianDB)

### Stream Metabolizer estimates ###
RussianDB$depth <- calc_depth(Q=u(RussianDB$discharge, "m^3 s^-1"), f=u(0.36))

latitude <- c(38.6329653)
longitude <- c(-122.8555494) 

RussianDB$solar.time <- calc_solar_time(RussianDB$timestamp, longitude)
RussianDB$light <- calc_light(RussianDB$solar.time,
                              latitude,
                              longitude,
                              max.PAR = u(2326, "umol m^-2 s^-1"),
                              attach.units = is.unitted(RussianDB$solar.time))

RussianDB$DO.sat <- calc_DO_sat(RussianDB$tempC, 
                                RussianDB$pressure_millibar, sal=0) 
names(RussianDB)
# colnames(RussianDB)[7] <- "temp.water"
# colnames(RussianDB)[9] <- "DO.obs"

RussianDBQ <- subset(RussianDB, select= c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light", "discharge"))

# write.csv(RussianDBQ, paste0(outputDir, "/11463980RussianDBDat_Infill.csv")) 



## ---------------------------
# DRY C BLW LAMBERT BR NR GEYSERVILLE CA
## ---------------------------
dat8 <- read.csv("https://www.dropbox.com/s/1zm8kzx909m6jh5/nwis.11465420.csv?dl=1")
dat8$timestampPDT <- as.POSIXct(dat8$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
dat8$timestamp <- as.POSIXct(round_date(dat8$timestampPDT, hour,unit="5 minutes"))
summary(dat8)

#Mesowest data near (38.65324355, -122.9272191)
KO69 <- read_csv("https://www.dropbox.com/s/zu85xhoqrdqd9rm/KO69.2020-10-01.csv?dl=1")
summary(KO69)

##
# Fix timestamp (round to nearest 5 minute):
KO69$timestampUTC <- as.POSIXct(KO69$Date_Time, format = "%Y-%m-%dT%H:%M", tz="UTC")
KO69$timestamp <- as.POSIXct(round_date(KO69$timestampUTC, hour,unit="5 minutes"))

names(KO69)
# Merge data files 
DryCBLW <- left_join(dat8, KO69[c("timestamp", "air_temp_set_1", "pressure_set_1d")],
                     by = c("timestamp" = "timestamp"))
summary(DryCBLW)

DryCBLW <- subset(DryCBLW, timestamp >= as.POSIXct('2015-10-01 00:53:00') & 
                    timestamp <= as.POSIXct('2020-10-01 00:00:00'))
range(DryCBLW$timestampPDT)
summary(DryCBLW)

## 
## Infill missing data: 
##   

# Infill some data with rolling averages:
ggplot(DryCBLW, aes(x=timestampPDT, y=(discharge))) +
  geom_point(alpha = 0.7)  +
  theme_classic()

# Water Temp
DryCBLW$tempCMean <- rollapply(DryCBLW$tempC, width=500,
                                 FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                 by.column=TRUE, partial=TRUE, fill=NA, align="center")

DryCBLW$tempC<- as.numeric(ifelse(is.na(DryCBLW$tempC),
                                    (DryCBLW$tempCMean),
                                    (DryCBLW$tempC)))

# Discharge:
DryCBLW$dischargeMean <- rollapply(DryCBLW$discharge, width=500,
                                     FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                     by.column=TRUE, partial=TRUE, fill=NA, align="center")

DryCBLW$discharge <- as.numeric(ifelse(is.na(DryCBLW$discharge), 
                                         (DryCBLW$dischargeMean),
                                         (DryCBLW$discharge)))
# Dissolved oxygen:
DryCBLW$DOmgLMean <- rollapply(DryCBLW$DOmgL, width=500,
                                 FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                 by.column=TRUE, partial=TRUE, fill=NA, align="center")

DryCBLW$DOmgL<- as.numeric(ifelse(is.na(DryCBLW$DOmgL), 
                                    (DryCBLW$DOmgLMean),
                                    (DryCBLW$DOmgL)))

# Pressure:
DryCBLW$pressureMean <- rollapply(DryCBLW$pressure_set_1d, width=500,
                                    FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                    by.column=TRUE, partial=TRUE, fill=NA, align="center")

DryCBLW$pressure_set_1d <- as.numeric(ifelse(is.na(DryCBLW$pressure_set_1d), 
                                               (DryCBLW$pressureMean),
                                               (DryCBLW$pressure_set_1d)))

DryCBLW$pressure_millibar <- c(DryCBLW$pressure_set_1d * 0.01)
summary(DryCBLW)

### Stream Metabolizer estimates ###
DryCBLW$depth <- calc_depth(Q=u(DryCBLW$discharge, "m^3 s^-1"), f=u(0.36))

latitude <- c(38.65324355)
longitude <- c(-122.9272191) 

DryCBLW$solar.time <- calc_solar_time(DryCBLW$timestamp, longitude)
DryCBLW$light <- calc_light(DryCBLW$solar.time,
                            latitude,
                            longitude,
                            max.PAR = u(2326, "umol m^-2 s^-1"),
                            attach.units = is.unitted(DryCBLW$solar.time))

DryCBLW$DO.sat <- calc_DO_sat(DryCBLW$tempC, 
                              DryCBLW$pressure_millibar, sal=0) 
names(DryCBLW)
# colnames(DryCBLW)[7] <- "temp.water"
# colnames(DryCBLW)[9] <- "DO.obs"

DryCBLWQ <- subset(DryCBLW, select= c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light", "discharge"))

# write.csv(DryCBLWQ, paste0(outputDir, "/11465240DryCBLWDat_Infill.csv")) 



## ---------------------------
# RUSSIAN R NR GUERNEVILLE CA
## ---------------------------
dat9 <- read.csv("https://www.dropbox.com/s/ich3wz5hy6bbasl/nwis.11467000.csv?dl=1")
dat9$timestampPDT <- as.POSIXct(dat9$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
dat9$timestamp <- as.POSIXct(round_date(dat9$timestampPDT, hour,unit="5 minutes"))
summary(dat9)

#Mesowest data near (38.50852336, -122.9277737)
KO69 <- read_csv("https://www.dropbox.com/s/zu85xhoqrdqd9rm/KO69.2020-10-01.csv?dl=1")
summary(KO69)

##
# Fix timestamp (round to nearest 5 minute):
KO69$timestampUTC <- as.POSIXct(KO69$Date_Time, format = "%Y-%m-%dT%H:%M", tz="UTC")
KO69$timestamp <- as.POSIXct(round_date(KO69$timestampUTC, hour,unit="5 minutes"))

names(KO69)
# Merge data files 
RussianNRG <- left_join(dat9, KO69[c("timestamp", "air_temp_set_1", "pressure_set_1d")],
                        by = c("timestamp" = "timestamp"))
summary(RussianNRG)

RussianNRG <- subset(RussianNRG, timestamp >= as.POSIXct('2015-10-01 00:53:00') & 
                       timestamp <= as.POSIXct('2020-10-01 00:00:00'))
range(RussianNRG$timestamp)
summary(RussianNRG)

## 
## Infill missing data: 
##   

# Infill some data with rolling averages:
ggplot(RussianNRG, aes(x=timestampPDT, y=(discharge))) +
  geom_point(alpha = 0.7)  +
  theme_classic()

# Water Temp
RussianNRG$tempCMean <- rollapply(RussianNRG$tempC, width=500,
                               FUN=function(x) mean(x, na.rm=TRUE), by=1,
                               by.column=TRUE, partial=TRUE, fill=NA, align="center")

RussianNRG$tempC<- as.numeric(ifelse(is.na(RussianNRG$tempC),
                                  (RussianNRG$tempCMean),
                                  (RussianNRG$tempC)))

# Discharge:
RussianNRG$dischargeMean <- rollapply(RussianNRG$discharge, width=500,
                                   FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                   by.column=TRUE, partial=TRUE, fill=NA, align="center")

RussianNRG$discharge <- as.numeric(ifelse(is.na(RussianNRG$discharge), 
                                       (RussianNRG$dischargeMean),
                                       (RussianNRG$discharge)))
# Dissolved oxygen:
RussianNRG$DOmgLMean <- rollapply(RussianNRG$DOmgL, width=500,
                               FUN=function(x) mean(x, na.rm=TRUE), by=1,
                               by.column=TRUE, partial=TRUE, fill=NA, align="center")

RussianNRG$DOmgL<- as.numeric(ifelse(is.na(RussianNRG$DOmgL), 
                                  (RussianNRG$DOmgLMean),
                                  (RussianNRG$DOmgL)))

# Pressure:
RussianNRG$pressureMean <- rollapply(RussianNRG$pressure_set_1d, width=500,
                                  FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                  by.column=TRUE, partial=TRUE, fill=NA, align="center")

RussianNRG$pressure_set_1d <- as.numeric(ifelse(is.na(RussianNRG$pressure_set_1d), 
                                             (RussianNRG$pressureMean),
                                             (RussianNRG$pressure_set_1d)))

RussianNRG$pressure_millibar <- c(RussianNRG$pressure_set_1d * 0.01)
summary(RussianNRG)

### Stream Metabolizer estimates ###
RussianNRG$depth <- calc_depth(Q=u(RussianNRG$discharge, "m^3 s^-1"), f=u(0.36))

latitude <- c(38.50852336)
longitude <- c(-122.9277737) 

RussianNRG$solar.time <- calc_solar_time(RussianNRG$timestamp, longitude)
RussianNRG$light <- calc_light(RussianNRG$solar.time,
                               latitude,
                               longitude,
                               max.PAR = u(2326, "umol m^-2 s^-1"),
                               attach.units = is.unitted(RussianNRG$solar.time))

RussianNRG$DO.sat <- calc_DO_sat(RussianNRG$tempC, 
                                 RussianNRG$pressure_millibar, sal=0) 
names(RussianNRG)
summary(RussianNRG)
# colnames(RussianNRG)[5] <- "temp.water"
# colnames(RussianNRG)[9] <- "DO.obs"

RussianNRGQ <- subset(RussianNRG, select= c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light", "discharge"))

# write.csv(RussianNRGQ, paste0(outputDir, "/11467000RussianNRGDat_Infill.csv")) 



## ---------------------------
# WHITE RIVER AT R STREET NEAR AUBURN, WA 
## ---------------------------
dat10 <- read.csv("https://www.dropbox.com/s/mk87oydbu601zqo/nwis.12100490.csv?dl=1")
dat10$timestampPDT <- as.POSIXct(dat10$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
dat10$timestamp <- as.POSIXct(round_date(dat10$timestampPDT, hour,unit="5 minutes"))
summary(dat10)

#Mesowest data near (47.27482295, -122.2078967)
KRNT <- read.csv("https://www.dropbox.com/s/amdmo1obr0upo46/KRNT.2020-10-01.csv?dl=1")
summary(KRNT)
##
# Fix timestamp (round to nearest 5 minute):
KRNT$timestampUTC <- as.POSIXct(KRNT$Date_Time, format = "%Y-%m-%dT%H:%M", tz="UTC")
KRNT$timestamp <- as.POSIXct(round_date(KRNT$timestampUTC, hour,unit="5 minutes"))


names(KRNT)
# Merge data files 
WhiteA <- left_join(dat10, KRNT[c("timestamp", "air_temp_set_1", "pressure_set_1d")],
                    by = c("timestamp" = "timestamp"))

WhiteA <- subset(WhiteA, timestamp >= as.POSIXct('2015-10-01 00:53:00') & 
                   timestamp <= as.POSIXct('2020-10-01 00:00:00'))
range(WhiteA$timestamp)
summary(WhiteA)
## 
## Infill missing data: 
##   

# Infill some data with rolling averages:
ggplot(WhiteA, aes(x=timestampPDT, y=(tempC))) +
  geom_point(alpha = 0.7)  +
  theme_classic()

# Water Temp
WhiteA$tempCMean <- rollapply(WhiteA$tempC, width=500,
                                  FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                  by.column=TRUE, partial=TRUE, fill=NA, align="center")

WhiteA$tempC<- as.numeric(ifelse(is.na(WhiteA$tempC),
                                     (WhiteA$tempCMean),
                                     (WhiteA$tempC)))

# Discharge:
WhiteA$dischargeMean <- rollapply(WhiteA$discharge, width=500,
                                      FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                      by.column=TRUE, partial=TRUE, fill=NA, align="center")

WhiteA$discharge <- as.numeric(ifelse(is.na(WhiteA$discharge), 
                                          (WhiteA$dischargeMean),
                                          (WhiteA$discharge)))
# Dissolved oxygen:
WhiteA$DOmgLMean <- rollapply(WhiteA$DOmgL, width=500,
                                  FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                  by.column=TRUE, partial=TRUE, fill=NA, align="center")

WhiteA$DOmgL<- as.numeric(ifelse(is.na(WhiteA$DOmgL), 
                                     (WhiteA$DOmgLMean),
                                     (WhiteA$DOmgL)))

# Pressure:
WhiteA$pressureMean <- rollapply(WhiteA$pressure_set_1d, width=500,
                                     FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                     by.column=TRUE, partial=TRUE, fill=NA, align="center")

WhiteA$pressure_set_1d <- as.numeric(ifelse(is.na(WhiteA$pressure_set_1d), 
                                                (WhiteA$pressureMean),
                                                (WhiteA$pressure_set_1d)))

WhiteA$pressure_millibar <- c(WhiteA$pressure_set_1d * 0.01)
summary(WhiteA)

### Stream Metabolizer estimates ###
WhiteA$depth <- calc_depth(Q=u(WhiteA$discharge, "m^3 s^-1"), f=u(0.36))

latitude <- c(47.27482295)
longitude <- c(-122.2078967) 

WhiteA$solar.time <- calc_solar_time(WhiteA$timestamp, longitude)
WhiteA$light <- calc_light(WhiteA$solar.time,
                           latitude,
                           longitude,
                           max.PAR = u(2326, "umol m^-2 s^-1"),
                           attach.units = is.unitted(WhiteA$solar.time))

WhiteA$DO.sat <- calc_DO_sat(WhiteA$tempC, 
                             WhiteA$pressure_millibar, sal=0) 
names(WhiteA)
# colnames(WhiteA)[7] <- "temp.water"
# colnames(WhiteA)[9] <- "DO.obs"

WhiteAQ <- subset(WhiteA, select= c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light", "discharge"))

# write.csv(WhiteAQ, paste0(outputDir, "/12100490WhiteADat_Infill.csv")) 



## ---------------------------
# CLARKS CREEK AT TACOMA ROAD NEAR PUYALLUP, WA 
## ---------------------------
dat11 <- read.csv("https://www.dropbox.com/s/s17a6c7hof5ukhs/nwis.12102075.csv?dl=1")
dat11$timestampPDT <- as.POSIXct(dat11$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
dat11$timestamp <- as.POSIXct(round_date(dat11$timestampPDT, hour,unit="5 minutes"))
summary(dat11)

#Mesowest data near (47.19760016, -122.337343)
KTCM <- read_csv("https://www.dropbox.com/s/r358tpw65u3qj5q/KTCM.2020-10-01.csv?dl=1")
summary(KTCM)
##
# Fix timestamp (round to nearest 5 minute):
KTCM$timestampUTC <- as.POSIXct(KTCM$Date_Time, format = "%Y-%m-%dT%H:%M", tz="UTC")
KTCM$timestamp <- as.POSIXct(round_date(KTCM$timestampUTC, hour,unit="5 minutes"))

names(KTCM)
# Merge data files 
ClarksC <- left_join(dat11, KTCM[c("timestamp", "air_temp_set_1", "pressure_set_1d")],
                     by = c("timestamp" = "timestamp"))

ClarksC <- subset(ClarksC, timestamp >= as.POSIXct('2015-10-01 00:53:00') & 
                    timestamp <= as.POSIXct('2020-10-01 00:00:00'))
range(ClarksC$timestamp)
summary(ClarksC)

## 
## Infill missing data: 
##   

# Infill some data with rolling averages:
ggplot(ClarksC, aes(x=timestampPDT, y=(tempC))) +
  geom_point(alpha = 0.7)  +
  theme_classic()

# Water Temp
ClarksC$tempCMean <- rollapply(ClarksC$tempC, width=500,
                              FUN=function(x) mean(x, na.rm=TRUE), by=1,
                              by.column=TRUE, partial=TRUE, fill=NA, align="center")

ClarksC$tempC<- as.numeric(ifelse(is.na(ClarksC$tempC),
                                 (ClarksC$tempCMean),
                                 (ClarksC$tempC)))

# Discharge:
ClarksC$dischargeMean <- rollapply(ClarksC$discharge, width=500,
                                  FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                  by.column=TRUE, partial=TRUE, fill=NA, align="center")

ClarksC$discharge <- as.numeric(ifelse(is.na(ClarksC$discharge), 
                                      (ClarksC$dischargeMean),
                                      (ClarksC$discharge)))
# Dissolved oxygen:
ClarksC$DOmgLMean <- rollapply(ClarksC$DOmgL, width=500,
                              FUN=function(x) mean(x, na.rm=TRUE), by=1,
                              by.column=TRUE, partial=TRUE, fill=NA, align="center")

ClarksC$DOmgL<- as.numeric(ifelse(is.na(ClarksC$DOmgL), 
                                 (ClarksC$DOmgLMean),
                                 (ClarksC$DOmgL)))
# Pressure:
ClarksC$pressureMean <- rollapply(ClarksC$pressure_set_1d, width=500,
                                 FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                 by.column=TRUE, partial=TRUE, fill=NA, align="center")

ClarksC$pressure_set_1d <- as.numeric(ifelse(is.na(ClarksC$pressure_set_1d), 
                                            (ClarksC$pressureMean),
                                            (ClarksC$pressure_set_1d)))

ClarksC$pressure_millibar <- c(ClarksC$pressure_set_1d * 0.01)
summary(ClarksC)

### Stream Metabolizer estimates ###
ClarksC$depth <- calc_depth(Q=u(ClarksC$discharge, "m^3 s^-1"), f=u(0.36))

latitude <- c(47.19760016)
longitude <- c(-122.337343) 

ClarksC$solar.time <- calc_solar_time(ClarksC$timestamp, longitude)
ClarksC$light <- calc_light(ClarksC$solar.time,
                            latitude,
                            longitude,
                            max.PAR = u(2326, "umol m^-2 s^-1"),
                            attach.units = is.unitted(ClarksC$solar.time))

ClarksC$DO.sat <- calc_DO_sat(ClarksC$tempC, 
                              ClarksC$pressure_millibar, sal=0) 
names(ClarksC)
# colnames(ClarksC)[7] <- "temp.water"
# colnames(ClarksC)[9] <- "DO.obs"

ClarksCQ <- subset(ClarksC, select= c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light", "discharge"))

# write.csv(ClarksCQ, paste0(outputDir, "/12102075ClarksCDat_Infill.csv")) 



## ---------------------------
# FANNO CREEK AT DURHAM, OR
## ---------------------------
dat13 <- read.csv("https://www.dropbox.com/s/vhteu1aox79mcr0/nwis.14206950.csv?dl=1")
dat13$timestampPDT <- as.POSIXct(dat13$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
dat13$timestamp <- as.POSIXct(round_date(dat13$timestampPDT, hour,unit="5 minutes"))
summary(dat13)

#Mesowest data near (45.403452, -122.7548185)
KHIO <- read_csv("https://www.dropbox.com/s/e9rcch76fb2z7ji/KHIO.2020-10-01.csv?dl=1")
##
# Fix timestamp (round to nearest 5 minute):
KHIO$timestampUTC <- as.POSIXct(KHIO$Date_Time, format = "%Y-%m-%dT%H:%M", tz="UTC")
KHIO$timestamp <- as.POSIXct(round_date(KHIO$timestampUTC, hour,unit="5 minutes"))

# Merge data files 
Fanno <- left_join(dat13, KHIO[c("timestamp", "air_temp_set_1", "pressure_set_1d")],
                   by = c("timestamp" = "timestamp"))

Fanno <- subset(Fanno, timestamp >= as.POSIXct('2015-10-01 00:53:00') & 
                  timestamp <= as.POSIXct('2020-10-01 00:00:00'))
range(Fanno$timestamp)
summary(Fanno)

## 
## Infill missing data: 
##   

# Infill some data with rolling averages:
ggplot(Fanno, aes(x=timestampPDT, y=(tempC))) +
  geom_point(alpha = 0.7)  +
  theme_classic()

# Water Temp
Fanno$tempCMean <- rollapply(Fanno$tempC, width=500,
                               FUN=function(x) mean(x, na.rm=TRUE), by=1,
                               by.column=TRUE, partial=TRUE, fill=NA, align="center")

Fanno$tempC<- as.numeric(ifelse(is.na(Fanno$tempC),
                                  (Fanno$tempCMean),
                                  (Fanno$tempC)))

# Discharge:
Fanno$discharge1 <- as.numeric(Fanno$discharge)
Fanno$dischargeMean <- rollapply(Fanno$discharge1, width=500,
                                   FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                   by.column=TRUE, partial=TRUE, fill=NA, align="center")

Fanno$discharge1 <- as.numeric(ifelse(is.na(Fanno$discharge1), 
                                       (Fanno$dischargeMean),
                                       (Fanno$discharge1)))
# Dissolved oxygen:
Fanno$DOmgLMean <- rollapply(Fanno$DOmgL, width=500,
                               FUN=function(x) mean(x, na.rm=TRUE), by=1,
                               by.column=TRUE, partial=TRUE, fill=NA, align="center")

Fanno$DOmgL<- as.numeric(ifelse(is.na(Fanno$DOmgL), 
                                  (Fanno$DOmgLMean),
                                  (Fanno$DOmgL)))

# Pressure:
Fanno$pressureMean <- rollapply(Fanno$pressure_set_1d, width=500,
                                  FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                  by.column=TRUE, partial=TRUE, fill=NA, align="center")

Fanno$pressure_set_1d <- as.numeric(ifelse(is.na(Fanno$pressure_set_1d), 
                                             (Fanno$pressureMean),
                                             (Fanno$pressure_set_1d)))
summary(Fanno)
Fanno$pressure_millibar <- c(Fanno$pressure_set_1d * 0.01)

### Stream Metabolizer estimates ###
Fanno$depth <- calc_depth(Q=u(Fanno$discharge1, "m^3 s^-1"), f=u(0.36))

latitude <- c(45.40352)
longitude <- c(-122.75481) 

Fanno$solar.time <- calc_solar_time(Fanno$timestamp, longitude)
Fanno$light <- calc_light(Fanno$solar.time,
                          latitude,
                          longitude,
                          max.PAR = u(2326, "umol m^-2 s^-1"),
                          attach.units = is.unitted(Fanno$solar.time))

Fanno$DO.sat <- calc_DO_sat(Fanno$tempC, 
                            Fanno$pressure_millibar, sal=0) 
names(Fanno)
# colnames(Fanno)[7] <- "temp.water"
# colnames(Fanno)[9] <- "DO.obs"

FannoQ <- subset(Fanno, select= c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light", "discharge"))

# write.csv(FannoQ, paste0(outputDir, "/14206950FannoDat_Infill.csv")) 


## ---------------------------
# MCKENZIE RIVER NEAR VIDA, OR
## ---------------------------
dat12 <- read.csv("https://www.dropbox.com/s/oq5h9m1mb6pcr5y/nwis.14162500.csv?dl=1")
dat12$timestampPDT <- as.POSIXct(dat12$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
dat12$timestamp <- as.POSIXct(round_date(dat12$timestampPDT, hour,unit="5 minutes"))
summary(dat12)

#Mesowest data near (44.1248487, -122.47062)
K6S2 <- read_csv("https://www.dropbox.com/s/qcax8ll7758ahw8/K6S2.2020-10-01.csv?dl=1")
summary(K6S2)
##
# Fix timestamp (round to nearest 5 minute):
K6S2$timestampUTC <- as.POSIXct(K6S2$Date_Time, format = "%Y-%m-%dT%H:%M", tz="UTC")
K6S2$timestamp <- as.POSIXct(round_date(K6S2$timestampUTC, hour,unit="5 minutes"))

names(K6S2)
# Merge data files 
Mckenzie <- left_join(dat12, K6S2[c("timestamp", "air_temp_set_1", "pressure_set_1d")],
                      by = c("timestamp" = "timestamp"))

Mckenzie <- subset(Mckenzie, timestamp >= as.POSIXct('2015-10-01 00:53:00') & 
                     timestamp<= as.POSIXct('2020-10-01 00:00:00'))
range(Mckenzie$timestamp)
summary(Mckenzie)

## 
## Infill missing data: 
##   *** Missing 2019 and 2020 pressure data

# Infill some data with rolling averages:
ggplot(Mckenzie, aes(x=timestampPDT, y=(pressure_set_1d))) +
  geom_point(alpha = 0.7)  +
  theme_classic()

# Water Temp
Mckenzie$tempCMean <- rollapply(Mckenzie$tempC, width=500,
                             FUN=function(x) mean(x, na.rm=TRUE), by=1,
                             by.column=TRUE, partial=TRUE, fill=NA, align="center")

Mckenzie$tempC<- as.numeric(ifelse(is.na(Mckenzie$tempC),
                                (Mckenzie$tempCMean),
                                (Mckenzie$tempC)))

# Discharge:
Mckenzie$discharge1 <- as.numeric(Mckenzie$discharge)
Mckenzie$dischargeMean <- rollapply(Mckenzie$discharge1, width=500,
                                 FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                 by.column=TRUE, partial=TRUE, fill=NA, align="center")

Mckenzie$discharge1 <- as.numeric(ifelse(is.na(Mckenzie$discharge1), 
                                      (Mckenzie$dischargeMean),
                                      (Mckenzie$discharge1)))
# Dissolved oxygen:
Mckenzie$DOmgLMean <- rollapply(Mckenzie$DOmgL, width=10000,
                             FUN=function(x) mean(x, na.rm=TRUE), by=1,
                             by.column=TRUE, partial=TRUE, fill=NA, align="center")

Mckenzie$DOmgL<- as.numeric(ifelse(is.na(Mckenzie$DOmgL), 
                                (Mckenzie$DOmgLMean),
                                (Mckenzie$DOmgL)))

# Pressure:
Mckenzie$pressureMean <- rollapply(Mckenzie$pressure_set_1d, width=10000,
                                FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                by.column=TRUE, partial=TRUE, fill=NA, align="center")

Mckenzie$pressure_set_1d <- as.numeric(ifelse(is.na(Mckenzie$pressure_set_1d), 
                                           (Mckenzie$pressureMean),
                                           (Mckenzie$pressure_set_1d)))

summary(Mckenzie)
Mckenzie$pressure_millibar <- c(Mckenzie$pressure_set_1d * 0.01)

### Stream Metabolizer estimates ###
Mckenzie$depth <- calc_depth(Q=u(Mckenzie$discharge1, "m^3 s^-1"), f=u(0.36))

latitude <- c(44.1248487)
longitude <- c(-122.47062) 

Mckenzie$solar.time <- calc_solar_time(Mckenzie$timestampPDT, longitude)
Mckenzie$light <- calc_light(Mckenzie$solar.time,
                      latitude,
                      longitude,
                      max.PAR = u(2326, "umol m^-2 s^-1"),
                      attach.units = is.unitted(Mckenzie$solar.time))

Mckenzie$DO.sat <- calc_DO_sat(Mckenzie$tempC, 
                               Mckenzie$pressure_millibar, sal=0) 
names(Mckenzie)
colnames(Mckenzie)[7] <- "discharge"
colnames(Mckenzie)[5] <- "temp.water"
colnames(Mckenzie)[9] <- "DO.obs"

MckenzieQ <- subset(Mckenzie, select= c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light", "discharge"))

# write.csv(MckenzieQ, paste0(outputDir, "/14162500MckenzieDat_Infill.csv")) 





## ---------------------------
# CLACKAMAS RIVER AT ESTACADA, OR
## ---------------------------
dat14 <- read.csv("https://www.dropbox.com/s/575866cahvfh3j7/nwis.14210000.csv?dl=1")
dat14$timestampPDT <- as.POSIXct(dat14$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
dat14$timestamp <- as.POSIXct(round_date(dat14$timestampPDT, hour,unit="5 minutes"))
summary(dat14)

#Mesowest data near (45.30053, -122.35359)
D7564 <- read_csv("https://www.dropbox.com/s/jqe4oouq2l2eqas/D7564.2020-10-01.csv?dl=1")
summary(D7564)
##
# Fix timestamp (round to nearest 5 minute):
D7564$timestampUTC <- as.POSIXct(D7564$Date_Time, format = "%Y-%m-%dT%H:%M", tz="UTC")
D7564$timestamp <- as.POSIXct(round_date(D7564$timestampUTC, hour,unit="5 minutes"))

names(D7564)
# Merge data files 
Clackamas <- left_join(dat14, D7564[c("timestamp", "air_temp_set_1", "pressure_set_1d")],
                       by = c("timestamp" = "timestamp"))

Clackamas <- subset(Clackamas, timestamp >= as.POSIXct('2015-10-01 00:53:00') & 
                      timestamp <= as.POSIXct('2020-10-01 00:00:00'))
range(Clackamas$timestamp)
summary(Clackamas)
## 
## Infill missing data: 
##   *** Missing 2019 and 2020 pressure data

# Infill some data with rolling averages:
ggplot(Clackamas, aes(x=timestampPDT, y=(tempC))) +
  geom_point(alpha = 0.7)  +
  theme_classic()

# Water Temp
Clackamas$tempCMean <- rollapply(Clackamas$tempC, width=500,
                                FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                by.column=TRUE, partial=TRUE, fill=NA, align="center")

Clackamas$tempC<- as.numeric(ifelse(is.na(Clackamas$tempC),
                                   (Clackamas$tempCMean),
                                   (Clackamas$tempC)))

# Discharge:
Clackamas$dischargeMean <- rollapply(Clackamas$discharge, width=500,
                                    FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                    by.column=TRUE, partial=TRUE, fill=NA, align="center")

Clackamas$discharge <- as.numeric(ifelse(is.na(Clackamas$discharge), 
                                         (Clackamas$dischargeMean),
                                         (Clackamas$discharge)))
# Dissolved oxygen:
Clackamas$DOmgLMean <- rollapply(Clackamas$DOmgL, width=500,
                                FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                by.column=TRUE, partial=TRUE, fill=NA, align="center")

Clackamas$DOmgL<- as.numeric(ifelse(is.na(Clackamas$DOmgL), 
                                   (Clackamas$DOmgLMean),
                                   (Clackamas$DOmgL)))

# Pressure:
Clackamas$pressureMean <- rollapply(Clackamas$pressure_set_1d, width=500,
                                   FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                   by.column=TRUE, partial=TRUE, fill=NA, align="center")

Clackamas$pressure_set_1d <- as.numeric(ifelse(is.na(Clackamas$pressure_set_1d), 
                                              (Clackamas$pressureMean),
                                              (Clackamas$pressure_set_1d)))

Clackamas$pressure_millibar <- c(Clackamas$pressure_set_1d * 0.01)
summary(Clackamas)

### Stream Metabolizer estimates ###
Clackamas$depth <- calc_depth(Q=u(Clackamas$discharge, "m^3 s^-1"), f=u(0.36))

latitude <- c(45.30053)
longitude <- c(-122.35359) 

Clackamas$solar.time <- calc_solar_time(Clackamas$timestamp, longitude)
Clackamas$light <- calc_light(Clackamas$solar.time,
                              latitude,
                              longitude,
                              max.PAR = u(2326, "umol m^-2 s^-1"),
                              attach.units = is.unitted(Clackamas$solar.time))

Clackamas$DO.sat <- calc_DO_sat(Clackamas$tempC, 
                                Clackamas$pressure_millibar, sal=0) 
names(Clackamas)
# colnames(Clackamas)[7] <- "temp.water"
# colnames(Clackamas)[9] <- "DO.obs"

ClackamasQ <- subset(Clackamas, select= c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light", "discharge"))

# write.csv(ClackamasQ, paste0(outputDir, "/14210000ClackamasDat_Infill.csv")) 



## ---------------------------
# CLACKAMAS RIVER NEAR OREGON CITY, OR 
## ---------------------------
dat15 <- read.csv("https://www.dropbox.com/s/mlv21xdjirevuaq/nwis.14211010.csv?dl=1")
dat15$timestampPDT <- as.POSIXct(dat15$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
dat15$timestamp <- as.POSIXct(round_date(dat15$timestampPDT, hour,unit="5 minutes"))
summary(dat15)

#Mesowest data near (45.3792874, -122.5773134)
ORCO3 <- read_csv("https://www.dropbox.com/s/tmiyxm8fl1gl55j/ORCO3.2020-10-01.csv?dl=1")
summary(ORCO3)
##
# Fix timestamp (round to nearest 5 minute):
ORCO3$timestampUTC <- as.POSIXct(ORCO3$Date_Time, format = "%Y-%m-%dT%H:%M", tz="UTC")
ORCO3$timestamp <- as.POSIXct(round_date(ORCO3$timestampUTC, hour,unit="5 minutes"))

names(ORCO3)
# Merge data files 
ClackamasOC <- left_join(dat15, ORCO3[c("timestamp", "air_temp_set_1", "altimeter_set_1d")],
                         by = c("timestamp" = "timestamp"))

ClackamasOC <- subset(ClackamasOC, timestamp >= as.POSIXct('2015-10-01 00:53:00') & 
                        timestamp <= as.POSIXct('2020-09-20 00:00:00'))
range(ClackamasOC$timestamp)
summary(ClackamasOC)

## 
## Infill missing data: 
##   *** Missing 2019 and 2020 pressure data

# Infill some data with rolling averages:
ggplot(ClackamasOC, aes(x=timestampPDT, y=(pressure_set_1d))) +
  geom_point(alpha = 0.7)  +
  theme_classic()

# Water Temp
ClackamasOC$tempCMean <- rollapply(ClackamasOC$tempC, width=500,
                                 FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                 by.column=TRUE, partial=TRUE, fill=NA, align="center")

ClackamasOC$tempC<- as.numeric(ifelse(is.na(ClackamasOC$tempC),
                                    (ClackamasOC$tempCMean),
                                    (ClackamasOC$tempC)))

# Discharge:
ClackamasOC$dischargeMean <- rollapply(ClackamasOC$discharge, width=500,
                                     FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                     by.column=TRUE, partial=TRUE, fill=NA, align="center")

ClackamasOC$discharge <- as.numeric(ifelse(is.na(ClackamasOC$discharge), 
                                         (ClackamasOC$dischargeMean),
                                         (ClackamasOC$discharge)))
# Dissolved oxygen:
ClackamasOC$DOmgLMean <- rollapply(ClackamasOC$DOmgL, width=500,
                                 FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                 by.column=TRUE, partial=TRUE, fill=NA, align="center")

ClackamasOC$DOmgL<- as.numeric(ifelse(is.na(ClackamasOC$DOmgL), 
                                    (ClackamasOC$DOmgLMean),
                                    (ClackamasOC$DOmgL)))

# Pressure:
ClackamasOC$pressureMean <- rollapply(ClackamasOC$pressure_set_1d, width=500,
                                    FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                    by.column=TRUE, partial=TRUE, fill=NA, align="center")

ClackamasOC$pressure_set_1d <- as.numeric(ifelse(is.na(ClackamasOC$pressure_set_1d), 
                                               (ClackamasOC$pressureMean),
                                               (ClackamasOC$pressure_set_1d)))

ClackamasOC$pressure_millibar <- c(ClackamasOC$pressure_set_1d * 0.01)
summary(ClackamasOC)


### Stream Metabolizer estimates ###
ClackamasOC$depth <- calc_depth(Q=u(ClackamasOC$discharge, "m^3 s^-1"), f=u(0.36))

latitude <- c(45.3792874)
longitude <- c(-122.5773134) 

ClackamasOC$solar.time <- calc_solar_time(ClackamasOC$timestamp, longitude)
ClackamasOC$light <- calc_light(ClackamasOC$solar.time,
                                latitude,
                                longitude,
                                max.PAR = u(2326, "umol m^-2 s^-1"),
                                attach.units = is.unitted(ClackamasOC$solar.time))

ClackamasOC$DO.sat <- calc_DO_sat(ClackamasOC$tempC, 
                                  ClackamasOC$pressure_millibar, sal=0) 
names(ClackamasOC)
# colnames(ClackamasOC)[7] <- "temp.water"
# colnames(ClackamasOC)[11] <- "DO.obs"

ClackamasOCQ <- subset(ClackamasOC, select= c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light", "discharge"))

#write.csv(ClackamasOC, paste0(outputDir, "/14211010ClackamasOCDat_Infill.csv")) 











































## ---------------------------
# Minor metabolism exploration... just to check each dataset. 
# mixed success. 
ClackamasQ %>% unitted::v() %>%
  mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
  select(solar.time, starts_with('DO')) %>%
  gather(type, DO.value, starts_with('DO')) %>%
  mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
  ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() +
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)', light='PAR\n(umol m^-2 s^-1)')
ClackamasQ %>% unitted::v() %>%
  select(solar.time, depth, temp.water, light) %>%
  gather(type, value, depth, temp.water, light) %>%
  mutate(
    type=ordered(type, levels=c('depth','temp.water','light')),
    units=ordered(labels[type], unname(labels))) %>%
  ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() +
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')


mm <- metab(specs(mm_name('mle')), data=ClackamasQ, info='my info')
predict_metab(mm)
get_info(mm)
get_fitting_time(mm)

plot_metab_preds(mm)
