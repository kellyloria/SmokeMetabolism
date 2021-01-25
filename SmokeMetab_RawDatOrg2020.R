#' Data Organization for modeling metabolism for signals of smoke effect based on time window of interest
#' 
#' @description Only for raw downloaded USGS data and Mesowest data
#'@details Last update: 2020-20-24
#' See http://usgs-r.github.io/streamMetabolizer for vignettes on the web.

## ---------------------------
# File path setup:
if (dir.exists('/Users/kellyloria/Documents/UNR_2020/Fall2020Projects/')){
  inputDir<- '/Users/kellyloria/Documents/UNR_2020/Fall2020Projects/'
  outputDir<- '/Users/kellyloria/Documents/UNR_2020/Fall2020Projects/SmokeMetabolism/SmokeMetabData/RawDat'
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
library(dataRetrieval)

## ---------------------------
# SANTA MARGARITA R NR TEMECULA
## ---------------------------

startDate <- "2015-10-01"
endDate <- "2020-10-01"

##
#select sites:
### 4
siteNumber <- "11044000" 
Info <- readNWISsite(siteNumber)
parameterCd <- c("00060",  # discharge
                 "00010",
                 "00300")  # Water temp

# readNWISqw

dat1 <- readNWISdv(siteNumber,parameterCd,
                      "1980-09-30","2020-10-01")

dat1 <- readNWISqw(siteNumber,parameterCd,
                   "2015-10-01","2020-10-01",
                   reshape=TRUE)

summary(dat1)

dat2 <- subset(dat1, select= c(agency_cd, site_no, 
                               DO.sat, depth, temp.water, light, discharge))



# colnames(TemeculaDat)[4] <- "temp.water"
# colnames(TemeculaDat)[6] <- "discharge"
# colnames(TemeculaDat)[8] <- "DO.obs"
# colnames(TemeculaDat)[19] <- "DO.sat"

qwData4$dec_lat_va <- Info$dec_lat_va
qwData4$dec_long_va <- Info$dec_long_va
qwData4$huc_cd <- Info$huc_cd
qwData4$drain_area_va <- Info$drain_area_va
qwData4$project_no <- Info$project_no
qwData4$station_nm <- Info$station_nm
qwData4$state_cd <- Info$state_cd
qwData4$Site_Lab <- "E. Clear" 


head(qwData4)

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
dat1$timestamp <- format(dat1$timestampPDT, tz="America/Los_Angeles", usetz=T)
as.POSIXct(dat1$timestamp)

# left join in Mesowest data
clim2 < read.csv("https://www.dropbox.com/s/gs8w50soclliami/KPSP.2020-10-01.csv?dl=1")
clim2$timestampUTC <- as.POSIXct(clim2$date_time, format = "%Y-%m-%dT%H:%M", tz="UTC")

clim2$timestamp2 <- round_date(clim2$timestampUTC, hour,unit="5 minutes")
clim2$timestamp <- format(clim2$timestamp2 , tz="America/Los_Angeles", usetz=TRUE)
as.POSIXct(clim2$timestamp)

# Merge data files
TemeculaDat <- left_join(dat1, clim2[c("timestampUTC", "pressure_set_1d", "air_temp_set_1")],
                         by = c("timestampPDT" = "timestampUTC"))
summary(TemeculaDat)

# Date range for both data sets = 201510010001 to 202009100001
TemeculaDat <- subset(TemeculaDat, timestampPDT >= as.POSIXct('2015-10-01 00:53:00') & 
                        timestampPDT <= as.POSIXct('2020-09-20 00:00:00'))
range(TemeculaDat$timestampPDT)
TemeculaDat$pressure_set_1d[is.nan(TemeculaDat$pressure_set_1d)] <- NA
TemeculaDat$air_temp_set_1[is.nan(TemeculaDat$air_temp_set_1)] <- NA

TemeculaDat$pressure_millibar <- c(TemeculaDat$pressure_set_1d * 0.01)

ggplot(TemeculaDat, aes(x=timestampPDT, y=(discharge))) +
  geom_point(alpha = 0.7)  +
  theme_classic()

### Stream Metabolizer estimates ###
TemeculaDat$depth <- calc_depth(Q=u(TemeculaDat$discharge, "m^3 s^-1"), f=u(0.36))
latitude <- c(33.4739175)
longitude <- c(-117.1422536) 

TemeculaDat$solar.time <- calc_solar_time(TemeculaDat$timestampPDT, longitude)
TemeculaDat$light <- calc_light(TemeculaDat$solar.time,
                                latitude,
                                longitude,
                                max.PAR = u(2326, "umol m^-2 s^-1"),
                                attach.units = is.unitted(TemeculaDat$solar.time))

TemeculaDat$DO_sat <- calc_DO_sat(TemeculaDat$tempC, 
                                  TemeculaDat$pressure_millibar, sal=0) 

names(TemeculaDat)
# colnames(TemeculaDat)[4] <- "temp.water"
# colnames(TemeculaDat)[6] <- "discharge"
# colnames(TemeculaDat)[8] <- "DO.obs"
# colnames(TemeculaDat)[19] <- "DO.sat"

TemeculaDatQ <- subset(TemeculaDat, select= c(solar.time, DO.obs, DO.sat, depth, temp.water, light, discharge))

# write.csv(TemeculaDatQ, paste0(outputDir, "/11044000TemeculaDat_Raw.csv")) 


