#' Data Organization for modeling metabolism for signals of smoke effect 
#' Last update: 2020-12-15
 
#' See http://usgs-r.github.io/streamMetabolizer for vignettes on the web.

## ---------------------------
# File path setup:
if (dir.exists('/Users/kellyloria/Documents/UNR_2020/Fall2020Projects/')){
  inputDir<- '/Users/kellyloria/Documents/UNR_2020/Fall2020Projects/'
  outputDir<- '/Users/kellyloria/Dropbox/Smoke\ Metabolism/smoke\ metabolism\ exploration/R_output'
}
## ---------------------------
## Load packages:
library(StreamMetabolism)
library(streamMetabolizer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(rstan)
library(unitted)
library(lubridate)
library(tidyverse)

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
d2 <- mw(service = 'timeseries', stid ='KPSP', vars = c('air_temp', 'pressure'), start = '201510010001', end = '202009100001', jsonsimplify= TRUE)
# parsing the nested lists
clim2 <- data.frame( lapply( d2$STATION$OBSERVATIONS, unlist) )
clim2$timestampUTC <- as.POSIXct(clim2$date_time, format = "%Y-%m-%dT%H:%M", tz="UTC")
#?round_date
clim2$timestamp2 <- round_date(clim2$timestampUTC, hour,unit="5 minutes")
clim2$timestamp <- format(clim2$timestamp2 , tz="America/Los_Angeles", usetz=TRUE)
as.POSIXct(clim2$timestamp)


# Merge data files  ## ***NEED to double timestamps still! ***
TemeculaDat <- left_join(dat1, clim2[c("timestampUTC", "pressure_set_1d", "air_temp_set_1")],
                           by = c("timestampPDT" = "timestampUTC"))
summary(TemeculaDat)

# Date range for both data sets = 201510010001 to 202009100001
TemeculaDat <- subset(TemeculaDat, timestampPDT >= as.POSIXct('2015-10-01 00:53:00') & 
                        timestampPDT <= as.POSIXct('2020-09-10 00:00:00'))
range(TemeculaDat$timestampPDT)
TemeculaDat$pressure_set_1d[is.nan(TemeculaDat$pressure_set_1d)] <- NA
TemeculaDat$air_temp_set_1[is.nan(TemeculaDat$air_temp_set_1)] <- NA



# Notes on data:
    # Discharge 15700.00 max does seem to be real and occur multiple years in Jan and Feb
    # DO above 15 appear to only occur in may 2018 and 2019

## ---------------------------
# Method 1 infill, then flag

# infill some data with rolling averages:
  # Water Temp
TemeculaDat$tempCMean <- rollapply(TemeculaDat$tempC, width=5,
                                       FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                       by.column=TRUE, partial=TRUE, fill=NA, align="center")

TemeculaDat$tempC[is.na(TemeculaDat$tempC)] <- as.numeric(ifelse(is.na(TemeculaDat$tempC), 
                                                                         paste(TemeculaDat$tempCMean),
                                                                         paste(TemeculaDat$tempC)))

  # Discharge:
TemeculaDat$dischargeMean <- rollapply(TemeculaDat$discharge, width=5,
                              FUN=function(x) mean(x, na.rm=TRUE), by=1,
                              by.column=TRUE, partial=TRUE, fill=NA, align="center")

TemeculaDat$discharge[is.na(TemeculaDat$discharge)] <- as.numeric(ifelse(is.na(TemeculaDat$discharge), 
                                                                         paste(TemeculaDat$dischargeMean),
                                                                         paste(TemeculaDat$discharge)))
  # Dissolved oxygen:
TemeculaDat$DOmgLMean <- rollapply(TemeculaDat$DOmgL, width=5,
                                       FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                       by.column=TRUE, partial=TRUE, fill=NA, align="center")

TemeculaDat$DOmgL[is.na(TemeculaDat$DOmgL)] <- as.numeric(ifelse(is.na(TemeculaDat$DOmgL), 
                                                                         paste(TemeculaDat$DOmgLMean),
                                                                         paste(TemeculaDat$DOmgL)))

  # Pressure:
TemeculaDat$pressureMean <- rollapply(TemeculaDat$pressure_set_1d, width=4500,
                                       FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                       by.column=TRUE, partial=TRUE, fill=NA, align="center")

TemeculaDat$pressure_set_1d[is.na(TemeculaDat$pressure_set_1d)] <- as.numeric(ifelse(is.na(TemeculaDat$pressure_set_1d), 
                                                                         paste(TemeculaDat$pressureMean),
                                                                         paste(TemeculaDat$pressure_set_1d)))
summary(TemeculaDat)
  # need to convert pascals to millibar
TemeculaDat$pressure_millibar <- c(TemeculaDat$pressure_set_1d * 0.01)
library(zoo)
###
# # Flag potential outliers:
# TemeculaDatQ = TemeculaDat %>% 
#   mutate(date=lubridate::date(timestampPDT))%>%
#   mutate(hour=lubridate::hour(timestampPDT))%>%
#   arrange(date, timestampPDT) %>%
#   group_by(date, hour) %>% 
#     # summarize for flow
#   mutate(
#     mnQ= rollapplyr(discharge, 30, FUN = mean, partial = T),
#     sdQ = rollapplyr(discharge, 30, FUN = sd, partial = T)) %>%
#   mutate(
#     loQ=c(mnQ- (3*sdQ)), 
#     hiQ=c(mnQ+ (3*sdQ)))%>%
#   # summarize for DO
#   mutate(mnDO=rollapply(DOmgL, width = 30, FUN = mean, partial = T),           # also filter out the NAs and >35s if you wanted to always have 15 values in your rolling window after removing bad values
#          sdDO=rollapply(DOmgL, width = 30, FUN = sd, partial = T)) %>%
#   mutate(
#     loDO=mnDO- (3*sdDO),
#     hiDO=mnDO+ (3*sdDO))%>%
#   # name flags for outliers:
#   mutate(
#     flag_discharge=
#       case_when( 
#         discharge<=loQ&!is.na(loQ) ~ 'o',
#         discharge>=hiQ&!is.na(hiQ) ~ 'o',
#         discharge<=(mean(discharge) + (sd(discharge)*3)) ~ 'p', 
#         TRUE ~ 'n')) %>%
#   mutate(
#     flag_DO=
#       case_when( 
#         DOmgL<=loDO&!is.na(loDO) ~ 'o',
#         DOmgL>=hiDO&!is.na(hiDO) ~ 'o',
#         DOmgL>= c(13.10366) ~ 'p',
#         DOmgL<= c(4.352055) ~ 'p',
#         TRUE ~ 'n'))

# 3.Check the flag 
p <- ggplot(TemeculaDatQ, aes(x=timestampPDT, y=(discharge), colour =as.factor(flag_discharge))) +
  geom_point(alpha = 0.7)  +
  theme_classic()


p <- ggplot(TemeculaDatQ1, aes(x=timestampPDT, y=(DOmgL), colour =as.factor(flag_DO))) +
  geom_point(alpha = 0.7)  +
  theme_classic()
# 
# 
# TemeculaDatQ1 <- subset(TemeculaDatQ, flag_discharge == "n" &
#                           flag_DO == "n")
# summary(night_check$light) # feel like solar 
# 
# TemeculaDatQ1$depth <- calc_depth(Q=u(TemeculaDatQ1$discharge, "m^3 s^-1"), f=u(0.36))
# #33.4739175,-117.1422536
# latitude <- c(33.4739175)
# longitude <- c(-117.1422536) 
# 
# TemeculaDatQ1$solar.time <- calc_solar_time(TemeculaDatQ1$timestampPDT, longitude)
# TemeculaDatQ1$light <- calc_light(TemeculaDatQ1$solar.time,
#                                 latitude,
#                                 longitude,
#                                 max.PAR = u(2326, "umol m^-2 s^-1"),
#                                 attach.units = is.unitted(TemeculaDatQ1$solar.time)
# )
# 
# TemeculaDatQ1$DO_sat <- calc_DO_sat(TemeculaDatQ1$tempC, 
#                                     TemeculaDatQ1$pressure_millibar, sal=0) 
# 
# names(TemeculaDatQ1)
# 
# colnames(TemeculaDatQ1)[4] <- "temp.water"
# colnames(TemeculaDatQ1)[6] <- "discharge"
# colnames(TemeculaDatQ1)[8] <- "DO.obs"
# colnames(TemeculaDatQ1)[35] <- "DO.sat"
# 
# TemCAdat <- subset(TemeculaDatQ1, select= c(solar.time, DO.obs, DO.sat, depth, temp.water, light, discharge))
# 
# 
# dim(TemeculaDatQ1)
# 
# # Minor metabolism exploration...
# TemCAdat %>% unitted::v() %>%
#   mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
#   select(solar.time, starts_with('DO')) %>%
#   gather(type, DO.value, starts_with('DO')) %>%
#   mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
#   ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() +
#   facet_grid(units ~ ., scale='free_y') + theme_bw() +
#   scale_color_discrete('variable')
# 
# labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)', light='PAR\n(umol m^-2 s^-1)')
# TemCAdat %>% unitted::v() %>%
#   select(solar.time, depth, temp.water, light) %>%
#   gather(type, value, depth, temp.water, light) %>%
#   mutate(
#     type=ordered(type, levels=c('depth','temp.water','light')),
#     units=ordered(labels[type], unname(labels))) %>%
#   ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() +
#   facet_grid(units ~ ., scale='free_y') + theme_bw() +
#   scale_color_discrete('variable')
# 
# TemCAdat2 <- na.omit(TemCAdat)
# 
# ?calc_DO_sat
# ?as.POSIXct
# 
# mm <- metab(specs(mm_name('mle')), data=TemCAdat2, info='my info')
# predict_metab(mm)
# get_info(mm)
# get_fitting_time(mm)
# 
# BW_ModelDOPlot <- plot_DO_preds(predict_DO(mm))
# #?plot_DO_preds
# BW_ModelMetabPlot <- plot_metab_preds(mm)



## ---------------------------
# Method 2 NA remove:
summary(TemeculaDat)
TemeculaDatNA <- na.omit(TemeculaDat)
summary(TemeculaDatNA)
dim(TemeculaDatNA)

TemeculaDatNA$pressure_millibar <- c(TemeculaDatNA$pressure_set_1d * 0.01)

p <- ggplot(TemeculaDatNA, aes(x=timestampPDT, y=(DOmgL))) +
  geom_point(alpha = 0.7)  +
  theme_classic()

TemeculaDatNA$depth <- calc_depth(Q=u(TemeculaDatNA$discharge, "m^3 s^-1"), f=u(0.36))
#33.4739175,-117.1422536
latitude <- c(33.4739175)
longitude <- c(-117.1422536) 

TemeculaDatNA$solar.time <- calc_solar_time(TemeculaDatNA$timestampPDT, longitude)
TemeculaDatNA$light <- calc_light(TemeculaDatNA$solar.time,
                                  latitude,
                                  longitude,
                                  max.PAR = u(2326, "umol m^-2 s^-1"),
                                  attach.units = is.unitted(TemeculaDatNA$solar.time))

TemeculaDatNA$DO_sat <- calc_DO_sat(TemeculaDatNA$tempC, 
                                    TemeculaDatNA$pressure_millibar, sal=0) 

names(TemeculaDatNA)

colnames(TemeculaDatNA)[4] <- "temp.water"
colnames(TemeculaDatNA)[6] <- "discharge"
colnames(TemeculaDatNA)[8] <- "DO.obs"
colnames(TemeculaDatNA)[19] <- "DO.sat"

TemCAdat_NA <- subset(TemeculaDatNA, select= c(solar.time, DO.obs, DO.sat, depth, temp.water, light, discharge))

# Minor metabolism exploration...
TemCAdat_NA %>% unitted::v() %>%
  mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
  select(solar.time, starts_with('DO')) %>%
  gather(type, DO.value, starts_with('DO')) %>%
  mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
  ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() +
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)', light='PAR\n(umol m^-2 s^-1)')
TemCAdat_NA %>% unitted::v() %>%
  select(solar.time, depth, temp.water, light) %>%
  gather(type, value, depth, temp.water, light) %>%
  mutate(
    type=ordered(type, levels=c('depth','temp.water','light')),
    units=ordered(labels[type], unname(labels))) %>%
  ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() +
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

#TemCAdat2 <- na.omit(TemCAdat)

?calc_DO_sat
?as.POSIXct

mm <- metab(specs(mm_name('mle')), data=TemCAdat_NA, info='my info')
predict_metab(mm)
get_info(mm)
get_fitting_time(mm)

BW_ModelDOPlot <- plot_DO_preds(predict_DO(mm))
#?plot_DO_preds
BW_ModelMetabPlot <- plot_metab_preds(mm)


##
## ---------------------------
# Method 3 coarse flag then infill:
# Flag potential outliers:

TemeculaDatQ_2 <- subset(TemeculaDatQ, DOmgL < 12.00988 ) 
TemeculaDatQ_3 <- subset(TemeculaDatQ_2, DOmgL > 5.466516) 
TemeculaDatQ_4 <- subset(TemeculaDatQ_3, discharge < 799.9161) 

p <- ggplot(TemeculaDatQ_4, aes(x=timestampPDT, y=(DOmgL)), color=flag_DO) +
  geom_point(alpha = 0.7)  +
  theme_classic()

p <- ggplot(TemeculaDatQ_4, aes(x=timestampPDT, y=(discharge)), color=flag_discharge) +
  geom_point(alpha = 0.7)  +
  theme_classic()

# Date range for both data sets = 201510010001 to 202009100001
TemeculaDatQ_4 <- subset(TemeculaDatQ_4, timestampPDT >= as.POSIXct('2015-10-01 00:53:00') & 
                        timestampPDT <= as.POSIXct('2020-09-10 00:00:00'))
range(TemeculaDatQ_4$timestampPDT)
TemeculaDatQ_4$pressure_set_1d[is.nan(TemeculaDatQ_4$pressure_set_1d)] <- NA
TemeculaDatQ_4$air_temp_set_1[is.nan(TemeculaDatQ_4$air_temp_set_1)] <- NA

## ---------------------------
# Infill, then flag

# infill some data with rolling averages:
# Water Temp
TemeculaDatQ_4$tempCMean <- rollapply(TemeculaDatQ_4$tempC, width=5,
                                   FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                   by.column=TRUE, partial=TRUE, fill=NA, align="center")

TemeculaDatQ_4$tempC[is.na(TemeculaDatQ_4$tempC)] <- as.numeric(ifelse(is.na(TemeculaDatQ_4$tempC), 
                                                                 paste(TemeculaDatQ_4$tempCMean),
                                                                 paste(TemeculaDatQ_4$tempC)))

# Discharge:
TemeculaDatQ_4$dischargeMean <- rollapply(TemeculaDatQ_4$discharge, width=5,
                                       FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                       by.column=TRUE, partial=TRUE, fill=NA, align="center")

TemeculaDatQ_4$discharge[is.na(TemeculaDatQ_4$discharge)] <- as.numeric(ifelse(is.na(TemeculaDatQ_4$discharge), 
                                                                         paste(TemeculaDatQ_4$dischargeMean),
                                                                         paste(TemeculaDatQ_4$discharge)))
# Dissolved oxygen:
TemeculaDatQ_4$DOmgLMean <- rollapply(TemeculaDatQ_4$DOmgL, width=5,
                                   FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                   by.column=TRUE, partial=TRUE, fill=NA, align="center")

TemeculaDatQ_4$DOmgL[is.na(TemeculaDatQ_4$DOmgL)] <- as.numeric(ifelse(is.na(TemeculaDatQ_4$DOmgL), 
                                                                 paste(TemeculaDatQ_4$DOmgLMean),
                                                                 paste(TemeculaDatQ_4$DOmgL)))

# Pressure:
TemeculaDatQ_4$pressureMean <- rollapply(TemeculaDatQ_4$pressure_set_1d, width=4500,
                                      FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                      by.column=TRUE, partial=TRUE, fill=NA, align="center")

TemeculaDatQ_4$pressure_set_1d[is.na(TemeculaDatQ_4$pressure_set_1d)] <- as.numeric(ifelse(is.na(TemeculaDatQ_4$pressure_set_1d), 
                                                                                     paste(TemeculaDatQ_4$pressureMean),
                                                                                     paste(TemeculaDatQ_4$pressure_set_1d)))
summary(TemeculaDatQ_4)
# need to convert pascals to millibar
TemeculaDatQ_4$pressure_millibar <- c(TemeculaDatQ_4$pressure_set_1d * 0.01)

TemeculaDatQ_4$depth <- calc_depth(Q=u(TemeculaDatQ_4$discharge, "m^3 s^-1"), f=u(0.36))
#33.4739175,-117.1422536
latitude <- c(33.4739175)
longitude <- c(-117.1422536) 

TemeculaDatQ_4$solar.time <- calc_solar_time(TemeculaDatQ_4$timestampPDT, longitude)
TemeculaDatQ_4$light <- calc_light(TemeculaDatQ_4$solar.time,
                                  latitude,
                                  longitude,
                                  max.PAR = u(2326, "umol m^-2 s^-1"),
                                  attach.units = is.unitted(TemeculaDatQ_4$solar.time))

TemeculaDatQ_4$DO_sat <- calc_DO_sat(TemeculaDatQ_4$tempC, 
                                    TemeculaDatQ_4$pressure_millibar, sal=0) 

names(TemeculaDatQ_4)

colnames(TemeculaDatQ_4)[4] <- "temp.water"
colnames(TemeculaDatQ_4)[6] <- "discharge"
colnames(TemeculaDatQ_4)[8] <- "DO.obs"
colnames(TemeculaDatQ_4)[35] <- "DO.sat"

TemCAdat_Q5 <- subset(TemeculaDatQ_4, select= c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light", "discharge"))

# Minor metabolism exploration...
TemCAdat_Q5 %>% unitted::v() %>%
  mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
  select(solar.time, starts_with('DO')) %>%
  gather(type, DO.value, starts_with('DO')) %>%
  mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
  ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() +
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)', light='PAR\n(umol m^-2 s^-1)')
TemCAdat_Q5 %>% unitted::v() %>%
  select(solar.time, depth, temp.water, light) %>%
  gather(type, value, depth, temp.water, light) %>%
  mutate(
    type=ordered(type, levels=c('depth','temp.water','light')),
    units=ordered(labels[type], unname(labels))) %>%
  ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() +
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

mm <- metab(specs(mm_name('mle')), data=TemCAdat_NA, info='my info')
predict_metab(mm)
get_info(mm)
get_fitting_time(mm)

BW_ModelMetabPlot <- plot_metab_preds(mm)


















dat2 <- read.csv("https://www.dropbox.com/s/eochnr2z8v3qsip/nwis.11128500.csv?dl=1")
summary(dat2) # SANTA YNEZ R A SOLVANG CA




# write.csv(dat1, paste0(outputDir, "TESTFILE.csv")) # complied data file of all DO sensors along buoy line