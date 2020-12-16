#' Data Organization for modeling metabolism for signals of smoke effect 
#' Last update: 2020-12-15
 
 
#' See http://usgs-r.github.io/streamMetabolizer for vignettes on the web.

#' Create a folder for data output:
#' https://www.dropbox.com/sh/li9v3dyg0jakc7r/AADJGVzv4Io21MVu9sy_sVraa?dl=1

# File path setup:

if (dir.exists('/Users/kellyloria/Documents/UNR_2020/Fall2020Projects/')){
  inputDir<- '/Users/kellyloria/Documents/UNR_2020/Fall2020Projects/'
  outputDir<- '/Users/kellyloria/Dropbox/Smoke\ Metabolism/smoke\ metabolism\ exploration/R_output'
}

library(lubridate)
library(tidyverse)

#' Link from dropbox
# https://www.dropbox.com/s/ujyad748ka7flli/nwis.11044000.csv?dl=0
# change dl=0 to dl=1
dat1 <- read.csv( "https://www.dropbox.com/s/ujyad748ka7flli/nwis.11044000.csv?dl=1" )
summary(dat1) 
# SANTA MARGARITA R NR TEMECULA CA https://waterdata.usgs.gov/nwis/inventory/?site_no=11044000&agency_cd=USGS
range(dat1$datetime)
dat1$timestampPDT <- as.POSIXct(dat1$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
dat1$timestamp <- format(dat1$timestampPDT, tz="America/Los_Angeles", usetz=T)
as.POSIXct(dat1$timestamp)

#left join in meso west data
d2 <- mw(service = 'timeseries', stid ='KPSP', vars = c('air_temp', 'pressure'), start = '201510010001', end = '202009100001', jsonsimplify= TRUE)
# parsing the nested lists
clim2 <- data.frame( lapply( d2$STATION$OBSERVATIONS, unlist) )
clim2$timestampUTC <- as.POSIXct(clim2$date_time, format = "%Y-%m-%dT%H:%M", tz="UTC")
#?round_date
clim2$timestamp2 <- round_date(clim2$timestampUTC, hour,unit="5 minutes")
clim2$timestamp <- format(clim2$timestamp2 , tz="America/Los_Angeles", usetz=TRUE)
as.POSIXct(clim2$timestamp)


# Merge data files  ## NEED to double timestamps still!
TemeculaDat <- left_join(dat1, clim2[c("timestampUTC", "pressure_set_1d", "air_temp_set_1")],
                           by = c("timestampPDT" = "timestampUTC"))
summary(TemeculaDat)

# Date range for both datasets = 201510010001 to 202009100001
TemeculaDat <- subset(TemeculaDat, timestampPDT >= as.POSIXct('2015-10-01 00:53:00') & 
                        timestampPDT <= as.POSIXct('2020-09-10 00:00:00'))
range(TemeculaDat$timestampPDT)
TemeculaDat$pressure_set_1d[is.nan(TemeculaDat$pressure_set_1d)] <- NA
TemeculaDat$air_temp_set_1[is.nan(TemeculaDat$air_temp_set_1)] <- NA

# infill some means

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
TemeculaDat$pressureMean <- rollapply(TemeculaDat$pressure_set_1d, width=5000,
                                       FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                       by.column=TRUE, partial=TRUE, fill=NA, align="center")

TemeculaDat$pressure_set_1d[is.na(TemeculaDat$pressure_set_1d)] <- as.numeric(ifelse(is.na(TemeculaDat$pressure_set_1d), 
                                                                         paste(TemeculaDat$pressureMean),
                                                                         paste(TemeculaDat$pressure_set_1d)))
summary(TemeculaDat)

#need to convert pascals to millibar
TemeculaDat$pressure_millibar <- c(TemeculaDat$pressure_set_1d * 0.01)

?calc_depth
TemeculaDat$depth <- calc_depth(Q=u(TemeculaDat$discharge, "m^3 s^-1"), f=u(0.36))
#33.4739175,-117.1422536
latitude <- c(33.4739175)
longitude <- c(-117.1422536) 

TemeculaDat$solar.time <- calc_solar_time(TemeculaDat$timestampPDT, longitude)
TemeculaDat$light <- calc_light(TemeculaDat$solar.time,
                                  latitude,
                                  longitude,
                                  max.PAR = u(2326, "umol m^-2 s^-1"),
                                  attach.units = is.unitted(TemeculaDat$solar.time)
)
#plot(TemeculaDat$timestampPDT, TemeculaDat$light)
TemeculaDat$DO_sat <- calc_DO_sat(TemeculaDat$tempC, 
                                  TemeculaDat$pressure_millibar, sal=0) 


# data selection 
metab_inputs('mle', 'data')

###
# Get data in correct name and column form
names(TemeculaDat)

# colnames(TemeculaDat)[4] <- "temp.water"
# colnames(TemeculaDat)[6] <- "discharge"
# colnames(TemeculaDat)[8] <- "DO.obs"
# colnames(TemeculaDat)[23] <- "DO.sat"

TemCAdat <- subset(TemeculaDat, select= c(solar.time, DO.obs, DO.sat, depth, temp.water, light, discharge))


# write.csv(TemCAdat, paste0(outputDir, "/11044000Temecula_MetabDat.csv")) # complied data file 


TemCAdat %>% unitted::v() %>%
  mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
  select(solar.time, starts_with('DO')) %>%
  gather(type, DO.value, starts_with('DO')) %>%
  mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
  ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() +
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)', light='PAR\n(umol m^-2 s^-1)')
TemCAdat %>% unitted::v() %>%
  select(solar.time, depth, temp.water, light) %>%
  gather(type, value, depth, temp.water, light) %>%
  mutate(
    type=ordered(type, levels=c('depth','temp.water','light')),
    units=ordered(labels[type], unname(labels))) %>%
  ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() +
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')


?calc_DO_sat
?as.POSIXct

mm <- metab(specs(mm_name('mle')), data=TemCAdat, info='my info')
predict_metab(mm)
get_info(mm)
get_fitting_time(mm)

BW_ModelDOPlot <- plot_DO_preds(predict_DO(mm))
#?plot_DO_preds
BW_ModelMetabPlot <- plot_metab_preds(mm)



?strptime










dat2 <- read.csv("https://www.dropbox.com/s/eochnr2z8v3qsip/nwis.11128500.csv?dl=1")
summary(dat2) # SANTA YNEZ R A SOLVANG CA




# write.csv(dat1, paste0(outputDir, "TESTFILE.csv")) # complied data file of all DO sensors along buoy line