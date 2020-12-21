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
library(zoo)

## ---------------------------
## Load packages:
dat2 <- read.csv("https://www.dropbox.com/s/eochnr2z8v3qsip/nwis.11128500.csv?dl=1")
summary(dat2) # SANTA YNEZ R A SOLVANG CA

dat2$timestampPDT <- as.POSIXct(dat2$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
summary(dat2)

# near (34.58498706, -120.1445926)

KIZAdat <- read_csv("SM_mesonetDat/KIZAdat.csv")
#Parsed with column specification:
  cols(
    Station_ID = col_character(),
    Date_Time = col_datetime(format = ""),
    altimeter_set_1_pascals = col_double(),
    air_temp_set_1_C = col_double(),
    pressure_set_1d_pascals = col_double()
  )

summary(KIZAdat)  
# Fix timestamp (round to nearest 5 minute):
KIZAdat$timestamp <- round_date(KIZAdat$Date_Time, hour,unit="5 minutes")
#KIZAdat$timestamp1 <- format(KIZAdat$timestamp , tz="America/Los_Angeles", usetz=TRUE)
KIZAdat$timestamp2 <- as.POSIXct(KIZAdat$timestamp1, format = "%Y-%m-%d %H:%M", tz="America/Los_Angeles", usetz=TRUE)

names(KIZAdat)
# Merge data files 
SantaYnezDat <- left_join(dat2, KIZAdat[c("timestamp2", "altimeter_set_1_pascals", "air_temp_set_1_C")],
                         by = c("timestampPDT" = "timestamp2"))
summary(SantaYnezDat)

SantaYnezDat <- subset(SantaYnezDat, timestampPDT >= as.POSIXct('2015-10-01 00:53:00') & 
                        timestampPDT <= as.POSIXct('2020-09-10 00:00:00'))
range(SantaYnezDat$timestampPDT)
summary(SantaYnezDat)

# # # Flag potential outliers:
SantaYnezDatQ = SantaYnezDat %>%
  mutate(date=lubridate::date(timestampPDT))%>%
  mutate(hour=lubridate::hour(timestampPDT))%>%
  arrange(date, timestampPDT) %>%
  group_by(date) %>%
    # summarize for flow
  mutate(
    mnQ= rollapplyr(discharge, 30, FUN = mean, partial = T),
    sdQ = rollapplyr(discharge, 30, FUN = sd, partial = T)) %>%
  mutate(
    loQ=c(mnQ- (3*sdQ)),
    hiQ=c(mnQ+ (3*sdQ)))%>%
  # summarize for DO
  mutate(mnDO=rollapply(DOmgL, width = 30, FUN = mean, partial = T),           # also filter out the NAs and >35s if you wanted to always have 15 values in your rolling window after removing bad values
         sdDO=rollapply(DOmgL, width = 30, FUN = sd, partial = T)) %>%
  mutate(
    loDO=mnDO- (3*sdDO),
    hiDO=mnDO+ (3*sdDO))%>%
  # name flags for outliers:
  mutate(
    flag_discharge=
      case_when(
        discharge>loQ&!is.na(loQ) ~ 'o',
        discharge<hiQ&!is.na(hiQ) ~ 'o',
        TRUE ~ 'n')) %>%
  mutate(
    flag_DO=
      case_when(
        DOmgL<=loDO&!is.na(loDO) ~ 'o',
        DOmgL>=hiDO&!is.na(hiDO) ~ 'o',
        TRUE ~ 'n'))

# plot to check flags 
p <- ggplot(SantaYnezDat, aes(x=timestampPDT, y=(discharge))) +
  geom_point(alpha = 0.7)  +
  theme_classic()

p <- ggplot(SantaYnezDatQ, aes(x=timestampPDT, y=(DOmgL), color=flag_DO)) +
  geom_point(alpha = 0.7)  +
  theme_classic()

pb <- ggplot(SantaYnezDatQ, aes(x=timestampPDT, y=(discharge), color=flag_discharge)) +
  geom_point(alpha = 0.7)  +
  theme_classic()

# Remove outliers:
SantaYnezDatQ1 <- subset(SantaYnezDatQ, flag_discharge == "n" & flag_DO == "n")
summary(SantaYnezDatQ1)

pb <- ggplot(SantaYnezDatQ1, aes(x=timestampPDT, y=(discharge), color=flag_discharge)) +
  geom_point(alpha = 0.7)  +
  theme_classic()

SantaYnezDatQ2 <- subset(SantaYnezDatQ1, discharge < 1000)
summary(SantaYnezDatQ2)

pb <- ggplot(SantaYnezDatQ2, aes(x=timestampPDT, y=(discharge), color=flag_discharge)) +
  geom_point(alpha = 0.7)  +
  theme_classic()

###
###
# Fill in NA's
SantaYnezDatQ2$tempC[is.nan(SantaYnezDatQ2$tempC)] <- NA
# SantaYnezDatQ2$pressure_set_1d[is.nan(TemeculaDatQ_4$pressure_set_1d)] <- NA
# TemeculaDatQ_4$air_temp_set_1[is.nan(TemeculaDatQ_4$air_temp_set_1)] <- NA

SantaYnezDatQ2$tempCMean <- rollapply(SantaYnezDatQ2$tempC, width=5,
                                   FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                   by.column=TRUE, partial=TRUE, fill=NA, align="center")

SantaYnezDatQ2$tempC[is.na(SantaYnezDatQ2$tempC)] <- as.numeric(ifelse(is.na(SantaYnezDatQ2$tempC), 
                                                                 paste(SantaYnezDatQ2$tempCMean),
                                                                 paste(SantaYnezDatQ2$tempC)))

# Discharge:
SantaYnezDatQ2$dischargeMean <- rollapply(SantaYnezDatQ2$discharge, width=5,
                                       FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                       by.column=TRUE, partial=TRUE, fill=NA, align="center")

SantaYnezDatQ2$discharge[is.na(SantaYnezDatQ2$discharge)] <- as.numeric(ifelse(is.na(SantaYnezDatQ2$discharge), 
                                                                         paste(SantaYnezDatQ2$dischargeMean),
                                                                         paste(SantaYnezDatQ2$discharge)))
# Dissolved oxygen:
SantaYnezDatQ2$DOmgLMean <- rollapply(SantaYnezDatQ2$DOmgL, width=5,
                                   FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                   by.column=TRUE, partial=TRUE, fill=NA, align="center")

SantaYnezDatQ2$DOmgL[is.na(SantaYnezDatQ2$DOmgL)] <- as.numeric(ifelse(is.na(SantaYnezDatQ2$DOmgL), 
                                                                 paste(SantaYnezDatQ2$DOmgLMean),
                                                                 paste(SantaYnezDatQ2$DOmgL)))

# Pressure:
SantaYnezDatQ2$pressureMean <- rollapply(SantaYnezDatQ2$altimeter_set_1_pascals, width=1500,
                                      FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                      by.column=TRUE, partial=TRUE, fill=NA, align="center")

SantaYnezDatQ2$altimeter_set_1_pascals[is.na(SantaYnezDatQ2$altimeter_set_1_pascals)] <- as.numeric(ifelse(is.na(SantaYnezDatQ2$altimeter_set_1_pascals), 
                                                                                     paste(SantaYnezDatQ2$pressureMean),
                                                                                     paste(SantaYnezDatQ2$altimeter_set_1_pascals)))
summary(SantaYnezDatQ2)
# need to convert pascals to millibar
SantaYnezDatQ2$pressure_millibar <- c(SantaYnezDatQ2$altimeter_set_1_pascals * 0.01)

# Estimate metabolism parameters:
SantaYnezDatQ2$depth <- calc_depth(Q=u(SantaYnezDatQ2$discharge, "m^3 s^-1"), f=u(0.36))

#(34.58498706, -120.1445926)
latitude <- c(34.58498706)
longitude <- c(-120.1445926) 

SantaYnezDatQ2$solar.time <- calc_solar_time(SantaYnezDatQ2$timestampPDT, longitude)
SantaYnezDatQ2$light <- calc_light(SantaYnezDatQ2$solar.time,
                                   latitude,
                                   longitude,
                                   max.PAR = u(2326, "umol m^-2 s^-1"),
                                   attach.units = is.unitted(SantaYnezDatQ2$solar.time))

SantaYnezDatQ2$DO_sat <- calc_DO_sat(SantaYnezDatQ2$tempC, 
                                     SantaYnezDatQ2$pressure_millibar, sal=0) 

names(SantaYnezDatQ2)

colnames(SantaYnezDatQ2)[7] <- "temp.water"
colnames(SantaYnezDatQ2)[5] <- "discharge"
colnames(SantaYnezDatQ2)[9] <- "DO.obs"
colnames(SantaYnezDatQ2)[34] <- "DO.sat"

SantaYnezDatQ3 <- subset(SantaYnezDatQ2, select= c("solar.time", "DO.obs", "DO.sat", "depth", "temp.water", "light", "discharge"))

# Minor metabolism exploration...
SantaYnezDatQ3 %>% unitted::v() %>%
  mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
  select(solar.time, starts_with('DO')) %>%
  gather(type, DO.value, starts_with('DO')) %>%
  mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
  ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() +
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)', light='PAR\n(umol m^-2 s^-1)')
SantaYnezDatQ3 %>% unitted::v() %>%
  select(solar.time, depth, temp.water, light) %>%
  gather(type, value, depth, temp.water, light) %>%
  mutate(
    type=ordered(type, levels=c('depth','temp.water','light')),
    units=ordered(labels[type], unname(labels))) %>%
  ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() +
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

summary(SantaYnezDatQ3)

mm <- metab(specs(mm_name('mle')), data=SantaYnezDatQ3, info='my info')
predict_metab(mm)
get_info(mm)
get_fitting_time(mm)

BW_ModelMetabPlot <- plot_metab_preds(mm)


## ---------------------------
## Load packages:
dat3 <- read.csv("https://www.dropbox.com/s/d1taddci0q0lvqk/nwis.11273400.csv?dl=1")
summary(dat3)

dat3$timestampPDT <- as.POSIXct(dat3$datetime, format = "%m/%d/%Y %H:%M", tz="America/Los_Angeles", usetz=TRUE)
summary(dat3)

KMOD <- read_csv("SM_mesonetDat/KMOD.csv")
  cols(
    Station_ID = col_character(),
    Date_Time = col_datetime(format = ""),
    altimeter_set_1 = col_double(),
    air_temp_set_1 = col_double(),
    pressure_set_1d = col_double()
  )

summary(KMOD)

##
# Fix timestamp (round to nearest 5 minute):
KMOD$timestamp <- round_date(KMOD$Date_Time, hour,unit="5 minutes")
KMOD$timestamp1 <- format(KMOD$timestamp , tz="America/Los_Angeles", usetz=TRUE)
KMOD$timestamp2 <- as.POSIXct(KMOD$timestamp1, format = "%Y-%m-%d %H:%M", tz="America/Los_Angeles", usetz=TRUE)

names(KIZAdat)
# Merge data files 
SanJoaquin <- left_join(dat3, KMOD[c("timestamp2", "air_temp_set_1", "pressure_set_1d")],
                          by = c("timestampPDT" = "timestamp2"))
summary(SanJoaquin)

SanJoaquin <- subset(SanJoaquin, timestampPDT >= as.POSIXct('2015-10-01 00:53:00') & 
                         timestampPDT <= as.POSIXct('2020-09-10 00:00:00'))
range(SanJoaquin$timestampPDT)
summary(SanJoaquin)

# # # Flag potential outliers:
SanJoaquinDatQ = SanJoaquin %>%
  mutate(date=lubridate::date(timestampPDT))%>%
  mutate(hour=lubridate::hour(timestampPDT))%>%
  arrange(date, timestampPDT) %>%
  group_by(date) %>%
  # summarize for flow
  mutate(
    mnQ= rollapplyr(discharge, 30, FUN = mean, partial = T),
    sdQ = rollapplyr(discharge, 30, FUN = sd, partial = T)) %>%
  mutate(
    loQ=c(mnQ- (3*sdQ)),
    hiQ=c(mnQ+ (3*sdQ)))%>%
  # summarize for DO
  mutate(mnDO=rollapply(DOmgL, width = 30, FUN = mean, partial = T),           # also filter out the NAs and >35s if you wanted to always have 15 values in your rolling window after removing bad values
         sdDO=rollapply(DOmgL, width = 30, FUN = sd, partial = T)) %>%
  mutate(
    loDO=mnDO- (3*sdDO),
    hiDO=mnDO+ (3*sdDO))%>%
  # name flags for outliers:
  mutate(
    flag_discharge=
      case_when(
        discharge>loQ&!is.na(loQ) ~ 'o',
        discharge<hiQ&!is.na(hiQ) ~ 'o',
        TRUE ~ 'n')) %>%
  mutate(
    flag_DO=
      case_when(
        DOmgL<=loDO&!is.na(loDO) ~ 'o',
        DOmgL>=hiDO&!is.na(hiDO) ~ 'o',
        TRUE ~ 'n'))

# plot to check flags 
p <- ggplot(SanJoaquinDatQ, aes(x=timestampPDT, y=(discharge))) +
  geom_point(alpha = 0.7)  +
  theme_classic()

p <- ggplot(SantaYnezDatQ, aes(x=timestampPDT, y=(DOmgL), color=flag_DO)) +
  geom_point(alpha = 0.7)  +
  theme_classic()

pb <- ggplot(SanJoaquinDatQ, aes(x=timestampPDT, y=(discharge), color=flag_discharge)) +
  geom_point(alpha = 0.7)  +
  theme_classic()

# Remove outliers:
SantaYnezDatQ1 <- subset(SantaYnezDatQ, flag_discharge == "n" & flag_DO == "n")
summary(SantaYnezDatQ1)

pb <- ggplot(SantaYnezDatQ1, aes(x=timestampPDT, y=(discharge), color=flag_discharge)) +
  geom_point(alpha = 0.7)  +
  theme_classic()

SantaYnezDatQ2 <- subset(SantaYnezDatQ1, discharge < 1000)
summary(SantaYnezDatQ2)

pb <- ggplot(SantaYnezDatQ2, aes(x=timestampPDT, y=(discharge), color=flag_discharge)) +
  geom_point(alpha = 0.7)  +
  theme_classic()


  