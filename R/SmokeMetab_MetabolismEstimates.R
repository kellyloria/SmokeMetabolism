## Smoke stream metabolism - J.R. Blaszczak

## Install streammetabolizer if not already done
## http://usgs-r.github.io/streamMetabolizer/articles/installation.html
#remotes::install_github('appling/unitted')
#remotes::install_github("USGS-R/streamMetabolizer")
#devtools::find_rtools()

## Load all packages
lapply(c("plyr","dplyr","ggplot2","cowplot","tidyverse",
         "streamMetabolizer","rstan","lubridate","parallel"), require, character.only=T) #any FALSE need to be installed
theme_set(theme_bw())

#####################
## Import data
#####################
getwd()
setwd("C:/Users/jblaszczak/Dropbox/Smoke Metabolism/smoke metabolism exploration/R_output/InfilledDat")

TS_files <- list.files()

TS_dat <- ldply(TS_files, function(filename) {
  dum = read.csv(filename, header = TRUE)
  dum$ID <- filename
  dum$NWIS_site <- substr(filename, 0, 9)
  return(dum)
})
head(TS_dat)

#######################
## Format
######################
## Subset
TS_dat <- TS_dat[,c("solar.time","DO.obs","DO.sat","depth","temp.water","light","discharge","NWIS_site")]
## Change to date
TS_dat$solar.time <- as.POSIXct(as.character(TS_dat$solar.time), format="%Y-%m-%d %H:%M:%S", tz="UTC")
## Split into list
TS_list <- split(TS_dat, TS_dat$NWIS_site)
## Get rid of NWIS_site column
TS_list <- lapply(TS_list, function(x) subset(x, select=-c(NWIS_site)))

## Subset data to only after 06-01-2020
TS_list <- lapply(TS_list, function(x) subset(x, solar.time > as.POSIXct("2020-06-01 00:00:00")))
lapply(TS_list, function(x) tail(x)) ## TS ends 2020-09-09


########################
## Visualize
########################
dat <- TS_list[[1]]

dat %>% unitted::v() %>%
  mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
  select(solar.time, starts_with('DO')) %>%
  gather(type, DO.value, starts_with('DO')) %>%
  mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
  ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() + 
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)', light='PAR\n(umol m^-2 s^-1)', discharge='Q\n(cms)')
dat %>% unitted::v() %>%
  select(solar.time, depth, temp.water, light, discharge) %>%
  gather(type, value, depth, temp.water, light, discharge) %>%
  mutate(
    type=ordered(type, levels=c('depth','temp.water','light','discharge')),
    units=ordered(labels[type], unname(labels))) %>%
  ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() + 
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')





########################
## Set model and specs
#######################
## Set bayes model
bayes_name_new <- mm_name(type='bayes', pool_K600="none",
                          err_obs_iid=TRUE, err_proc_iid = TRUE,
                          ode_method = "trapezoid", deficit_src='DO_mod', engine='stan')
## set model specs
bayes_specs <- specs(bayes_name_new)
bayes_specs ## modify specs next time

######################
## Run model
#####################
## data without discharge for no pooling
TS_list_noQ <- lapply(TS_list, function(x) subset(x, select=-c(discharge)))

## test run
mm <- metab(bayes_specs, data=na.omit(TS_list_noQ[[1]]))
## run model with every site
mm_list <- lapply(TS_list_noQ, function(x) metab(bayes_specs, data=na.omit(x)))

########################
## Inspect model
######################
predict_metab(mm)
plot_metab_preds(mm)

























