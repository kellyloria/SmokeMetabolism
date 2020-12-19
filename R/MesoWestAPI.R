#https://api.synopticdata.com/v2/stations/latest?radius=KHOU,50&limit=10&vars=air_temp&within=100&token=2c9fda8f7b564ae185110f7c0e972cf2

library(devtools)
install_github('fickse/mesowest')
library(mesowest)
?mesowest

# requestToken(apikey = "KymyZRyqaXp7uDR73MYYMN15e15NiDHfMnfL7HsxPl")
mw(service = 'metadata', complete=1, state='CA', county='Riverside') # near (33.4739175,-117.1422536)
# KPSP
# request hourly timeseries data)
d <- mw(service = 'timeseries', stid ='CI062', vars = c('air_temp', 'pressure_1500_meter'), start = '201510010001', end = '202009100001', jsonsimplify= TRUE)
# parsing the nested lists
clim <- data.frame( lapply( d$STATION$OBSERVATIONS, unlist) )

d2 <- mw(service = 'timeseries', stid ='KPSP', vars = c('air_temp', 'pressure'), start = '201510010001', end = '202009100001', jsonsimplify= TRUE)
#?mw
# parsing the nested lists
clim2 <- data.frame( lapply( d2$STATION$OBSERVATIONS, unlist) )
