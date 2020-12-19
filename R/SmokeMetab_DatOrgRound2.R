dat2 <- read.csv("https://www.dropbox.com/s/eochnr2z8v3qsip/nwis.11128500.csv?dl=1")
summary(dat2) # SANTA YNEZ R A SOLVANG CA

# near (34.58498706, -120.1445926)
mw(service = 'metadata', complete=1, state='CA', county='Santa Barbara County') # near (33.4739175,-117.1422536)
# KPSP
?mw

# https://download-1.s3-us-west-2.amazonaws.com/8r2bwp41xj80clxn3tpsfg99lgtxu0plvd1gdyjzp244f2y12169vsu3y77nwodexo3nju/KIZA.2020-09-10.csv
climdat2 <- read.csv("/SM_mesonetDat/KIZA.2020-09-10.csv")
summary(climdat2)

getwd()