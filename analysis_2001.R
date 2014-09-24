# load 2001 address data
geokodok <- read.table("../data/2001_geokodok_tabbed.TXT",
                       sep="\t",
                       stringsAsFactors=FALSE,
                       header=TRUE)
# load original dataset:
dat <- read.csv("../data/dat_w_add_2001.csv",
                stringsAsFactors=FALSE)
# load segregation functions:
source("indices.R")
source("segr_inds.R")

# merge geocodes:
fulldat <- merge(geokodok, dat,
                 by.x=c("KSH_kód", "Számlálókörzet"),
                 by.y=c("kshterkod","SZLOK"),
                 all.y=TRUE)
fulldat$area <- 1
# elkészítem a terület változót, hogy egyedi legyen:
fulldat$reg <- paste(fulldat$KSH_kód,
                     fulldat$Számlálókörzet,
                     sep="_")
# kiszedem az ismétrődőket, egyáltalán ne legyenek benne:
fulldat <- fulldat[!duplicated(fulldat$reg),]
# Let's roll
D <- segr_inds(fulldat,
               region="reg",
               county="KSH_kód",
               groups=names(fulldat)[grepl("status", names(fulldat))],
               indices=ls()[grepl("index", ls())],
               lng=as.name("EOV_Y"),
               lat=as.name("EOV_X"),
               mode="realistic",
               center_mode="GeoCentrum",
               area=as.name("area"),
               b=0.5)


# keep on testing
dat <- data.frame(minor=1:5)
dat$major <- dat$minor*10
dat$lng <- dat$minor*100
dat$lat <- dat$major*10
dat$area <- dat$minor*1000


