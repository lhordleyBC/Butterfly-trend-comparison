# Compare land cover of WCBS and tBMS squares 
# Using UKCEH Land Cover Map data 20215

library(ggplot2)
library(data.table)
library(dplyr)
library(giscoR)
library(tidyverse)
library(sf)
library(igr)
library(plotrix)
library(rnrfa)
library(gridExtra)
library(ggpubr)
library(raster)

## LAND COVER ANALYSIS
# Read in and filter data
bmsdata <- readRDS("Data/UKBMS/sp_weeklycount_1976-2024_zerofilled.rds")
# Filter to 2009 to 2024 
bmsdata <- bmsdata[YEAR >= 2009]
# Create survey ID column
bmsdata <- bmsdata %>% dplyr::mutate(SURVEY = ifelse(SITENO<50000, "tBMS", "WCBS"))
# filter by species
splist <- read.csv(file="Data/UKBMS/Species_list_GAM.csv", header=TRUE)
# filter bms data by species to get list of sites
bmsdata <- bmsdata[bmsdata$SPECIES %in% splist$SPECIES,]
length(unique(bmsdata$SITENO)) # 5,760
## read in site info
sitedata <- read.csv("Data/UKBMS/UKBMS_sites_location.csv", header=TRUE)
# filter this by sites included 
sitedata <- sitedata[sitedata$SITENO %in% bmsdata$SITENO,]
# remove year - not needed here
sitedata$YEAR <- NULL
sitedata <- unique(sitedata) # 5,760
length(unique(sitedata$MONAD)) # 5,168

length(unique(sitedata$SITENO[sitedata$SURVEY=="tBMS"])) # 3,558 tBMS sites from 2009-2024 for 23 species
length(unique(sitedata$SITENO[sitedata$SURVEY=="WCBS"])) # 2,202 WCBS squares from 2009-2024 for 23 species

## remove Channel Island sites - no land cover data
sitedata <- sitedata[!grepl("WV", sitedata$GRIDREF),]
sitedata <- sitedata[!grepl("WA", sitedata$GRIDREF),]
# 5,699 used in landscape and protected area analysis

#### GREAT BRITAIN
## Add in LCM 2015 data - GB ##
lcm_gb = raster::stack("Data/Land cover/LCM2015_GB_1km_percent_cover_aggregate_class.tif")
crs(lcm_gb) ## projection
res(lcm_gb)
#plot(lcm_gb)

# remove NI sites from sitedata
sitedata_gb <- sitedata %>% filter(nchar(sitedata$MONAD)==6)
# create unique monad dataset
monad_gb <- unique(sitedata_gb[,c("MONAD", "EAST_1K", "NORTH_1K", "SURVEY")]) # 5119
# extract LCM for each 1km square
bms_lcm_gb <- raster::extract(lcm_gb, monad_gb[,2:3], df=T)
bms_lcm_gb <- cbind(bms_lcm_gb, monad_gb)
bms_lcm_gb$ID <- NULL
bms_lcm_gb <- bms_lcm_gb[, c(11:14,1:10)]
bms_lcm_gb$tot <- rowSums(bms_lcm_gb[,5:14])
# remove sites with all zeros (i.e. no land cover data - 4 IoM sites)
bms_lcm_gb <- bms_lcm_gb[bms_lcm_gb$tot != 0, ]

#### NORTHERN IRELAND
## Add in LCM 2015 data 
lcm_ni = raster::stack("Data/Land cover/lcm2015_ni_1km_percent_cover_aggregate_class.tif")
crs(lcm_ni) ## projection
res(lcm_ni)
#plot(lcm_ni)

# keep NI sites from sitedata
sitedata_ni <- sitedata %>% filter(nchar(sitedata$MONAD)==5)
monad_ni <- unique(sitedata_ni[,c("MONAD", "EAST_1K", "NORTH_1K", "SURVEY")]) # 117 

bms_lcm_ni <- raster::extract(lcm_ni, monad_ni[,2:3], df=T)
bms_lcm_ni <- cbind(bms_lcm_ni, monad_ni)
bms_lcm_ni$ID <- NULL
bms_lcm_ni <- bms_lcm_ni[, c(11:14,1:10)]
bms_lcm_ni$tot <- rowSums(bms_lcm_ni[,5:14]) # there are no sites with tot = 0 (i.e. no land cover data)

## bind the two datasets together
# change column names first so they match 
colnames(bms_lcm_gb)[5:14] <- c("PERC_BROADLEAF", "PERC_CONIFEROUS", "PERC_ARABLE", "PERC_IMPROVED", "PERC_SEMINATURAL", "PERC_BOG",
                                "PERC_SALTWATER", "PERC_FRESHWATER", "PERC_COASTAL", "PERC_URBAN")
colnames(bms_lcm_ni)[5:14] <- c("PERC_BROADLEAF", "PERC_CONIFEROUS", "PERC_ARABLE", "PERC_IMPROVED", "PERC_SEMINATURAL", "PERC_BOG",
                                "PERC_SALTWATER", "PERC_FRESHWATER", "PERC_COASTAL", "PERC_URBAN")
bms_lcm <- rbind(bms_lcm_gb, bms_lcm_ni) 
bms_lcm <- na.omit(bms_lcm) # remove NAs - channel island and isle of man sites
length(unique(bms_lcm$MONAD)) # 5164 1km squares
# 5236 rows - some duplicates where there is a WCBS and tBMS transect within a 1km square
setDT(bms_lcm)[, .(count = uniqueN(MONAD)), by = SURVEY] # tBMS = 3034, WCBS = 2198

# All UKBMS 1km squares - i.e. remove duplcates where there is a WCBS and tBMS transect within a 1km squares
# NI is ok - there are 117 unique 1 km squares already
# So just repeat for GB
length(unique(monad_gb$MONAD)) # 5051 (but 5119 rows)
monad_gb$SURVEY <- NULL
monad_gb <- unique(monad_gb) # 5051 rows
ukbms_lcm_gb <- raster::extract(lcm_gb, monad_gb[,2:3], df=T)
ukbms_lcm_gb <- cbind(ukbms_lcm_gb, monad_gb)
ukbms_lcm_gb <- ukbms_lcm_gb[,-1]
ukbms_lcm_gb <- ukbms_lcm_gb[, c(11:13,1:10)]
ukbms_lcm_gb$tot <- rowSums(ukbms_lcm_gb[,4:13])
# remove sites with all zeros (i.e. no land cover data - 4 IoM sites)
ukbms_lcm_gb <- ukbms_lcm_gb[ukbms_lcm_gb$tot != 0, ]
colnames(ukbms_lcm_gb)[4:13] <- c("PERC_BROADLEAF", "PERC_CONIFEROUS", "PERC_ARABLE", "PERC_IMPROVED", "PERC_SEMINATURAL", "PERC_BOG",
                                "PERC_SALTWATER", "PERC_FRESHWATER", "PERC_COASTAL", "PERC_URBAN")
bms_lcm_ni$SURVEY <- NULL
ukbms_lcm <- rbind(ukbms_lcm_gb, bms_lcm_ni) 
ukbms_lcm <- na.omit(ukbms_lcm) # remove NAs - channel island and isle of man sites

## All UK (GB + NI separately) LCM 

lcm_gb_df <- raster::as.data.frame(lcm_gb, xy=TRUE)
lcm_gb_df <- lcm_gb_df[ rowSums(lcm_gb_df[,3:12]) > 0, ]
# 242,677 1km squares 

lcm_ni_df <- raster::as.data.frame(lcm_ni, xy=TRUE)
lcm_ni_df <- lcm_ni_df[ rowSums(lcm_ni_df[,3:12]) > 0, ]
# 14,867 1km squares 

colnames(lcm_gb_df)[3:12] <- c("PERC_BROADLEAF", "PERC_CONIFEROUS", "PERC_ARABLE", "PERC_IMPROVED", "PERC_SEMINATURAL", "PERC_BOG",
                                "PERC_SALTWATER", "PERC_FRESHWATER", "PERC_COASTAL", "PERC_URBAN")
colnames(lcm_ni_df)[3:12] <- c("PERC_BROADLEAF", "PERC_CONIFEROUS", "PERC_ARABLE", "PERC_IMPROVED", "PERC_SEMINATURAL", "PERC_BOG",
                                "PERC_SALTWATER", "PERC_FRESHWATER", "PERC_COASTAL", "PERC_URBAN")

# bind together
lcm_gb_ni <- rbind(lcm_gb_df, lcm_ni_df)
# save the coordinates of the squares to use for protected area too


## Calculate mean and 2*SE for UK, tBMS and WCBS separately 
# remove land cover types not interested in first
bms_lcm$tot <- NULL
bms_lcm <- subset(bms_lcm, select=-c(PERC_SALTWATER,PERC_FRESHWATER))
bms_lcm_l <- melt(setDT(bms_lcm), id.vars = c("SURVEY", "MONAD", "EAST_1K", "NORTH_1K"), variable.name = "LAND_COVER") # change to long format
sum_bms_lcm <- bms_lcm_l %>% group_by(SURVEY, LAND_COVER) %>% summarise(mean=mean(value), se=2*(plotrix::std.error(value)))

ukbms_lcm$tot <- NULL
ukbms_lcm <- subset(ukbms_lcm, select=-c(PERC_SALTWATER,PERC_FRESHWATER))
ukbms_lcm$SURVEY <- "Combined"
ukbms_lcm_l <- melt(setDT(ukbms_lcm), id.vars = c("SURVEY", "MONAD", "EAST_1K", "NORTH_1K"), variable.name = "LAND_COVER") # change to long format
sum_ukbms_lcm <- ukbms_lcm_l %>% group_by(SURVEY, LAND_COVER) %>% summarise(mean=mean(value), se=2*(plotrix::std.error(value)))

lcm_gb_ni <- subset(lcm_gb_ni, select=-c(PERC_SALTWATER,PERC_FRESHWATER))
lcm_gb_ni$SURVEY <- "UK"
lcm_gb_ni_l <- melt(setDT(lcm_gb_ni), id.vars = c("x", "y", "SURVEY"), variable.name = "LAND_COVER") # change to long format
sum_gb_ni_lcm <- lcm_gb_ni_l %>% group_by(SURVEY, LAND_COVER) %>% summarise(mean=mean(value), se=2*(plotrix::std.error(value)))

sum_lcm <- rbind(sum_bms_lcm, sum_ukbms_lcm, sum_gb_ni_lcm)

## plot result

sum_lcm$SURVEY <- factor(sum_lcm$SURVEY, levels = c("tBMS", "WCBS", "Combined", "UK"))
land_cover_labs <- c("PERC_BROADLEAF"="Broadleaf \nwoodland", "PERC_CONIFEROUS"="Coniferous \nwoodland", "PERC_ARABLE"="Arable",
                    "PERC_IMPROVED"="Improved \ngrassland", "PERC_SEMINATURAL"="Semi-natural \ngrassland",
                    "PERC_BOG"="Mountain, heath \nand bog", "PERC_COASTAL"="Coastal", "PERC_URBAN"="Built-up areas \nand gardens")

lcm <- ggplot(data=sum_lcm)+
  facet_wrap(~LAND_COVER, scales="free", labeller=as_labeller(land_cover_labs))+
  geom_point(mapping=aes(x=SURVEY, y=mean, colour=SURVEY), size=4)+
  geom_linerange(aes(x=SURVEY, ymin=mean-se, ymax=mean+se, colour=SURVEY), lwd=1)+
  ylab("Percent coverage")+
  xlab("Survey")+
  theme_minimal()+
  theme(text=element_text(size=30),
        legend.title=element_blank(),
        legend.key.size = unit(1.5, "cm"),
        legend.key = element_rect(color = NA, fill = NA),
        axis.text.x=element_text(angle=45, size=18, vjust=0.5))
lcm
# save plot
ggsave(lcm, file="Output/Figures/Figure2.png", height=12, width=16)

# save table of mean + SE
write.csv(sum_lcm, file="Output/Landscape analysis/LCM_mean_SE.csv", row.names=FALSE)
sum_lcm <- read.csv("Output/Landscape analysis/LCM_mean_SE.csv", header=TRUE)






#############################
## same land cover analyses as above, but with 21 target classes (rather than 10 aggregate classes)

## LAND COVER ANALYSIS
# Read in and filter data
bmsdata <- readRDS("Data/UKBMS/sp_weeklycount_1976-2024_zerofilled.rds")
# Filter to 2009 to 2024 
bmsdata <- bmsdata[YEAR >= 2009]
# Create survey ID column
bmsdata <- bmsdata %>% dplyr::mutate(SURVEY = ifelse(SITENO<50000, "tBMS", "WCBS"))
# filter by species
splist <- read.csv(file="Data/UKBMS/Species_list_GAM.csv", header=TRUE)
# filter bms data by species to get list of sites
bmsdata <- bmsdata[bmsdata$SPECIES %in% splist$SPECIES,]
length(unique(bmsdata$SITENO)) # 5,760
## read in site info
sitedata <- read.csv("Data/UKBMS/UKBMS_sites_location.csv", header=TRUE)
# filter this by sites included 
sitedata <- sitedata[sitedata$SITENO %in% bmsdata$SITENO,]
# remove year - not needed here
sitedata$YEAR <- NULL
sitedata <- unique(sitedata) # 5,760
## remove Channel Island sites - no land cover data
sitedata <- sitedata[!grepl("WV", sitedata$GRIDREF),]
sitedata <- sitedata[!grepl("WA", sitedata$GRIDREF),]
# 5,699 used in landscape and protected area analysis
length(unique(sitedata$MONAD)) # 5,168

#### GREAT BRITAIN
## Add in LCM 2015 data - GB ##
library(raster)
lcm_gb = stack("Data/Land cover/LCM2015_GB_1km_percent_cover_target_class.tif")
crs(lcm_gb) ## projection
res(lcm_gb)
#plot(lcm_gb)

# remove NI sites from sitedata
sitedata_gb <- sitedata %>% filter(nchar(sitedata$MONAD)==6)
# create unique monad dataset
monad_gb <- unique(sitedata_gb[,c("MONAD", "EAST_1K", "NORTH_1K", "SURVEY")]) # 5119
# extract LCM for each 1km square
bms_lcm_gb <- raster::extract(lcm_gb, monad_gb[,2:3], df=T)
bms_lcm_gb <- cbind(bms_lcm_gb, monad_gb)
bms_lcm_gb <- bms_lcm_gb[,-1]
bms_lcm_gb <- bms_lcm_gb[, c(22:25,1:21)]
bms_lcm_gb$tot <- rowSums(bms_lcm_gb[,5:25])
# remove sites with all zeros (i.e. no land cover data - 4 IoM sites)
bms_lcm_gb <- bms_lcm_gb[bms_lcm_gb$tot != 0, ]
bms_lcm_gb$tot <- NULL

#### NORTHERN IRELAND
## Add in LCM 2015 data 
lcm_ni = stack("Data/Land cover/lcm2015_ni_1km_percent_cover_target_class.tif")
crs(lcm_ni) ## projection
res(lcm_ni)
#plot(lcm_ni)

# keep NI sites from sitedata
sitedata_ni <- sitedata %>% filter(nchar(sitedata$MONAD)==5)
monad_ni <- unique(sitedata_ni[,c("MONAD", "EAST_1K", "NORTH_1K", "SURVEY")]) # 117

bms_lcm_ni <- raster::extract(lcm_ni, monad_ni[,2:3], df=T)
bms_lcm_ni <- cbind(bms_lcm_ni, monad_ni)
bms_lcm_ni <- bms_lcm_ni[,-1]
bms_lcm_ni <- bms_lcm_ni[, c(22:25,1:21)]
bms_lcm_ni$tot <- rowSums(bms_lcm_ni[,5:25])
bms_lcm_ni$tot <- NULL

## bind the two datasets together
# change column names first so they match 
colnames(bms_lcm_gb)[5:25] <- c("PERC_BROADLEAF", "PERC_CONIFEROUS", "PERC_ARABLE", "PERC_IMPROVED", "PERC_NEUTRAL", "PERC_CALCAREOUS",
                                "PERC_ACID", "PERC_FEN_MARSH_SWAMP", "PERC_HEATHER", "PERC_HEATHER_GRASSLAND", "PERC_BOG", "PERC_INLAND_ROCK",
                                "PERC_SALTWATER", "PERC_FRESHWATER", "PERC_SUPRA_LITTORAL_ROCK", "PERC_SUPRA_LITTORAL_SEDIMENT", "PERC_LITTORAL_ROCK",
                                "PERC_LITTORAL_SEDIMENT", "PERC_SALTMARSH", "PERC_URBAN", "PERC_SUBURBAN")
colnames(bms_lcm_ni)[5:25] <- c("PERC_BROADLEAF", "PERC_CONIFEROUS", "PERC_ARABLE", "PERC_IMPROVED", "PERC_NEUTRAL", "PERC_CALCAREOUS",
                                "PERC_ACID", "PERC_FEN_MARSH_SWAMP", "PERC_HEATHER", "PERC_HEATHER_GRASSLAND", "PERC_BOG", "PERC_INLAND_ROCK",
                                "PERC_SALTWATER", "PERC_FRESHWATER", "PERC_SUPRA_LITTORAL_ROCK", "PERC_SUPRA_LITTORAL_SEDIMENT", "PERC_LITTORAL_ROCK",
                                "PERC_LITTORAL_SEDIMENT", "PERC_SALTMARSH", "PERC_URBAN", "PERC_SUBURBAN")
## NOTE: BROADLEAF, CONIFER, ARABLE, AND IMPROVED GRASSLAND CATEGORIES ARE THE SAME
# not interested in coastal - just look at urban, mountain/heath/bog and semi-natural grassland
bms_lcm <- rbind(bms_lcm_gb, bms_lcm_ni) 
bms_lcm <- na.omit(bms_lcm) # remove NAs - channel island and isle of man sites
length(unique(bms_lcm$MONAD)) # 5164 1km squares

# All UKBMS 1km squares - i.e. remove duplcates where there is a WCBS and tBMS transect within a 1km squares
# NI is ok - there are 115 unique 1 km squares alread
# So just repeat for GB
length(unique(monad_gb$MONAD)) # 5051 (but 5119 rows)
monad_gb$SURVEY <- NULL
monad_gb <- unique(monad_gb) # 5051 rows
ukbms_lcm_gb <- raster::extract(lcm_gb, monad_gb[,2:3], df=T)
ukbms_lcm_gb <- cbind(ukbms_lcm_gb, monad_gb)
ukbms_lcm_gb <- ukbms_lcm_gb[,-1]
ukbms_lcm_gb <- ukbms_lcm_gb[, c(22:24,1:21)]
ukbms_lcm_gb$tot <- rowSums(ukbms_lcm_gb[,4:24])
# remove sites with all zeros (i.e. no land cover data - 4 IoM sites)
ukbms_lcm_gb <- ukbms_lcm_gb[ukbms_lcm_gb$tot != 0, ]
ukbms_lcm_gb$tot <- NULL

colnames(ukbms_lcm_gb)[4:24] <- c("PERC_BROADLEAF", "PERC_CONIFEROUS", "PERC_ARABLE", "PERC_IMPROVED", "PERC_NEUTRAL", "PERC_CALCAREOUS",
                                "PERC_ACID", "PERC_FEN_MARSH_SWAMP", "PERC_HEATHER", "PERC_HEATHER_GRASSLAND", "PERC_BOG", "PERC_INLAND_ROCK",
                                "PERC_SALTWATER", "PERC_FRESHWATER", "PERC_SUPRA_LITTORAL_ROCK", "PERC_SUPRA_LITTORAL_SEDIMENT", "PERC_LITTORAL_ROCK",
                                "PERC_LITTORAL_SEDIMENT", "PERC_SALTMARSH", "PERC_URBAN", "PERC_SUBURBAN")
bms_lcm_ni$SURVEY <- NULL
ukbms_lcm <- rbind(ukbms_lcm_gb, bms_lcm_ni) 
ukbms_lcm <- na.omit(ukbms_lcm) # remove NAs - channel island and isle of man sites


## All UK (GB + NI separately) LCM 

lcm_gb_df <- raster::as.data.frame(lcm_gb, xy=TRUE)
lcm_gb_df <- lcm_gb_df[ rowSums(lcm_gb_df[,3:23]) > 0, ]
# 242,677 1km squares 

lcm_ni_df <- raster::as.data.frame(lcm_ni, xy=TRUE)
lcm_ni_df <- lcm_ni_df[ rowSums(lcm_ni_df[,3:23]) > 0, ]
# 14,867 1km squares 

colnames(lcm_gb_df)[3:23] <- c("PERC_BROADLEAF", "PERC_CONIFEROUS", "PERC_ARABLE", "PERC_IMPROVED", "PERC_NEUTRAL", "PERC_CALCAREOUS",
                               "PERC_ACID", "PERC_FEN_MARSH_SWAMP", "PERC_HEATHER", "PERC_HEATHER_GRASSLAND", "PERC_BOG", "PERC_INLAND_ROCK",
                               "PERC_SALTWATER", "PERC_FRESHWATER", "PERC_SUPRA_LITTORAL_ROCK", "PERC_SUPRA_LITTORAL_SEDIMENT", "PERC_LITTORAL_ROCK",
                               "PERC_LITTORAL_SEDIMENT", "PERC_SALTMARSH", "PERC_URBAN", "PERC_SUBURBAN")
colnames(lcm_ni_df)[3:23] <- c("PERC_BROADLEAF", "PERC_CONIFEROUS", "PERC_ARABLE", "PERC_IMPROVED", "PERC_NEUTRAL", "PERC_CALCAREOUS",
                               "PERC_ACID", "PERC_FEN_MARSH_SWAMP", "PERC_HEATHER", "PERC_HEATHER_GRASSLAND", "PERC_BOG", "PERC_INLAND_ROCK",
                               "PERC_SALTWATER", "PERC_FRESHWATER", "PERC_SUPRA_LITTORAL_ROCK", "PERC_SUPRA_LITTORAL_SEDIMENT", "PERC_LITTORAL_ROCK",
                               "PERC_LITTORAL_SEDIMENT", "PERC_SALTMARSH", "PERC_URBAN", "PERC_SUBURBAN")

# bind together
lcm_gb_ni <- rbind(lcm_gb_df, lcm_ni_df)
# save the coordinates of the squares to use for protected area too


## Calculate mean and 2*SE for UK, tBMS and WCBS separately 
# remove land cover types not interested in first
bms_lcm <- bms_lcm[ , !names(bms_lcm) %in% c("PERC_BROADLEAF", "PERC_CONIFEROUS", "PERC_ARABLE", "PERC_IMPROVED", "PERC_SALTWATER", 
                                             "PERC_FRESHWATER", "PERC_SUPRA_LITTORAL_ROCK", "PERC_SUPRA_LITTORAL_SEDIMENT", "PERC_LITTORAL_ROCK",
                                             "PERC_LITTORAL_SEDIMENT", "PERC_SALTMARSH")] 
bms_lcm_l <- melt(setDT(bms_lcm), id.vars = c("SURVEY", "MONAD", "EAST_1K", "NORTH_1K"), variable.name = "LAND_COVER") # change to long format
sum_bms_lcm <- bms_lcm_l %>% group_by(SURVEY, LAND_COVER) %>% summarise(mean=mean(value), se=2*(plotrix::std.error(value)))

ukbms_lcm <- ukbms_lcm[ , !names(ukbms_lcm) %in% c("PERC_BROADLEAF", "PERC_CONIFEROUS", "PERC_ARABLE", "PERC_IMPROVED", "PERC_SALTWATER", 
                                             "PERC_FRESHWATER", "PERC_SUPRA_LITTORAL_ROCK", "PERC_SUPRA_LITTORAL_SEDIMENT", "PERC_LITTORAL_ROCK",
                                             "PERC_LITTORAL_SEDIMENT", "PERC_SALTMARSH")] 
ukbms_lcm$SURVEY <- "Combined"
ukbms_lcm_l <- melt(setDT(ukbms_lcm), id.vars = c("SURVEY", "MONAD", "EAST_1K", "NORTH_1K"), variable.name = "LAND_COVER") # change to long format
sum_ukbms_lcm <- ukbms_lcm_l %>% group_by(SURVEY, LAND_COVER) %>% summarise(mean=mean(value), se=2*(plotrix::std.error(value)))


lcm_gb_ni <- lcm_gb_ni[ , !names(lcm_gb_ni) %in% c("PERC_BROADLEAF", "PERC_CONIFEROUS", "PERC_ARABLE", "PERC_IMPROVED", "PERC_SALTWATER", 
                                                   "PERC_FRESHWATER", "PERC_SUPRA_LITTORAL_ROCK", "PERC_SUPRA_LITTORAL_SEDIMENT", "PERC_LITTORAL_ROCK",
                                                   "PERC_LITTORAL_SEDIMENT", "PERC_SALTMARSH")] 
lcm_gb_ni$SURVEY <- "UK"
lcm_gb_ni_l <- melt(setDT(lcm_gb_ni), id.vars = c("x", "y", "SURVEY"), variable.name = "LAND_COVER") # change to long format
sum_gb_ni_lcm <- lcm_gb_ni_l %>% group_by(SURVEY, LAND_COVER) %>% summarise(mean=mean(value), se=2*(plotrix::std.error(value)))

sum_lcm <- rbind(sum_bms_lcm, sum_ukbms_lcm, sum_gb_ni_lcm)

## plot result

sum_lcm$SURVEY <- factor(sum_lcm$SURVEY, levels = c("tBMS", "WCBS", "Combined", "UK"))
land_cover_labs <- c("PERC_NEUTRAL"="Neutral \ngrassland", "PERC_CALCAREOUS"="Calcareous \ngrassland", "PERC_ACID"="Acid \ngrassland",
                     "PERC_FEN_MARSH_SWAMP"="Fen, marsh \nand swamp", "PERC_HEATHER"="Heather", "PERC_HEATHER_GRASSLAND"="Heather \ngrassland",
                     "PERC_BOG"="Bog", "PERC_INLAND_ROCK"="Inland rock", "PERC_URBAN"="Urban", "PERC_SUBURBAN"="Suburban")
sum_lcm <- sum_lcm[sum_lcm$LAND_COVER %in% c("PERC_NEUTRAL", "PERC_CALCAREOUS", "PERC_ACID", "PERC_FEN_MARSH_SWAMP"), ]

lcm <- ggplot(data=sum_lcm)+
  facet_wrap(~LAND_COVER, scales="free", labeller=as_labeller(land_cover_labs))+
  geom_point(mapping=aes(x=SURVEY, y=mean, colour=SURVEY), size=4)+
  geom_linerange(aes(x=SURVEY, ymin=mean-se, ymax=mean+se, colour=SURVEY), lwd=1)+
  ylab("Percent coverage")+
  xlab("Survey")+
  theme_minimal()+
  theme(text=element_text(size=24),
        legend.title=element_blank(),
        legend.key.size = unit(1.5, "cm"),
        legend.key = element_rect(color = NA, fill = NA),
        axis.text.x=element_text(angle=45, size=18, vjust=0.5))
lcm
# save plot
ggsave(lcm, file="Output/Figures/FigureS2.png", height=8, width=10)

# save table of mean + SE
write.csv(sum_lcm, file="Output/Landscape analysis/LCM_mean_SE_grassland_bog_urban.csv", row.names=FALSE)

