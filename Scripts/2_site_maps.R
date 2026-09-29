## create map of sites used in analysis (including creating 10km grid references)
## calculate proportion of WCBS vs tBMS in each country and region

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

# Read in UKBMS data
bmsdata <- readRDS("Data/UKBMS/sp_weeklycount_1976-2024_zerofilled.rds")
# Filter to 2009 to 2024 
bmsdata <- bmsdata[YEAR >= 2009]
# Create survey ID column
bmsdata <- bmsdata %>% dplyr::mutate(SURVEY = ifelse(SITENO<50000, "tBMS", "WCBS"))

## Add in UKBMS location data for each site
# take columns of interest
sitedata <- unique(bmsdata[,c("SITENO", "YEAR", "SURVEY")]) 
length(unique(sitedata$SITENO)) # 5,833 sites
# read in site location data
site_info <- read.csv("Data/UKBMS/site_location_data_2024.csv", header=TRUE)
sitedata <- merge(x = sitedata, y = site_info[ , c("SITENO", "GRIDREF", "EAST", "NORTH", "LONGITUDE", "LATITUDE")], by = "SITENO", all.x=TRUE)
# 10 WCBS sites without any location data - read a separate file in to get the lat and lon of these
channel_isle_info <- read.csv("Data/UKBMS/Channel_Isles_location_data.csv")
channel_isle_info <- channel_isle_info %>% dplyr::select(-Site.Name)
names(channel_isle_info) <- c("SITENO", "LATITUDE", "LONGITUDE") 
sitedata_nas <- sitedata %>% filter(if_any(everything(), is.na))
sitedata_nas <- dplyr::select(sitedata_nas, -c('LATITUDE', 'LONGITUDE'))
sitedata_nas2 <- merge(sitedata_nas, channel_isle_info, by="SITENO", all.x=TRUE)
sitedata_nas2$GRIDREF <- c("WV38", "WV27", "WV37", "WV38", "WV27", "WV27", "WV27", "WV55", "WV54", "WV54") # obtained these online
sitedata <- na.omit(sitedata)
sitedata <- dplyr::bind_rows(sitedata, sitedata_nas2)
# this is waaaay too long but it works
length(unique(sitedata$SITENO)) # 5,833 sites

## Create 10k grid references
sitedata <- sitedata %>% mutate(HECTAD=case_when(
  nchar(GRIDREF)==4 ~ GRIDREF,
  nchar(GRIDREF)==5 ~ paste(substr(GRIDREF, 1,2), substr(GRIDREF, 4,4), sep=""),
  nchar(GRIDREF)==6 ~ paste(substr(GRIDREF, 1,3), substr(GRIDREF, 5,5), sep=""),
  nchar(GRIDREF)==7 ~ paste(substr(GRIDREF, 1,2), substr(GRIDREF, 5,5), sep=""),
  nchar(GRIDREF)==8 ~ paste(substr(GRIDREF, 1,3), substr(GRIDREF, 6,6), sep=""),
  nchar(GRIDREF)==10 ~ paste(substr(GRIDREF, 1,3), substr(GRIDREF, 7,7), sep=""),
  nchar(GRIDREF)==12 ~ paste(substr(GRIDREF, 1,3), substr(GRIDREF, 8,8), sep=""),
))

## Create 1k grid references (if already 1k, stay the same, if longer/finer scale, make it shorter)
sitedata <- sitedata %>% mutate(MONAD=case_when(
  nchar(GRIDREF)==4 ~ GRIDREF, # these are some channel island gridrefs - don't need to be finer scale
  nchar(GRIDREF)==5 ~ GRIDREF,
  nchar(GRIDREF)==6 ~ GRIDREF,
  nchar(GRIDREF)==7 ~ paste(substr(GRIDREF, 1,3), substr(GRIDREF, 5,6), sep=""),
  nchar(GRIDREF)==8 ~ paste(substr(GRIDREF, 1,4), substr(GRIDREF, 6,7), sep=""),
  nchar(GRIDREF)==10 ~ paste(substr(GRIDREF, 1,4), substr(GRIDREF, 7,8), sep=""),
  nchar(GRIDREF)==12 ~ paste(substr(GRIDREF, 1,4), substr(GRIDREF, 8,9), sep=""),
))
length(unique(sitedata$MONAD)) # 5,270 unique 1 km squares 

# first take Channel Island grid refs separately
# these don't convert using osg_parse, so do it manually
# only need the 10k lat/lon to plot map below
# as these and IoM sites are removed for landscape + protected area analysis
all_sites_ci <-  sitedata[grep('^W', sitedata$HECTAD),]
all_sites_ci$EAST_10K <- NA
all_sites_ci$NORTH_10K <- NA
all_sites_ci$EAST_1K <- NA
all_sites_ci$NORTH_1K <- NA
ci_lat_lon <- read.csv(file="Data/Channel_islands_hectad_lat_lon.csv", header=TRUE)
all_sites_ci <- merge(all_sites_ci, ci_lat_lon, by="HECTAD")
all_sites_ci$lon_1k <- NA
all_sites_ci$lat_1k <- NA

###### GB #######
all_sites_gb <-  sitedata[!grep('^W', sitedata$HECTAD),] # remove channel island sites
all_sites_gb <- all_sites_gb %>% filter(nchar(all_sites_gb$MONAD)==6) # keep only 6 figure GRs (i.e. remove NI sites)
length(unique(all_sites_gb$MONAD)) # 5,093 unique 1 km squares in GB
# change 10 km grid refs to easting/northing 
gb_east_north <- as.data.frame(osg_parse(grid_refs=all_sites_gb$HECTAD))
# move easting northing to middle of hectad
gb_east_north$easting <- gb_east_north$easting + 5000
gb_east_north$northing <- gb_east_north$northing + 5000
colnames(gb_east_north) <- c("EAST_10K", "NORTH_10K")
all_sites_gb <- cbind(all_sites_gb, gb_east_north)
# change 1 km grid refs to easting/northing 
gb_east_north <- as.data.frame(osg_parse(grid_refs=all_sites_gb$MONAD))
# move easting northing to middle of hectad
gb_east_north$easting <- gb_east_north$easting + 500
gb_east_north$northing <- gb_east_north$northing + 500
colnames(gb_east_north) <- c("EAST_1K", "NORTH_1K")
all_sites_gb <- cbind(all_sites_gb, gb_east_north)
# change 10 km easting/northing to lat/ lon
gb_lat_lon <- all_sites_gb %>%
  st_as_sf(coords = c("EAST_10K", "NORTH_10K"), crs = 27700) %>%
  st_transform(4326) %>%
  st_coordinates() %>%
  as_tibble()
colnames(gb_lat_lon) <- c("lon_10k","lat_10k")
all_sites_gb <- cbind(all_sites_gb, gb_lat_lon)
# change 1 km easting/northing to lat/ lon
gb_lat_lon <- all_sites_gb %>%
  st_as_sf(coords = c("EAST_1K", "NORTH_1K"), crs = 27700) %>%
  st_transform(4326) %>%
  st_coordinates() %>%
  as_tibble()
colnames(gb_lat_lon) <- c("lon_1k","lat_1k")
all_sites_gb <- cbind(all_sites_gb, gb_lat_lon)

###### Nothern Ireland ###### 
all_sites_ni <- sitedata %>% filter(nchar(sitedata$MONAD)==5)
# change 10 km grid refs to easting/northing 
ni_east_north <- as.data.frame(igr_to_ig(all_sites_ni$HECTAD, centroids = TRUE))
colnames(ni_east_north) <- c("EAST_10K", "NORTH_10K")
all_sites_ni <- cbind(all_sites_ni, ni_east_north)
ni_east_north <- as.data.frame(igr_to_ig(all_sites_ni$MONAD, centroids = TRUE))
colnames(ni_east_north) <- c("EAST_1K", "NORTH_1K")
all_sites_ni <- cbind(all_sites_ni, ni_east_north)
## 10 km grid reference to lat/lon
ni_hectads <- data.frame(igr=all_sites_ni$HECTAD)
ni_lat_lon <- st_igr_as_sf(ni_hectads, "igr", crs = 4326, centroids = TRUE)
ni_lat_lon <- ni_lat_lon %>%
  dplyr::mutate(lon_10k = sf::st_coordinates(.)[,1],
                lat_10k = sf::st_coordinates(.)[,2])
ni_lat_lon <- st_drop_geometry(ni_lat_lon)
ni_lat_lon <- ni_lat_lon[,-1]
all_sites_ni <- cbind(all_sites_ni, ni_lat_lon)
## 1 km grid reference to lat/lon
ni_monads <- data.frame(igr=all_sites_ni$MONAD)
ni_lat_lon <- st_igr_as_sf(ni_hectads, "igr", crs = 4326, centroids = TRUE)
ni_lat_lon <- ni_lat_lon %>%
  dplyr::mutate(lon_1k = sf::st_coordinates(.)[,1],
                lat_1k = sf::st_coordinates(.)[,2])
ni_lat_lon <- st_drop_geometry(ni_lat_lon)
ni_lat_lon <- ni_lat_lon[,-1]
all_sites_ni <- cbind(all_sites_ni, ni_lat_lon)

## Add datasets together
all_sites_final <- rbind(all_sites_gb, all_sites_ni, all_sites_ci) # 5,833 sites
length(unique(all_sites_final$MONAD)) # 5,270
length(unique(all_sites_final$SITENO[all_sites_final$SURVEY=="tBMS"])) # 3,624 tBMS sites from 2009-2024 for ALL species
length(unique(all_sites_final$SITENO[all_sites_final$SURVEY=="WCBS"])) # 2,209 WCBS squares from 2009-2024 for ALL species
# save file
write.csv(all_sites_final, file="Data/UKBMS/UKBMS_sites_location.csv", row.names=FALSE)

######################################################################################################
######################################################################################################

# Plot map of sites used in GAM analysis - need to filter by species
# Read in UKBMS data
bmsdata <- readRDS("Data/UKBMS/sp_weeklycount_1976-2024_zerofilled.rds")
# Filter to 2009 to 2024 
bmsdata <- bmsdata[YEAR >= 2009]
# Create survey ID column
bmsdata <- bmsdata %>% dplyr::mutate(SURVEY = ifelse(SITENO<50000, "tBMS", "WCBS"))
# Filter to species which are recorded at greater than 25 WCBS squares per year between 2009 and 2024
# Only based on WCBS as there are fewer sites recorded for WCBS
# These are the species used in the GAM model
sp_site_count <- bmsdata[COUNT>0, .(no_sites=uniqueN(SITENO)), by=c("COMMON_NAME","SPECIES", "YEAR", "SURVEY")]
sp_site_count2 <- sp_site_count %>% group_by(COMMON_NAME) %>% filter(no_sites > 25 & SURVEY=="WCBS") %>% filter(all(c(2009:2024) %in% YEAR))
length(unique(sp_site_count2$COMMON_NAME)) # 23 species
splist <- unique(sp_site_count2[c("COMMON_NAME", "SPECIES")])
write.csv(splist, file="Data/UKBMS/Species_list_GAM.csv", row.names=FALSE)
# filter bms data by species to get list of sites
bmsdata <- bmsdata[bmsdata$SPECIES %in% splist$SPECIES,]
length(unique(bmsdata$SITENO)) # 5,760

## read in site info
all_sites_final <- read.csv("Data/UKBMS/UKBMS_sites_location.csv", header=TRUE)
# filter this by sites included 
all_sites_final <- all_sites_final[all_sites_final$SITENO %in% bmsdata$SITENO,]
length(unique(all_sites_final$SITENO)) # 5,760

# remove year - not needed here
all_sites_final$YEAR <- NULL
all_sites_final <- unique(all_sites_final) # 5,760

## Count the number of sites within a hectad to plot on a map
hectad_data <- all_sites_final[,c("HECTAD", "SURVEY", "lon_10k", "lat_10k")]
hectad_count <- hectad_data %>% group_by(HECTAD, SURVEY, lon_10k, lat_10k) %>% dplyr::summarise(count=n())
# create groups for plotting
hectad_count <- hectad_count %>% mutate(group = case_when(count==1 ~ '1',
                                                          count==2 ~ '2',
                                                          count==3 ~ '3',
                                                          count==4 ~ '4',
                                                          count==5 ~ '5',
                                                          count==6 ~ '6',
                                                          count>=7 ~ '7 or more')) %>%
  group_by(group) 
hectad_count <- as.data.frame(hectad_count)
# tBMS
names <- c("GBR", "IMN", "GGY", "JEY")
UK <- gisco_get_countries(country = names, resolution = 1)
tbms_sites2 <- ggplot() + 
  geom_sf(data = UK, fill = "grey", alpha = 0.3) + 
  geom_point(data = hectad_count[hectad_count$SURVEY=="tBMS",], 
             aes(x = as.numeric(lon_10k), 
                 y = as.numeric(lat_10k), colour=group), shape=20, size=2) + 
  scale_color_viridis_d() +
  guides(col=guide_legend(title='Number of sites \nper 10km square',override.aes = list(size=10))) +
  theme_void() +
  theme(title = element_text(size = 20), legend.text = element_text(size=16), legend.position="bottom")
tbms_sites2
# WCBS
wcbs_sites2 <- ggplot() + 
  geom_sf(data = UK, fill = "grey", alpha = 0.3) + 
  geom_point(data = hectad_count[hectad_count$SURVEY=="WCBS",], 
             aes(x = as.numeric(lon_10k), 
                 y = as.numeric(lat_10k), colour=group), shape=20, size=2) + 
  scale_color_viridis_d() +
  guides(col=guide_legend(title='Number of sites \nper 10km square',override.aes = list(size=10))) +
  theme_void() +
  theme(title = element_text(size = 20), legend.text = element_text(size=16), legend.position="bottom")
wcbs_sites2

site_maps <- ggarrange(tbms_sites2, wcbs_sites2, common.legend=TRUE, legend="bottom", labels=c("(a)", "(b)"), 
                       font.label=list(color="black",size=16))
site_maps

ggsave(site_maps, file="Output/Figures/Figure1.png", height=22, width=20, units="cm")


######################################################################################################
######################################################################################################

### count the number of tBMS and WCBS sites in each country and each region in England based on the 23 species
# Read in UKBMS data
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
all_sites_final <- read.csv("Data/UKBMS/UKBMS_sites_location.csv", header=TRUE)
# filter this by sites included 
sitedata <- all_sites_final[all_sites_final$SITENO %in% bmsdata$SITENO,]
length(unique(sitedata$SITENO)) # 5,760

# first plot a graph showing the number of tBMS and WCBS sites over time
## count number of sites per year for each survey
no_sites <- sitedata %>% group_by(SURVEY, YEAR) %>% summarise(no_sites=n())
number_sites <- ggplot(data=no_sites, aes(x=YEAR, y=no_sites, fill=SURVEY))+
  geom_bar(stat="identity")+
  labs(y="Number of sites", x="Year")+
  theme_bw()+
  theme(legend.title=element_blank())
number_sites
ggsave(number_sites, file="Output/Figures/Number_sites.png", height=8, width=10, units="cm")
# WCBS remained more stable over time - tBMS has increased


### now get data prepared by county boundaries to plot a graph of proportion of WCBS sites per county
## remove Channel Island sites - they don't fall within the country boundaries
sitedata <- sitedata[!grepl("WV", sitedata$GRIDREF),]
sitedata <- sitedata[!grepl("WA", sitedata$GRIDREF),]

# first, use the country shapefile to categorise each grid-reference into a country
countries <- st_read("Data/Boundaries/CTRY_DEC_2024_UK_BFC.shp")
# Ensure CRS match
sitedata <- na.omit(sitedata)
sitedata <- st_as_sf(sitedata, coords=c("EAST", "NORTH"))
sitedata <- sitedata %>% st_set_crs(st_crs(countries))

sitedata_countries <- st_join(sitedata, countries, join = st_intersects)
# some of these have NAs because the grid reference is technically in the sea
# remove these and run through st_nearest_feature instead
sitedata_nas <- sitedata_countries[is.na(sitedata_countries$CTRY24NM),]
sitedata_nas <- sitedata_nas[,c(1:6,15)]
sitedata_countries_nas <- st_join(sitedata_nas, countries, join = st_nearest_feature)
sitedata_countries <- na.omit(sitedata_countries)
sitedata_countries <- rbind(sitedata_countries, sitedata_countries_nas)
# select columns of interest
sitedata_countries <- as.data.frame(sitedata_countries)
sitedata_countries <- sitedata_countries[,c(1:3,8)]

# then for each country and year, count the number of tBMS and WCBS sites separately
count_countries <- sitedata_countries %>% group_by(CTRY24NM, YEAR, SURVEY) %>% summarise(no_sites=n())
# for each country, calculate the average number of tBMS and WCBS sites over time
count_countries <- count_countries %>% group_by(CTRY24NM, SURVEY) %>% summarise(avg_no_sites=mean(no_sites))
# change from long to wide
count_countries <- as.data.frame(count_countries)
count_countries2 <- reshape(count_countries, idvar = "CTRY24NM", v.names= c("avg_no_sites"), timevar = "SURVEY", direction = "wide")
count_countries2$total_sites <- count_countries2$avg_no_sites.WCBS + count_countries2$avg_no_sites.tBMS
# change these into a proportion
count_countries2$prop_tBMS <- count_countries2$avg_no_sites.tBMS/count_countries2$total_sites
count_countries2$prop_WCBS <- count_countries2$avg_no_sites.WCBS/count_countries2$total_sites

## repeat for regions in England
sitedata_england <- sitedata_countries[sitedata_countries$CTRY24NM=="England",]

bmsdata <- readRDS("Data/UKBMS/Formatted/sp_weeklycount_1976-2024_zerofilled.rds")
bmsdata <- bmsdata[bmsdata$YEAR>=2009,]
bmsdata <- bmsdata %>% dplyr::mutate(SURVEY = ifelse(SITENO<50000, "tBMS", "WCBS"))

sitedata <- unique(bmsdata[,c("SITENO","YEAR", "SURVEY")])
site_info <- read.csv("Data/UKBMS/site_location_data_2024.csv", header=TRUE)
sitedata <- merge(x = sitedata, y = site_info[ , c("SITENO", "GRIDREF", "EAST", "NORTH", "LONGITUDE", "LATITUDE")], by = "SITENO", all.x=TRUE)

sitedata <- sitedata[!grepl("WV", sitedata$GRIDREF),]
sitedata <- sitedata[!grepl("WA", sitedata$GRIDREF),]

sitedata_england <- merge(sitedata_england, sitedata, by=c("SITENO", "SURVEY", "YEAR"), all.x=TRUE)
sitedata_england$CTRY24NM <- NULL

regions <- st_read("Data/Boundaries/RGN_DEC_2025_EN_BUC.shp")
# Ensure CRS match
sitedata_england <- na.omit(sitedata_england)
sitedata_england <- st_as_sf(sitedata_england, coords=c("EAST", "NORTH"))
sitedata_england <- sitedata_england %>% st_set_crs(st_crs(regions))

sitedata_regions <- st_join(sitedata_england, regions, join = st_intersects)
# some of these have NAs because the grid reference is technically in the sea
# remove these and run through st_nearest_feature instead
sitedata_nas <- sitedata_regions[is.na(sitedata_regions$RGN25NM),]
sitedata_nas <- sitedata_nas[,c(1:6,15)]
sitedata_regions_nas <- st_join(sitedata_nas, regions, join = st_nearest_feature)
sitedata_regions <- na.omit(sitedata_regions)
sitedata_regions <- rbind(sitedata_regions, sitedata_regions_nas)
# select columns of interest
sitedata_regions <- as.data.frame(sitedata_regions)
sitedata_regions <- sitedata_regions[,c(1:3,8)]

# then for each country and year, count the number of tBMS and WCBS sites separately
count_regions <- sitedata_regions %>% group_by(RGN25NM, YEAR, SURVEY) %>% summarise(no_sites=n())
# for each country, calculate the average number of tBMS and WCBS sites over time
count_regions <- count_regions %>% group_by(RGN25NM, SURVEY) %>% summarise(avg_no_sites=mean(no_sites))
# change from long to wide
count_regions <- as.data.frame(count_regions)
count_regions$geometry <- NULL
count_regions2 <- reshape(count_regions, idvar = "RGN25NM", v.names= c("avg_no_sites"), timevar = "SURVEY", direction = "wide")
count_regions2$total_sites <- count_regions2$avg_no_sites.WCBS + count_regions2$avg_no_sites.tBMS
# change these into a proportion
count_regions2$prop_tBMS <- count_regions2$avg_no_sites.tBMS/count_regions2$total_sites
count_regions2$prop_WCBS <- count_regions2$avg_no_sites.WCBS/count_regions2$total_sites

# take england out of countries one and add in regions
count_countries2 <- count_countries2[!count_countries2$CTRY24NM=="England",]
colnames(count_countries2)[1] <- "Area"
colnames(count_regions2)[1] <- "Area"
count_countries_regions <- rbind(count_countries2, count_regions2)

# merge in the boundaries - remove england and add regions
countries <- countries[!countries$CTRY24NM=="England",]
countries <- subset(countries, select = -c(CTRY24CD, CTRY24NMW))
colnames(countries)[1] <- "Area"
regions <- subset(regions, select = -c(RGN25CD, RGN25NMW))
colnames(regions)[1] <- "Area"
countries_regions <- rbind(countries, regions)
countries_regions <- merge(countries_regions, count_countries_regions, by="Area")

## plot map
WCBS_sites_map <- ggplot(data=countries_regions) + geom_sf(aes(fill=prop_WCBS))+
  scale_fill_continuous()+
  guides(fill=guide_legend(title="Proportion of \nWCBS sites"))+
  theme_bw()+
  theme(text = element_text(size = 16))
WCBS_sites_map
ggsave(WCBS_sites_map, file="Output/Figures/Proportion_WCBS_sites.png", height=22, width=28, units="cm")
