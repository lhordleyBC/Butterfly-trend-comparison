# Compare land cover of WCBS and tBMS squares 
# Uses UKCEH Land Cover Map data 20215

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

## Tidy up data and assign sites a 1k and 10k grid reference for plotting maps

# Read in UKBMS data
bmsdata <- readRDS("Data/UKBMS/Formatted/sp_weeklycount_1976-2024_zerofilled.rds")

# Filter to 2009 to 2024 only
sitedata <- bmsdata[YEAR >= 2009]
# Create survey ID column
sitedata <- sitedata %>% dplyr::mutate(SURVEY = ifelse(SITENO<50000, "tBMS", "WCBS"))
# filter sites to match those used in GAM analysis 
# Apr-Sept for tBMS sites, Jul-Aug for WCBS sites
sitedata <- sitedata %>% filter(case_when(
  SURVEY=="WCBS" ~ MONTH %in% 7:8,
  SURVEY=="tBMS" ~ MONTH %in% 4:9
))
# also filter to species of interest to get sites to match
splist <- read.csv(file="Data/UKBMS/Species_list_GAM.csv", header=TRUE)
sitedata <- sitedata[sitedata$COMMON_NAME %in% splist$COMMON_NAME,]
setDT(sitedata)[, .(count = uniqueN(SITENO)), by = SURVEY]
# 3558 for tBMS
# 2145 for WCBS

# take columns of interest
sitedata <- unique(sitedata[,c("SITENO","YEAR","SURVEY")])

length(unique(sitedata$SITENO)) # 5,703 sites total
length(unique(sitedata$SITENO[sitedata$SURVEY=="tBMS"])) # 3,558 tBMS sites between April and September from 2009-2024
length(unique(sitedata$SITENO[sitedata$SURVEY=="WCBS"])) # 2,145 WCBS squares in July and August from 2009-2024

# Read in site info
site_info <- read.csv("Data/UKBMS/site_location_data_2024.csv", header=TRUE)
sitedata <- merge(x = sitedata, y = site_info[ , c("SITENO", "GRIDREF", "EAST", "NORTH", "LONGITUDE", "LATITUDE")], by = "SITENO", all.x=TRUE)
# add in channel island one too
channel_isle_info <- read.csv("Data/UKBMS/Channel_Isles_location_data.csv")
channel_isle_info <- channel_isle_info %>% dplyr::select(-Site.Name)
names(channel_isle_info) <- c("SITENO", "LATITUDE", "LONGITUDE") 
sitedata_nas <- sitedata %>% filter(if_any(everything(), is.na))
sitedata_nas <- dplyr::select(sitedata_nas, -c('LATITUDE', 'LONGITUDE'))
sitedata_nas2 <- merge(sitedata_nas, channel_isle_info, by="SITENO", all.x=TRUE)
sitedata_nas2$GRIDREF <- c("WV38", "WV27", "WV37", "WV38", "WV27", "WV27", "WV27", "WV55", "WV54", "WV54")
sitedata <- na.omit(sitedata)
sitedata <- dplyr::bind_rows(sitedata, sitedata_nas2)
# this is waaaay too long but it works

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
# what about turnover?
library(codyn)
sitedata <- sitedata[,c(1:3)]
sitedata$ABUNDANCE <- 1

site_turnover <- turnover(df = sitedata, 
                         time.var = "YEAR", 
                         species.var = "SITENO", 
                         abundance.var = "ABUNDANCE", 
                         replicate.var = "SURVEY")
total_turnover <- ggplot(data=site_turnover, aes(x=YEAR, y=total, colour=SURVEY))+
  geom_line(lwd=1)+
  labs(x="Year", y="Total turnover")+
  scale_x_continuous(n.breaks = 10)+
  theme_bw()+
  theme(text = element_text(size = 16),
        legend.title=element_blank())
total_turnover
ggsave(total_turnover, file="Output/Figures/Total_turnover.png", height=8, width=15, units="cm")

site_turnover %>% group_by(SURVEY) %>% summarise(mean_turnover=mean(total))
# WCBS = 0.363
# tBMS = 0.235
# WCBS has higher turnover compared to tBMS
## this is driven more by sites disappearing in WCBS over time

site_appearance <- turnover(df = sitedata, 
                           time.var = "YEAR",
                           species.var = "SITENO",
                           abundance.var = "ABUNDANCE",
                           replicate.var = "SURVEY",
                           metric = "appearance")
site_appearance %>% group_by(SURVEY) %>% summarise(mean_appearance=mean(appearance))
# WCBS = 0.183
# tBMS = 0.143
## WCBS has slightly more sites appear over time
site_disappearance <- turnover(df = sitedata, 
                            time.var = "YEAR",
                            species.var = "SITENO",
                            abundance.var = "ABUNDANCE",
                            replicate.var = "SURVEY",
                            metric = "disappearance")
site_disappearance %>% group_by(SURVEY) %>% summarise(mean_appearance=mean(disappearance))
# WCBS = 0.179
# tBMS = 0.0919
## WCBS has more sites disappear over time compared to tBMS


###
# some gridrefs (tBMS) are 100m x 100m
unique(nchar(sitedata$GRIDREF)) # varying lengths: 0, 5, 6, 7, 8, 10, 12
# remove year - not needed
all_sites <- unique(sitedata[,c("SITENO", "GRIDREF", "EAST", "NORTH", "LATITUDE", "LONGITUDE", "SURVEY")])
unique(nchar(all_sites$GRIDREF)) # varying lengths: 5, 6, 7, 8, 10, 12
unique(nchar(all_sites$EAST)) # 5 and 6
# NAs are the 10 Channel Island sites with only lat/lon data

# try plotting a map - tBMS
names <- c("GBR", "IMN", "GGY", "JEY")
UK <- gisco_get_countries(country = names, resolution = 1)
tbms_sites <- ggplot() + 
  geom_sf(data = UK, fill = "grey", alpha = 0.3) + 
  geom_point(data = all_sites[all_sites$SURVEY=="tBMS",], 
             aes(x = as.numeric(LONGITUDE), 
                 y = as.numeric(LATITUDE)), shape=20, size=1) + 
  theme_void() +
  theme(title = element_text(size = 12))
tbms_sites
# WCBS
wcbs_sites <- ggplot() + 
  geom_sf(data = UK, fill = "grey", alpha = 0.3) + 
  geom_point(data = all_sites[all_sites$SURVEY=="WCBS",], 
             aes(x = as.numeric(LONGITUDE), 
                 y = as.numeric(LATITUDE)), shape=20, size=1) + 
  theme_void() +
  theme(title = element_text(size = 12))
wcbs_sites

# map for each year - not very clear
ggplot() + 
  facet_wrap(~YEAR)+
  geom_sf(data = UK, fill = "grey", alpha = 0.3) + 
  geom_point(data = sitedata, 
             aes(x = as.numeric(LONGITUDE), 
                 y = as.numeric(LATITUDE), colour=SURVEY), shape=20, size=1) + 
  theme_void() +
  theme(title = element_text(size = 12))

# number of sites over time
x <- sitedata %>%
  group_by(YEAR, SURVEY) %>%
  summarise(nsites=n())

ggplot(x, aes(x=YEAR, y=nsites, colour=SURVEY))+
  geom_point()+
  geom_line()+
  theme_minimal()
# WCBS stayed fairly steady from 2009, tBMS sites have increased exponentially since 2009 (from 941 to 2281)

## Add in 10k grid references
all_sites <- all_sites %>% mutate(HECTAD=case_when(
  nchar(GRIDREF)==4 ~ GRIDREF,
  nchar(GRIDREF)==5 ~ paste(substr(GRIDREF, 1,2), substr(GRIDREF, 4,4), sep=""),
  nchar(GRIDREF)==6 ~ paste(substr(GRIDREF, 1,3), substr(GRIDREF, 5,5), sep=""),
  nchar(GRIDREF)==7 ~ paste(substr(GRIDREF, 1,2), substr(GRIDREF, 5,5), sep=""),
  nchar(GRIDREF)==8 ~ paste(substr(GRIDREF, 1,3), substr(GRIDREF, 6,6), sep=""),
  nchar(GRIDREF)==10 ~ paste(substr(GRIDREF, 1,3), substr(GRIDREF, 7,7), sep=""),
  nchar(GRIDREF)==12 ~ paste(substr(GRIDREF, 1,3), substr(GRIDREF, 8,8), sep=""),
))

## Add in 1k grid references (if already 1k, stay the same, if longer/finer scale, make it shorter)
all_sites <- all_sites %>% mutate(MONAD=case_when(
  nchar(GRIDREF)==4 ~ GRIDREF, # these are some channel island gridrefs - don't need to be finer scale
  nchar(GRIDREF)==5 ~ GRIDREF,
  nchar(GRIDREF)==6 ~ GRIDREF,
  nchar(GRIDREF)==7 ~ paste(substr(GRIDREF, 1,3), substr(GRIDREF, 5,6), sep=""),
  nchar(GRIDREF)==8 ~ paste(substr(GRIDREF, 1,4), substr(GRIDREF, 6,7), sep=""),
  nchar(GRIDREF)==10 ~ paste(substr(GRIDREF, 1,4), substr(GRIDREF, 7,8), sep=""),
  nchar(GRIDREF)==12 ~ paste(substr(GRIDREF, 1,4), substr(GRIDREF, 8,9), sep=""),
))

# first take Channel Island grid refs separately
# these don't convert using osg_parse, so do it manually
# only need the 10k lat/lon to plot map below
# as these and IoM sites are removed for landscape + protected area analysis
all_sites_ci <-  all_sites[grep('^W', all_sites$HECTAD),]
all_sites_ci$EAST_1K <- NA
all_sites_ci$NORTH_1K <- NA
all_sites_ci$EAST_10K <- NA
all_sites_ci$NORTH_10K <- NA
ci_lat_lon <- read.csv(file="Data/Channel_islands_hectad_lat_lon.csv", header=TRUE)
all_sites_ci <- merge(all_sites_ci, ci_lat_lon, by="HECTAD")
all_sites_ci$lon_1k <- NA
all_sites_ci$lat_1k <- NA

## Split into GB and NI 
## GB
all_sites_gb <-  all_sites[!grep('^W', all_sites$HECTAD),] # remove channel island sites
all_sites_gb <- all_sites_gb %>% filter(nchar(all_sites_gb$MONAD)==6) # keep only 6 figure GRs (i.e. remove NI sites)

# change grid refs to easting/northing 
gb_east_north <- as.data.frame(osg_parse(grid_refs=all_sites_gb$HECTAD))
# move easting northing to middle of hectad
gb_east_north$easting <- gb_east_north$easting + 5000
gb_east_north$northing <- gb_east_north$northing + 5000
colnames(gb_east_north) <- c("EAST_10K", "NORTH_10K")
all_sites_gb <- cbind(all_sites_gb, gb_east_north)
# same for monad 
gb_east_north2 <- as.data.frame(osg_parse(grid_refs=all_sites_gb$MONAD))
# move easting northing to middle of hectad
gb_east_north2$easting <- gb_east_north2$easting + 500
gb_east_north2$northing <- gb_east_north2$northing + 500
colnames(gb_east_north2) <- c("EAST_1K", "NORTH_1K")
all_sites_gb <- cbind(all_sites_gb, gb_east_north2)

# change easting/northing to lat/ lon
gb_lat_lon <- all_sites_gb %>%
  st_as_sf(coords = c("EAST_10K", "NORTH_10K"), crs = 27700) %>%
  st_transform(4326) %>%
  st_coordinates() %>%
  as_tibble()
colnames(gb_lat_lon) <- c("lon_10k","lat_10k")
all_sites_gb <- cbind(all_sites_gb, gb_lat_lon)
# same for monad
gb_lat_lon2 <- all_sites_gb %>%
  st_as_sf(coords = c("EAST_1K", "NORTH_1K"), crs = 27700) %>%
  st_transform(4326) %>%
  st_coordinates() %>%
  as_tibble()
colnames(gb_lat_lon2) <- c("lon_1k","lat_1k")
all_sites_gb <- cbind(all_sites_gb, gb_lat_lon2)

#### Nothern Ireland
all_sites_ni <- all_sites %>% filter(nchar(all_sites$MONAD)==5)
# change grid refs to easting/northing 
ni_east_north <- as.data.frame(igr_to_ig(all_sites_ni$HECTAD, centroids = TRUE))
colnames(ni_east_north) <- c("EAST_10K", "NORTH_10K")
all_sites_ni <- cbind(all_sites_ni, ni_east_north)
# same for monad
ni_east_north <- as.data.frame(igr_to_ig(all_sites_ni$MONAD, centroids = TRUE))
colnames(ni_east_north) <- c("EAST_1K", "NORTH_1K")
all_sites_ni <- cbind(all_sites_ni, ni_east_north)

## grid reference to lat/lon
ni_hectads <- data.frame(igr=all_sites_ni$HECTAD)
ni_lat_lon <- st_igr_as_sf(ni_hectads, "igr", crs = 4326, centroids = TRUE)
ni_lat_lon <- ni_lat_lon %>%
  dplyr::mutate(lon_1k = sf::st_coordinates(.)[,1],
                lat_1k = sf::st_coordinates(.)[,2])
ni_lat_lon <- st_drop_geometry(ni_lat_lon)
ni_lat_lon <- ni_lat_lon[,-1]
all_sites_ni <- cbind(all_sites_ni, ni_lat_lon)
# same for monad
ni_monads <- data.frame(igr=all_sites_ni$MONAD)
ni_lat_lon2 <- st_igr_as_sf(ni_monads, "igr", crs = 4326, centroids = TRUE)
ni_lat_lon2 <- ni_lat_lon2 %>%
  dplyr::mutate(lon_10k = sf::st_coordinates(.)[,1],
                lat_10k = sf::st_coordinates(.)[,2])
ni_lat_lon2 <- st_drop_geometry(ni_lat_lon2)
ni_lat_lon2 <- ni_lat_lon2[,-1]
all_sites_ni <- cbind(all_sites_ni, ni_lat_lon2)

## Add datasets together
all_sites_ci <- all_sites_ci[,c(2:8,1,9,12:13,10,11,14:17)]
all_sites_final <- rbind(all_sites_gb, all_sites_ni, all_sites_ci) # 5,703 sites
length(unique(all_sites_final$SITENO[all_sites_final$SURVEY=="tBMS"])) # 3,558 tBMS sites in April-September from 2009-2024 
length(unique(all_sites_final$SITENO[all_sites_final$SURVEY=="WCBS"])) # 2,145 WCBS squares in July and August from 2009-2024
# save file
write.csv(all_sites_final, file="Data/UKBMS/WCBS_tBMS_sites_location.csv", row.names=FALSE)

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
                       font.label=list(color="black",size=20))
site_maps

ggsave(site_maps, file="Output/Figures/Site_maps.png", height=22, width=28, units="cm")


### count the number of tBMS and WCBS sites in each country and each region in England 
bmsdata <- readRDS("Data/UKBMS/Formatted/sp_weeklycount_1976-2024_zerofilled.rds")
bmsdata <- bmsdata[bmsdata$YEAR>=2009,]
bmsdata <- bmsdata %>% dplyr::mutate(SURVEY = ifelse(SITENO<50000, "tBMS", "WCBS"))

sitedata <- unique(bmsdata[,c("SITENO","YEAR", "SURVEY")])
site_info <- read.csv("Data/UKBMS/site_location_data_2024.csv", header=TRUE)
sitedata <- merge(x = sitedata, y = site_info[ , c("SITENO", "GRIDREF", "EAST", "NORTH", "LONGITUDE", "LATITUDE")], by = "SITENO", all.x=TRUE)

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
# try and plot as a map? 
# take england out of the country boundaries file and merge with regions
# then plot graph, maybe numbers showing proportion of either tBMS or WCBS? Or both if they can fit

WCBS_sites_map <- ggplot(data=countries_regions) + geom_sf(aes(fill=prop_WCBS))+
  scale_fill_continuous()+
  guides(fill=guide_legend(title="Proportion of \nWCBS sites"))+
  theme_bw()+
  theme(text = element_text(size = 16))
WCBS_sites_map
ggsave(WCBS_sites_map, file="Output/Figures/Proportion_WCBS_sites.png", height=22, width=28, units="cm")




######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

## LAND COVER ANALYSIS
all_sites_final <- read.csv(file="Data/UKBMS/WCBS_tBMS_sites_location.csv", header=TRUE)

#### GREAT BRITAIN
## Add in LCM 2015 data - GB ##
library(raster)
lcm_gb = stack("Data/Land cover/LCM2015_GB_1km_percent_cover_aggregate_class.tif")
crs(lcm_gb) ## projection
res(lcm_gb)
#plot(lcm_gb)

# remove NI sites from all_sites
all_sites_gb <- all_sites_final %>% filter(nchar(all_sites_final$MONAD)==6)
# create unique monad dataset
monad_gb <- unique(all_sites_gb[,c(9,12:13,7)]) # 5111
# extract LCM for each 1km square
bms_lcm_gb <- raster::extract(lcm_gb, monad_gb[,2:3], df=T)
bms_lcm_gb <- cbind(bms_lcm_gb, monad_gb)
bms_lcm_gb <- bms_lcm_gb[,-1]
bms_lcm_gb <- bms_lcm_gb[, c(11:14,1:10)]

#### NORTHERN IRELAND
## Add in LCM 2015 data 
lcm_ni = stack("Data/Land cover/lcm2015_ni_1km_percent_cover_aggregate_class.tif")
crs(lcm_ni) ## projection
res(lcm_ni)
#plot(lcm_ni)

# keep NI sites from all_sites
all_sites_ni <- all_sites_final %>% filter(nchar(all_sites_final$MONAD)==5)
monad_ni <- unique(all_sites_ni[,c(9,12:13,7)]) # 115 

bms_lcm_ni <- raster::extract(lcm_ni, monad_ni[,2:3], df=T)
bms_lcm_ni <- cbind(bms_lcm_ni, monad_ni)
bms_lcm_ni <- bms_lcm_ni[,-1]
bms_lcm_ni <- bms_lcm_ni[, c(11:14,1:10)]

## bind the two datasets together
# change column names first so they match 
colnames(bms_lcm_gb)[5:14] <- c("PERC_BROADLEAF", "PERC_CONIFEROUS", "PERC_ARABLE", "PERC_IMPROVED", "PERC_SEMINATURAL", "PERC_BOG",
                                "PERC_SALTWATER", "PERC_FRESHWATER", "PERC_COASTAL", "PERC_URBAN")
colnames(bms_lcm_ni)[5:14] <- c("PERC_BROADLEAF", "PERC_CONIFEROUS", "PERC_ARABLE", "PERC_IMPROVED", "PERC_SEMINATURAL", "PERC_BOG",
                                "PERC_SALTWATER", "PERC_FRESHWATER", "PERC_COASTAL", "PERC_URBAN")
bms_lcm <- rbind(bms_lcm_gb, bms_lcm_ni) 
bms_lcm <- na.omit(bms_lcm) # remove NAs - channel island and isle of man sites
length(unique(bms_lcm$MONAD)) # 5114 1km squares
# 5180 rows - some duplicates where there is a WCBS and tBMS transect within a 1km square
# setDT(bms_lcm)[, .(count = uniqueN(MONAD)), by = SURVEY]

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
bms_lcm <- bms_lcm[ , !names(bms_lcm) %in% c("PERC_SALTWATER","PERC_FRESHWATER")] 
bms_lcm_l <- melt(setDT(bms_lcm), id.vars = c("SURVEY", "MONAD", "EAST_1K", "NORTH_1K"), variable.name = "LAND_COVER") # change to long format
sum_bms_lcm <- bms_lcm_l %>% group_by(SURVEY, LAND_COVER) %>% summarise(mean=mean(value), se=2*(plotrix::std.error(value)))

lcm_gb_ni <- lcm_gb_ni[ , !names(lcm_gb_ni) %in% c("PERC_SALTWATER","PERC_FRESHWATER")] 
lcm_gb_ni$SURVEY <- "UK"
lcm_gb_ni_l <- melt(setDT(lcm_gb_ni), id.vars = c("x", "y", "SURVEY"), variable.name = "LAND_COVER") # change to long format
sum_gb_ni_lcm <- lcm_gb_ni_l %>% group_by(SURVEY, LAND_COVER) %>% summarise(mean=mean(value), se=2*(plotrix::std.error(value)))

sum_lcm <- rbind(sum_bms_lcm, sum_gb_ni_lcm)

# bind together raw data too for analysis
bms_lcm <- bms_lcm[,-c(1:3)]
lcm_gb_ni <- lcm_gb_ni[,-c(1:2)]
lcm_data <- rbind(bms_lcm, lcm_gb_ni)
write.csv(lcm_data, file="Output/Landscape analysis/LCM_raw_data.csv", row.names=FALSE)
lcm_data <- read.csv("Output/Landscape analysis/LCM_raw_data.csv", header=TRUE)
  
## plot result

sum_lcm$SURVEY <- factor(sum_lcm$SURVEY, levels = c("tBMS", "WCBS", "UK"))
land_cover_labs <- c("PERC_BROADLEAF"="Broadleaf \nwoodland", "PERC_CONIFEROUS"="Coniferous \nwoodland", "PERC_ARABLE"="Arable",
                    "PERC_IMPROVED"="Improved \ngrassland", "PERC_SEMINATURAL"="Semi-natural \ngrassland",
                    "PERC_BOG"="Mountain, heath \nand bog", "PERC_COASTAL"="Coastal", "PERC_URBAN"="Built-up areas \nand gardens")

lcm <- ggplot(data=sum_lcm)+
  facet_wrap(~LAND_COVER, scales="free", labeller=as_labeller(land_cover_labs))+
  geom_point(mapping=aes(x=SURVEY, y=mean, colour=SURVEY))+
  geom_linerange(aes(x=SURVEY, ymin=mean-se, ymax=mean+se, colour=SURVEY))+
  ylab("Percent coverage")+
  xlab("Survey")+
  theme_minimal()+
  theme(text=element_text(size=16),
        legend.title=element_blank())
lcm
# save plot
ggsave(lcm, file="Output/Figures/LCM_comparison.png", height=6, width=8)

# save table of mean + SE
write.csv(sum_lcm, file="Output/Landscape analysis/LCM_mean_SE.csv", row.names=FALSE)
sum_lcm <- read.csv("Output/Landscape analysis/LCM_mean_SE.csv", header=TRUE)


## check for significant differences - kruskall wallis test
library(rstatix)
lcm_data <- read.csv("Output/Landscape analysis/LCM_raw_data.csv", header=TRUE)
hist(lcm_data$PERC_BROADLEAF)
hist(lcm_data$PERC_CONIFEROUS)

lcm_data$SURVEY <- ordered(lcm_data$SURVEY,
                         levels = c("tBMS", "WCBS", "UK"))

# broadleaf woodland
kruskal.test(PERC_BROADLEAF ~ SURVEY, data = lcm_data) # significant overall
broadleaf_results <- lcm_data %>% wilcox_test(PERC_BROADLEAF ~ SURVEY, p.adjust.method = "bonferroni") # all significantly different from each other
ggplot(lcm_data, aes(x=SURVEY, y=PERC_BROADLEAF))+
  geom_boxplot()

# coniferous woodland
hist(lcm_data$PERC_CONIFEROUS)
kruskal.test(PERC_CONIFEROUS ~ SURVEY, data = lcm_data) # significant overall
coniferous_results <- lcm_data %>% wilcox_test(PERC_CONIFEROUS ~ SURVEY, p.adjust.method = "bonferroni") # tBMS and UK not sign. different (***ODD RESULT***)
ggplot(lcm_data, aes(x=SURVEY, y=PERC_CONIFEROUS))+
  geom_boxplot()

# arable
hist(lcm_data$PERC_ARABLE)
kruskal.test(PERC_ARABLE ~ SURVEY, data = lcm_data) # significant overall
arable_results <- lcm_data %>% wilcox_test(PERC_ARABLE ~ SURVEY, p.adjust.method = "bonferroni") # all significantly different from each other
ggplot(lcm_data, aes(x=SURVEY, y=PERC_ARABLE))+
  geom_boxplot()

# improved grassland
hist(lcm_data$PERC_IMPROVED)
kruskal.test(PERC_IMPROVED ~ SURVEY, data = lcm_data) # significant overall
impgrass_results <- lcm_data %>% wilcox_test(PERC_IMPROVED ~ SURVEY, p.adjust.method = "bonferroni") # all significantly different from each other
ggplot(lcm_data, aes(x=SURVEY, y=PERC_IMPROVED))+
  geom_boxplot()

# semi-natural grassland
hist(lcm_data$PERC_SEMINATURAL)
kruskal.test(PERC_SEMINATURAL ~ SURVEY, data = lcm_data) # significant overall
seminatural_results <- lcm_data %>% wilcox_test(PERC_SEMINATURAL ~ SURVEY, p.adjust.method = "bonferroni") # all significantly different from each other (tBMS and WCBS least significant)
ggplot(lcm_data, aes(x=SURVEY, y=PERC_SEMINATURAL))+
  geom_boxplot()

# mountain, heath and bog
kruskal.test(PERC_BOG ~ SURVEY, data = lcm_data) # significant overall
bog_results <- lcm_data %>% wilcox_test(PERC_BOG ~ SURVEY, p.adjust.method = "bonferroni") # all significantly different from each other (slightly odd - tBMS and WCBS similar?)
ggplot(lcm_data, aes(x=SURVEY, y=PERC_BOG))+
  geom_boxplot()

# coastal
kruskal.test(PERC_COASTAL ~ SURVEY, data = lcm_data) # significant overall
coastal_results <- lcm_data %>% wilcox_test(PERC_COASTAL ~ SURVEY, p.adjust.method = "bonferroni") # all significantly different from each other 
ggplot(lcm_data, aes(x=SURVEY, y=PERC_COASTAL))+
  geom_boxplot()

# urban
kruskal.test(PERC_URBAN ~ SURVEY, data = lcm_data) # significant overall
urban_results <- lcm_data %>% wilcox_test(PERC_URBAN ~ SURVEY, p.adjust.method = "bonferroni") # all significantly different from each other 
ggplot(lcm_data, aes(x=SURVEY, y=PERC_URBAN))+
  geom_boxplot()

lcm_long <- melt(setDT(lcm_data), id.vars = c("SURVEY"), variable.name = "LAND_COVER") # change to long format

tapply(lcm_data$PERC_BROADLEAF, lcm_data$SURVEY, median)
tapply(lcm_data$PERC_CONIFEROUS, lcm_data$SURVEY, median)
tapply(lcm_data$PERC_ARABLE, lcm_data$SURVEY, median)
tapply(lcm_data$PERC_IMPROVED, lcm_data$SURVEY, median)
tapply(lcm_data$PERC_SEMINATURAL, lcm_data$SURVEY, median)
tapply(lcm_data$PERC_BOG, lcm_data$SURVEY, median)
tapply(lcm_data$PERC_COASTAL, lcm_data$SURVEY, median)
tapply(lcm_data$PERC_URBAN, lcm_data$SURVEY, median)























#############################
## same land cover analyses as above, but with 21 target classes (rather than 10 aggregate classes)

## LAND COVER ANALYSIS
all_sites_final <- read.csv(file="Data/UKBMS/WCBS_tBMS_sites_location.csv", header=TRUE)

#### GREAT BRITAIN
## Add in LCM 2015 data - GB ##
library(raster)
lcm_gb = stack("Data/Land cover/LCM2015_GB_1km_percent_cover_target_class.tif")
crs(lcm_gb) ## projection
res(lcm_gb)
#plot(lcm_gb)

# remove NI sites from all_sites
all_sites_gb <- all_sites_final %>% filter(nchar(all_sites_final$MONAD)==6)
# create unique monad dataset
monad_gb <- unique(all_sites_gb[,c(9,12:13,7)]) # 5111
# extract LCM for each 1km square
bms_lcm_gb <- raster::extract(lcm_gb, monad_gb[,2:3], df=T)
bms_lcm_gb <- cbind(bms_lcm_gb, monad_gb)
bms_lcm_gb <- bms_lcm_gb[,-1]
bms_lcm_gb <- bms_lcm_gb[, c(22:25,1:21)]

#### NORTHERN IRELAND
## Add in LCM 2015 data 
lcm_ni = stack("Data/Land cover/lcm2015_ni_1km_percent_cover_target_class.tif")
crs(lcm_ni) ## projection
res(lcm_ni)
#plot(lcm_ni)

# keep NI sites from all_sites
all_sites_ni <- all_sites_final %>% filter(nchar(all_sites_final$MONAD)==5)
monad_ni <- unique(all_sites_ni[,c(9,12:13,7)]) # 115 

bms_lcm_ni <- raster::extract(lcm_ni, monad_ni[,2:3], df=T)
bms_lcm_ni <- cbind(bms_lcm_ni, monad_ni)
bms_lcm_ni <- bms_lcm_ni[,-1]
bms_lcm_ni <- bms_lcm_ni[, c(22:25,1:21)]

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
length(unique(bms_lcm$MONAD)) # 5114 1km squares
# 5180 rows - some duplicates where there is a WCBS and tBMS transect within a 1km square
# setDT(bms_lcm)[, .(count = uniqueN(MONAD)), by = SURVEY]

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

lcm_gb_ni <- lcm_gb_ni[ , !names(lcm_gb_ni) %in% c("PERC_BROADLEAF", "PERC_CONIFEROUS", "PERC_ARABLE", "PERC_IMPROVED", "PERC_SALTWATER", 
                                                   "PERC_FRESHWATER", "PERC_SUPRA_LITTORAL_ROCK", "PERC_SUPRA_LITTORAL_SEDIMENT", "PERC_LITTORAL_ROCK",
                                                   "PERC_LITTORAL_SEDIMENT", "PERC_SALTMARSH")] 
lcm_gb_ni$SURVEY <- "UK"
lcm_gb_ni_l <- melt(setDT(lcm_gb_ni), id.vars = c("x", "y", "SURVEY"), variable.name = "LAND_COVER") # change to long format
sum_gb_ni_lcm <- lcm_gb_ni_l %>% group_by(SURVEY, LAND_COVER) %>% summarise(mean=mean(value), se=2*(plotrix::std.error(value)))

sum_lcm <- rbind(sum_bms_lcm, sum_gb_ni_lcm)

## plot result

sum_lcm$SURVEY <- factor(sum_lcm$SURVEY, levels = c("tBMS", "WCBS", "UK"))
land_cover_labs <- c("PERC_NEUTRAL"="Neutral \ngrassland", "PERC_CALCAREOUS"="Calcareous \ngrassland", "PERC_ACID"="Acid \ngrassland",
                     "PERC_FEN_MARSH_SWAMP"="Fen, marsh \nand swamp", "PERC_HEATHER"="Heather", "PERC_HEATHER_GRASSLAND"="Heather \ngrassland",
                     "PERC_BOG"="Bog", "PERC_INLAND_ROCK"="Inland rock", "PERC_URBAN"="Urban", "PERC_SUBURBAN"="Suburban")
sum_lcm <- sum_lcm[sum_lcm$LAND_COVER=c("PERC_NEUTRAL")]
sum_lcm <- sum_lcm[sum_lcm$LAND_COVER %in% c("PERC_NEUTRAL", "PERC_CALCAREOUS", "PERC_ACID", "PERC_FEN_MARSH_SWAMP"), ]

lcm <- ggplot(data=sum_lcm)+
  facet_wrap(~LAND_COVER, scales="free", labeller=as_labeller(land_cover_labs))+
  geom_point(mapping=aes(x=SURVEY, y=mean, colour=SURVEY))+
  geom_linerange(aes(x=SURVEY, ymin=mean-se, ymax=mean+se, colour=SURVEY))+
  ylab("Percent coverage")+
  xlab("Survey")+
  theme_minimal()+
  theme(text=element_text(size=16),
        legend.title=element_blank())
lcm
# save plot
ggsave(lcm, file="Output/Figures/LCM_comparison2.png", height=6, width=6)

# save table of mean + SE
write.csv(sum_lcm, file="Output/Landscape analysis/LCM_mean_SE_grassland_bog_urban.csv", row.names=FALSE)

