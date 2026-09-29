# Compare protected area coverage of WCBMS and tBMS squares 
# Uses WDPA data

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
library(viridis)

# read in data
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
# remove IoM sites
sqIM <- sitedata[sitedata$NORTH_1K < 505000 & sitedata$NORTH_1K > 450000 &
                   sitedata$EAST_1K > 200000 & sitedata$EAST_1K < 250000,]
ggplot(sitedata, aes(EAST_1K, NORTH_1K))+
  geom_point()+coord_fixed()+
  geom_point(aes(EAST_1K, NORTH_1K), data = sqIM, color="red") # looks good
ggplot(sitedata[!sitedata$MONAD %in% sqIM$MONAD,], aes(EAST_1K, NORTH_1K))+
  geom_point()+coord_fixed()
## This removes the isle of man squares
sitedata <- sitedata[!sitedata$MONAD %in% sqIM$MONAD,]
length(unique(sitedata$MONAD)) # 5,164

# Read in protected area data
file_list <- list.files("Data/Protected areas", pattern = "*shp", full.names = TRUE)
shapefile_list <- lapply(file_list, read_sf)
protected <- do.call("rbind", shapefile_list)

summary(protected)
st_crs(protected) # WGS 84
#

# ggplot() + geom_sf(data=bms_sites_buffer[bms_sites_buffer$MONAD=="H8542",])

unique(protected$DESIG)
# Designations of interest: SSSIs, SACs, SPAs, Ramsar sites, NNRs and Local Nature Reserves 
# Based on Cooke et al. 2023
# Area of Special Scientific Interest (NI version of SSI)
# Site of Special Scientific Interest (SSSI)
# Special Protected Area (SPA)
# National Nature Reserve (NNR)
# Ramsar Site, Wetland of International Importance
# Local Nature Reserve
# Site of Special Scientific Interest (Uk)
# Nature Reserve (I think these are all RSPB, or similar)

protected <- protected[protected$DESIG == "Area of Special Scientific Interest" | protected$DESIG == "Site of Special Scientific Interest" |
                         protected$DESIG == "Special Protected Area" | protected$DESIG == "National Nature Reserve" |
                         protected$DESIG == "Ramsar Site, Wetland of International Importance" | protected$DESIG == "Local Nature Reserve" |
                         protected$DESIG == "Site of Special Scientific Interest (Uk)" | protected$DESIG == "Nature Reserve",]
protected<- protected[!protected$NAME %in% c("Larnaca Salt Lake", "Limassol Lake (Akrotiri)"),] 
# 31 variables, 10407 rows
st_write(protected, "Data/Protected areas/protected_areas_filter.gpkg")

protected <- st_union(protected) # large sfc_multipolygon (278.9mb)
st_write(protected, "Data/Protected areas/protected_areas_union.gpkg")
protected_transf <- st_transform(protected, crs = 27700)
st_write(protected_transf, "Data/Protected areas/protected_areas_union_transf.gpkg")


## GB
protected <- st_read("Data/Protected areas/protected_areas_union_transf.gpkg")
st_crs(protected)

bms_sites <- unique(sitedata[,c("MONAD", "NORTH_1K", "EAST_1K")]) # 5164
bms_sites_gb <- bms_sites %>% filter(nchar(bms_sites$MONAD)==6) # 5047 in GB
bms_sites_gb_sf <- st_as_sf(unique(bms_sites_gb, by=c("EAST_1K", "NORTH_1K")), coords = c("EAST_1K","NORTH_1K"), crs = 27700)

bms_sites_gb_buffer <- st_buffer(bms_sites_gb_sf, dist=500, endCapStyle = "SQUARE") # 500m is the radius = 1km square

st_crs(bms_sites_gb_buffer)
pa_bms_gb <- st_intersection(protected, bms_sites_gb_buffer) # 2518 before, should be 5047 

# Calculate area per 1km square/geometry
pa_bms_gb$area <- st_area(pa_bms_gb)
# Convert to df
pa_bms_gb_df <- setDT(st_drop_geometry(pa_bms_gb))
# Add easting northing info
pa_bms_gb_df <- merge(pa_bms_gb_df, bms_sites_gb, by = "MONAD", all.y=TRUE)

# Get the total area per 1km (since some will contain multiple geometries)
pa_bms_gb_df <- pa_bms_gb_df[, .(area = sum(area)), by = .(MONAD, EAST_1K, NORTH_1K)]

# Calculate proportion of each 1km that's protected
pa_bms_gb_df[, prop := as.numeric(area)/(1000*1000)]

# Replace NAs with zeros
pa_bms_gb_df[is.na(pa_bms_gb_df)] <- 0

## Repeat for NI sites
# EPSG:29903
protected <- st_read("Data/Protected areas/protected_areas_union.gpkg")
protected_ni <- st_transform(protected, crs = 29903)
st_crs(protected_ni) # ESPG 29903 projected CRS 

## Filter to NI BMS sites
bms_sites <- unique(sitedata[,c("MONAD", "NORTH_1K", "EAST_1K")])
bms_sites_ni <- bms_sites %>% filter(nchar(bms_sites$MONAD)==5) # 117
bms_sites_ni_sf <- st_as_sf(bms_sites_ni, coords = c("EAST_1K", "NORTH_1K"), crs = 29903, dim = "XY")
bms_sites_ni_buffer <- st_buffer(bms_sites_ni_sf, dist=500, endCapStyle = "SQUARE") # 500m is the radius = 1km square
st_crs(bms_sites_ni_buffer) # ESPG 29903 projected CRS

pa_bms_ni <- st_intersection(protected_ni, bms_sites_ni_buffer) 

# Calculate area per 1km square/geometry
pa_bms_ni$area <- st_area(pa_bms_ni)
# Convert to df
pa_bms_ni_df <- setDT(st_drop_geometry(pa_bms_ni))
# Add easting northing info
pa_bms_ni_df <- merge(pa_bms_ni_df, bms_sites_ni, by = "MONAD", all.y=TRUE)

# Get the total area per 1km (since some will contain multiple geometries)
pa_bms_ni_df <- pa_bms_ni_df[, .(area = sum(area)), by = .(MONAD, EAST_1K, NORTH_1K)]

# Calculate proportion of each 1km that's protected
pa_bms_ni_df[, prop := as.numeric(area)/(1000*1000)]

# Replace NAs with zeros
pa_bms_ni_df[is.na(pa_bms_ni_df)] <- 0

# put both datasets together
pa_coverage <- rbind(pa_bms_gb_df, pa_bms_ni_df)
pa_coverage$perc_cover <- pa_coverage$prop*100
write.csv(pa_coverage, file="Output/Landscape analysis/PA_coverage_bms_monads.csv", row.names=FALSE)

##

pa_coverage <- read.csv("Output/Landscape analysis/PA_coverage_bms_monads.csv", header=TRUE)
length(unique(pa_coverage$MONAD)) # 5164 unique squares 
# 68 1 km squares which have both a tBMS and WCBS site within them 
setDT(pa_coverage)[, .(count = uniqueN(MONAD)), by = SURVEY]
# 3034 tBMS
# 2198 WCBS
# 5164 unique monads (both tBMS and WCBS together)

# plot area
pa_area_sum <- pa_coverage %>% group_by(SURVEY) %>% summarise(mean=mean(area), se=2*(plotrix::std.error(area)))
pa_area_sum$SURVEY <- factor(pa_area_sum$SURVEY, levels = c("tBMS", "WCBS", "Combined"))

pa_area_p <- ggplot(data=pa_area_sum)+
  geom_point(mapping=aes(x=SURVEY, y=mean, colour=SURVEY), size=4)+
  geom_linerange(aes(x=SURVEY, ymin=mean-se, ymax=mean+se, colour=SURVEY), lwd=1)+
  ylab("Area covered by Protected Areas (m^2)")+
  # scale_colour_manual(values=c("#F8766D", "#7CAE00"))+
  theme(strip.text = element_text(size = 10),
        legend.title=element_blank())+
  theme_classic()
pa_area_p

# plot % cover
pa_perc_sum <- pa_coverage %>% group_by(SURVEY) %>% summarise(mean=mean(perc_cover), se=2*(plotrix::std.error(perc_cover)))
pa_perc_sum$SURVEY <- factor(pa_perc_sum$SURVEY, levels = c("tBMS", "WCBS", "Combined"))

pa_perc_p <- ggplot(data=pa_perc_sum)+
  geom_point(mapping=aes(x=SURVEY, y=mean, colour=SURVEY), size=4)+
  geom_linerange(aes(x=SURVEY, ymin=mean-se, ymax=mean+se, colour=SURVEY), lwd=1)+
  ylab("Protected Area \npercent coverage")+
  scale_colour_manual(values=c("#F8766D", "#7CAE00", "#00BFC4"))+
  xlab("Survey")+
  theme_minimal()+
  theme(text=element_text(size=16),
        legend.title=element_blank(),
        legend.position = "none")
pa_perc_p
ggsave(pa_perc_p, file="Output/Figures/Figure3.png", height=5, width=6)
# 21.6% average for tBMS squares and 7.65% for WCBMS squares
# 15.6% for WCBS and tBMS sites together

