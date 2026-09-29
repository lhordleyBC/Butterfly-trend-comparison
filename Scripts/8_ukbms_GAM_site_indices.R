## Here we calculate site indices for tBMS and WCBS separately using the flight curves created previously

# load packages
library(rbms)
library(dplyr)
library(lubridate)
library(data.table)
library(ggrepel)
library(basicPlotteR)
library(parallel)

#######################################
## From flight curve to site indices ##
#######################################

pheno <- readRDS("Data/UKBMS/Flight curves/All_UKBMS_weekly_anchor_flight_curves.rds")
pheno_orangetip <- readRDS("Data/UKBMS/Flight curves/All_UKBMS_weekly_anchor_flight_curves_OT.rds")
names(pheno_orangetip)[names(pheno_orangetip) == 'species'] <- 'sp'
pheno_orangetip$SPECIES <- as.character(pheno_orangetip$SPECIES)
pheno_orangetip$sp <- as.character(pheno_orangetip$sp)
pheno <- bind_rows(pheno, pheno_orangetip)

bmsdata <- readRDS("Data/UKBMS/Formatted/sp_weeklycount_1976-2024_zerofilled.rds")
splist <- read.csv(file="Data/UKBMS/Species_list_GAM.csv", header=TRUE)

# focus on the 23 species filtered by
# min 26 squares per year between 2009 and 2024 between april and september
bmsdata <- bmsdata[bmsdata$SPECIES %in% splist$SPECIES,]
length(unique(bmsdata$SPECIES))
# create date column
bmsdata <- bmsdata %>% 
  mutate(DATE = make_date(day=DAY, month=MONTH, year=YEAR))
# initialise time-series with day-week-month-year info
ts_date <- rbms::ts_dwmy_table(InitYear = 2009, LastYear = 2024, WeekDay1 = 'monday')

## Create season_visit data frames
# tBMS
# add monitoring season - this is Apr-Sept for tBMS 
tbms_season <- rbms::ts_monit_season(ts_date, StartMonth = 4, EndMonth = 9, StartDay = 1, 
                                   EndDay = NULL, CompltSeason = TRUE, Anchor = TRUE, 
                                   AnchorLength = 2, AnchorLag = 2, TimeUnit = 'w')
tbmsdata <- bmsdata[SITENO < 50000]
tbms_count <- tbmsdata[,c("SITENO","DATE","SPECIES","DAY","MONTH","YEAR","COUNT")]
colnames(tbms_count)[1] <- "SITE_ID"
tbms_visit <- unique(tbmsdata[,c("SITENO", "DATE")])
colnames(tbms_visit)[1] <- "SITE_ID"
# Create tbms_season_visit based on tBMS April - September
tbms_season_visit <- rbms::ts_monit_site(tbms_season, tbms_visit)
length(unique(tbms_season_visit$SITE_ID)) # 3558

# WCBS
# add monitoring season - same as tBMS April to September (as often counts are outside July and August)
# and because we can account for phenology here, these counts can be included
wcbs_season <- rbms::ts_monit_season(ts_date, StartMonth = 4, EndMonth = 9, StartDay = 1, 
                                     EndDay = NULL, CompltSeason = TRUE, Anchor = TRUE, 
                                     AnchorLength = 2, AnchorLag = 2, TimeUnit = 'w')
wcbsdata <- bmsdata[SITENO >= 50000]
wcbs_count <- wcbsdata[,c("SITENO","DATE","SPECIES","DAY","MONTH","YEAR","COUNT")]
colnames(wcbs_count)[1] <- "SITE_ID"
wcbs_visit <- unique(wcbsdata[,c("SITENO", "DATE")])
colnames(wcbs_visit)[1] <- "SITE_ID"
# Create tbms_season_visit based on WCBMS July and August only
wcbs_season_visit <- rbms::ts_monit_site(wcbs_season, wcbs_visit)
length(unique(wcbs_season_visit$SITE_ID)) # 2202

# filter by count (>0), filter by month (exclude July and August)
# group by species and year
# count the number of rows
wcbsdata <- wcbsdata[wcbsdata$YEAR>=2009,]
test <- wcbsdata %>% group_by(YEAR, COMMON_NAME) %>% dplyr::filter(MONTH < 7 | MONTH > 8) %>% dplyr::filter(COUNT>0) %>% summarise(n_counts=n())
test2 <- test %>% group_by(COMMON_NAME) %>% summarise(mean_count_outside_season=mean(n_counts))
#### Calculate site indices for tBMS April - September #### 

pheno$SPECIES <- as.integer(pheno$SPECIES)
species <- unique(bmsdata$SPECIES) # 23 species
sindex_final <- NULL
for(i in species) {print(i)
  
  # extract site index for each site, year and species
  # impute_count() function uses the count data generated from ts_season_count() function and the flight curves
  # it looks for the phenology available to estimate and input missing values
  # imputation are made on a weekly basis
  ts_season_count <- rbms::ts_monit_count_site(tbms_season_visit, tbms_count, sp = i)
  tryCatch({
    impt_counts <- rbms::impute_count(ts_season_count=ts_season_count, ts_flight_curve=pheno, YearLimit= NULL, TimeUnit='w')
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  # change Inf to NaN - this issues comes up for Meadow brown where 16 individuals were recorded in May 2010, but the flight curve
  # for that date is zero (anchors might be causing this..)
  impt_counts[sapply(impt_counts, is.infinite)] <- NA
  # this produces a data.table with original counts and imputed counts over the monitoring season, total count per site and year,
  # and total proportion of flight curve covered by the visits and SINDEX = sum of both observed and imputed counts over sampling season
  
  # From the imputed count, the site index can be calculated for each site
  sindex_temp <- rbms::site_index(butterfly_count = impt_counts, MinFC=0.01) 
  sindex_temp$n_observed <- length(na.omit(impt_counts$COUNT[impt_counts$COUNT>0]))
  sindex_temp$SPECIES <- i
  sindex_final <- rbind(sindex_final, sindex_temp)
}

saveRDS(sindex_final,file="Output/GAM abundance analysis/Site_indices_tBMS.rds")

#### Calculate site indices for WCBS April to September #### 

pheno$SPECIES <- as.integer(pheno$SPECIES)
species <- unique(bmsdata$SPECIES) # 23 species
sindex_final <- NULL
for(i in species) {print(i)
  
  # extract site index for each site, year and species
  # impute_count() function uses the count data generated from ts_season_count() function and the flight curves
  # it looks for the phenology available to estimate and input missing values
  # imputation are made on a daily basis
  ts_season_count <- rbms::ts_monit_count_site(wcbs_season_visit, wcbs_count, sp = i)
  tryCatch({
    impt_counts <- rbms::impute_count(ts_season_count=ts_season_count, ts_flight_curve=pheno, YearLimit= NULL, TimeUnit='w')
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  # this produces a data.table with original counts and imputed counts over the monitoring season, total count per site and year,
  # and total proportion of flight curve covered by the visits and SINDEX = sum of both observed and imputed counts over sampling season
  
  # From the imputed count, the site index can be calculated for each site
  sindex_temp<- rbms::site_index(butterfly_count = impt_counts, MinFC=0.01) 
  sindex_temp$n_observed <- length(na.omit(impt_counts$COUNT[impt_counts$COUNT>0]))
  sindex_temp <- sindex_temp[is.finite(sindex_temp$SINDEX),]# remove inf values - where the proportion of flight curve is 0 even when there is a positive count made
  sindex_temp$SPECIES <- i
  sindex_final <- rbind(sindex_final, sindex_temp)
}

saveRDS(sindex_final,file="Output/GAM abundance analysis/Site_indices_WCBS.rds")





