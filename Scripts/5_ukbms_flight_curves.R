# Calculate flight curves using all UKBMS data (tBMS and WCBMS)
# Using package 'rbms'

# load packages
library(rbms)
library(dplyr)
library(lubridate)
library(ggplot2)
library(data.table)
library(ggrepel)
library(foreach)
library(doParallel)

bmsdata <- readRDS("Data/UKBMS/Formatted/sp_weeklycount_1976-2024_zerofilled.rds")

# filter to species of interest
splist <- read.csv(file="Data/UKBMS/Species_list_GAM.csv", header=TRUE)
bmsdata <- bmsdata[bmsdata$SPECIES %in% splist$SPECIES,]
length(unique(bmsdata$SPECIES))

#################################
## From counts to flight curve ##
#################################

## 1. Create visit data
bmsdata <- bmsdata %>% 
  mutate(DATE = make_date(day=DAY, month=MONTH, year=YEAR))
ukbms_visit <- unique(bmsdata[,c("SITENO", "DATE")])
colnames(ukbms_visit)[1] <- "SITE_ID"

## 2. Create count data
ukbms_count <- bmsdata[,c("SITENO","DATE","COMMON_NAME","DAY","MONTH","YEAR","COUNT")]
colnames(ukbms_count)[1] <- "SITE_ID"
colnames(ukbms_count)[3] <- "SPECIES"

## 3. Initialise a time-series
ts_date <- rbms::ts_dwmy_table(InitYear = 2009, LastYear = 2024, WeekDay1 = 'monday')

# 3a. Add monitoring data to the time-series

ts_season <- rbms::ts_monit_season(ts_date, StartMonth = 4, EndMonth = 9, StartDay = 1, 
                                   EndDay = NULL, CompltSeason = TRUE, Anchor = TRUE, 
                                   AnchorLength = 2, AnchorLag = 2, TimeUnit = 'w')
# 3b. Add site visits to the time-series
ts_season_visit <- rbms::ts_monit_site(ts_season, ukbms_visit)

## 4. Add observed counts and compute yearly flight curve for the data
nCores <- detectCores() - 1
cl <- makeCluster(nCores) # Define number of CPUs
registerDoParallel(cl)

#options(warn=2)
species <- unique(bmsdata$COMMON_NAME)
pheno_final <- NULL
start_time <- Sys.time() 
for(i in species){
print(i)
  ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, ukbms_count, sp = i)
  ts_flight_curve <- rbms::flight_curve(ts_season_count, NbrSample = 300, MinVisit = 5, MinOccur = 3, 
                                        MinNbrSite = 5, MaxTrial = 4, GamFamily = 'nb', SpeedGam = FALSE, 
                                        CompltSeason = TRUE, SelectYear = NULL, TimeUnit = 'w')
  
  # extract pheno object: contains the shape of annual flight curves, standardised to sum 1
  
  pheno <- ts_flight_curve$pheno
  pheno$sp <- i
  pheno_final <- rbind(pheno, pheno_final)
}
end_time <- Sys.time() # takes ~ 5.5 hours to calculate flight curves

saveRDS(pheno_final, file="Data/UKBMS/Flight curves/All_UKBMS_weekly_anchor_flight_curves.rds")
length(unique(pheno_final$SPECIES)) # 23 species

# plot flight curves
for (i in unique(pheno_final$SPECIES)) {
  print(i)
  flight_curve <- ggplot(na.omit(pheno_final[pheno_final$SPECIES==i,]), aes(x=trimWEEKNO, y=NM, colour=M_YEAR))+
    geom_line()+
    labs(x="Monitoring Week", y="Relative Abundance", colour="")+
    ggtitle(i)+
    #scale_x_continuous(breaks=seq(0,366, by=50)) +
    theme_classic()+
    theme(title = element_text(size = 14), legend.text=element_text(size=12),axis.text=element_text(size=12))
  ggsave(flight_curve, file=paste0("Graphs/Flight curves/Flight_curve_", i,".png"), width = 20, height = 15, units = "cm")
  Sys.sleep(2)
}
dev.off()





