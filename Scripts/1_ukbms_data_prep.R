# Prepare latest UKBMS data
# Raw UKBMS data can be found here: F:\Monitoring\Transect\Data\UKBMS Annual Data
# The R script ukbms_data_prep.R does some of the processing/formatting of raw data, adding zeros etc
# by merging count and visit data

# 1. Add zero counts to weekly counts based on visit data

library(data.table)

# Read in UKBMS weekly counts
counts <- fread(paste0("Data/UKBMS/sp_weeklycount_1976-2024.txt"))
counts <- counts[COUNT >= 0] # removes NAs 
counts[,START_TIME:=NULL]

# Visit data
visits <- fread(paste0("Data/UKBMS/VISIT_table_1976-2024_limited.txt"))
setnames(visits, "SITECODE", "SITENO") # change column name to match count data 

# Use collated indices to get species list
splist <- unlist(unique(fread(paste0("Data/UKBMS/Collated_indices_1976-2024.txt"))[,'SPECIES CODE']))

allcounts <- NULL
for(spp in splist){
  # Filter visit table to sites where given species has been counted
  visits_spp <- visits[SITENO %in% counts[SPECIES == spp]$SITENO] # should be the same number of sites as counts_spp - but will 
                                                                  # remove sites where species was counted but no visit data
  # Filter counts to given species
  counts_spp <- counts[SPECIES == spp]
  # Merge visits and counts - adds rows where visit was made but no count
  temp <- dplyr::left_join(visits_spp[,c("SITENO","DAY","MONTH","YEAR", "WEEKNO")],
                           counts_spp,
                           by = c("SITENO", "DAY", "MONTH", "YEAR", "WEEKNO"))
  temp$SPECIES <- spp
  temp$COMMON_NAME <- counts_spp$COMMON_NAME[1]
  # Add zeros where visit made but no count of that species
  temp[is.na(temp$COUNT),]$COUNT <- 0
  
  allcounts <- rbind(allcounts, temp)
}

# Save
saveRDS(allcounts, paste0("Data/UKBMS/sp_weeklycount_1976-2024_zerofilled.rds"))





