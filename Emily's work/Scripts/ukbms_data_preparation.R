# Prepare latest UKBMS data
# Raw UKBMS data can be found here: F:\Monitoring\Transect\Data\UKBMS Annual Data
# The R script ukbms_data_preparation.R does some of the processing/formatting of raw data, adding zeros etc
# by merging count and visit data

endyear <- 2023

# 1. Add zero counts to weekly counts based on visit data

# 2. Separate out the index only data from the site indices table (based on visit data)

# 1. ----

library(data.table)

# Read in UKBMS weekly counts
counts <- fread(paste0("Data to ", endyear, "/",
                                        "sp_weeklycount_1976-", endyear, ".txt"))
counts <- counts[COUNT >= 0]

# Visit data
visits <- fread(paste0("./Data to ", endyear, "/",
                                   "VISIT_table_1976-", endyear, ".txt"))
setnames(visits, "SITECODE", "SITENO")


# Use collated indices to get species list
splist <- unlist(unique(fread(paste0("./Data to ", endyear, "/Collated_indices_1976-", endyear, ".txt"))[,'SPECIES CODE']))

allcounts <- NULL
for(spp in splist){
  # Filter visit table to sites where given species has been counted
  visits_spp <- visits[SITENO %in% counts[SPECIES == spp]$SITENO]
  # Filter counts to given species
  counts_spp <- counts[SPECIES == spp]
  # Merge visits and counts
  temp <- dplyr::left_join(visits_spp[,c("SITENO","DAY","MONTH","YEAR","WEEKNO","START_TIME")],
                            counts_spp,
                           by = c("SITENO", "DAY", "MONTH", "YEAR", "WEEKNO", "START_TIME"))
  temp$SPECIES <- spp
  temp$COMMON_NAME <- counts_spp$COMMON_NAME[1]
  # Add zeros where visit made but no count of that species
  temp[is.na(temp$COUNT),]$COUNT <- 0

  allcounts <- rbind(allcounts, temp)
}

# Save (csv and rds)
#write.csv(allcounts, "./2018/Formatted/sp_weeklycount_1976-2018_zerofilled.csv")
dir.create(paste0("./Data to ", endyear, "/Formatted/"))
saveRDS(allcounts, paste0("./Data to ", endyear, "/Formatted/sp_weeklycount_1976-", endyear, "_zerofilled.rds"))

# 2. ----

# Index only data
sindices <- fread(paste0("./Data to ", endyear, "/",
                                     "site_indices_1976-", endyear, ".txt"))
colnames(sindices)[colnames(sindices)=="SITE CODE"] <- "SITENO"

# Get the index only data based on cases where there is no visit data
indexonly <- dplyr::anti_join(sindices, visits[,c("SITENO","DAY","MONTH","YEAR","WEEKNO","START_TIME")])


saveRDS(indexonly, paste0("./Data to ", endyear, "/Formatted/sp_indexonly_1976-", endyear, ".rds"))

