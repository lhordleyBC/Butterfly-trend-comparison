#### Calculating annual indices from site indices for tBMS and WCBS separately
## Then reading in bootstrapped annual indices from JASMIN and calculating confidence intervals

library(rbms)
library(tidyverse)
library(data.table)

#### tBMS 

sindex_tbms <- readRDS("Output/GAM abundance analysis/Site_indices_tBMS.rds")

sp_names <- read.csv("Data/UKBMS/Species_list_GAM.csv", header=TRUE)
sindex_tbms <- merge(sindex_tbms, sp_names, by="SPECIES")

species <- unique(sindex_tbms$SPECIES) # 23 species
co_index <- list()

for(i in species) { 
  print(paste("species", i))

  sindex_temp <- sindex_tbms[sindex_tbms$SPECIES==i,]
  co_index <- list()
  
  co_index <- collated_index(data = sindex_temp, s_sp = i, sindex_value = "SINDEX", glm_weights = TRUE, rm_zero = TRUE)
  
  co_index <- co_index$col_index
  co_index <- co_index[COL_INDEX > 0.0001 & COL_INDEX < 100000, ]
  co_index$SPECIES <- i
  
  # save result here
  saveRDS(co_index, file=paste0("Output/GAM abundance analysis/Annual_indices_tBMS/Annual_index_tBMS_spp_", i, ".rds"))
}


# Read in all annual index files - true values and each of the bootstrap files from JASMIN
tbms_files <- list.files(path = "Output/GAM abundance analysis/Annual_indices_tBMS/", pattern = "\\.rds$", full.names = TRUE)
tbms_annual_indices <- do.call("rbind", lapply(tbms_files, readRDS))

length(unique(tbms_annual_indices$SPECIES)) # 23
length(unique(tbms_annual_indices$BOOTi)) # 1001

## Calculate LCI for the true and bootstrap samples and save one file for all results
species <- unique(tbms_annual_indices$SPECIES)
annual_indices_final <- NULL
for(i in species){
  print(i)
  co_index_b <- tbms_annual_indices[tbms_annual_indices$SPECIES==i,]
  co_index_logInd <- co_index_b[BOOTi == 0, .(M_YEAR, COL_INDEX)][, log(COL_INDEX)/log(10), by = M_YEAR][, mean_logInd := mean(V1)]
  
  ## merge the mean log index with the full bootstrap dataset
  data.table::setnames(co_index_logInd, "V1", "logInd"); setkey(co_index_logInd, M_YEAR); setkey(co_index_b, M_YEAR)
  co_index_b <- merge(co_index_b, co_index_logInd, all.x = TRUE)
  
  data.table::setkey(co_index_b, BOOTi, M_YEAR)
  co_index_b[ , boot_logInd := log(COL_INDEX)/log(10)]
  b1 <- data.table(M_YEAR = co_index_b$M_YEAR, LCI = 2 + co_index_b$boot_logInd - co_index_b$mean_logInd)
  b1$BOOTi <- rep(0:1000, each = 16)
  b1 <- merge(b1, co_index_b, by=c("BOOTi", "M_YEAR"))
  b1$sp <- i
  annual_indices_final <- rbind(annual_indices_final, b1)
}

saveRDS(annual_indices_final, file="Output/GAM abundance analysis/Annual_indices_tBMS.rds")


#### WCBS 

sindex_wcbs <- readRDS("Output/GAM abundance analysis/Site_indices_WCBS.rds")
sp_names <- read.csv("Data/UKBMS/Species_list_GAM.csv", header=TRUE)
sindex_wcbs <- merge(sindex_wcbs, sp_names, by="SPECIES")

species <- unique(sindex_wcbs$SPECIES) # 23 species
co_index <- list()

for(i in species) { 
  print(paste("species", i))
  
  sindex_temp <- sindex_wcbs[sindex_wcbs$SPECIES==i,]
  co_index <- list()
  
  co_index <- collated_index(data = sindex_temp, s_sp = i, sindex_value = "SINDEX", glm_weights = TRUE, rm_zero = TRUE)
  
  co_index <- co_index$col_index
  co_index <- co_index[COL_INDEX > 0.0001 & COL_INDEX < 100000, ]
  co_index$SPECIES <- i
  
  # save result here
  saveRDS(co_index, file=paste0("Output/GAM abundance analysis/Annual_indices_WCBS/Annual_index_WCBS_spp_", i, ".rds"))
}

# Read in all annual index files - true values and each of the bootstrap files from JASMIN
wcbs_files <- list.files(path = "Output/GAM abundance analysis/Annual_indices_WCBS/", pattern = "\\.rds$", full.names = TRUE)
wcbs_annual_indices <- do.call("rbind", lapply(wcbs_files, readRDS))

length(unique(wcbs_annual_indices$SPECIES)) # 23
length(unique(wcbs_annual_indices$BOOTi)) # 1001

## Calculate LCI for the true and bootstrap samples and save one file for all results
species <- unique(wcbs_annual_indices$SPECIES)
annual_indices_final <- NULL
for(i in species){
  print(i)
  co_index_b <- wcbs_annual_indices[wcbs_annual_indices$SPECIES==i,]
co_index_logInd <- co_index_b[BOOTi == 0, .(M_YEAR, COL_INDEX)][, log(COL_INDEX)/log(10), by = M_YEAR][, mean_logInd := mean(V1)]

## merge the mean log index with the full bootstrap dataset
data.table::setnames(co_index_logInd, "V1", "logInd"); setkey(co_index_logInd, M_YEAR); setkey(co_index_b, M_YEAR)
co_index_b <- merge(co_index_b, co_index_logInd, all.x = TRUE)

data.table::setkey(co_index_b, BOOTi, M_YEAR)
co_index_b[ , boot_logInd := log(COL_INDEX)/log(10)]
b1 <- data.table(M_YEAR = co_index_b$M_YEAR, LCI = 2 + co_index_b$boot_logInd - co_index_b$mean_logInd)
b1$BOOTi <- rep(0:1000, each = 16)
b1 <- merge(b1, co_index_b, by=c("BOOTi", "M_YEAR"))
b1$sp <- i
annual_indices_final <- rbind(annual_indices_final, b1)
}

saveRDS(annual_indices_final, file="Output/GAM abundance analysis/Annual_indices_WCBS.rds")


