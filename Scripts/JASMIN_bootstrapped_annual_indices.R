## Bootstrapping confidence intervals for annual indices on JASMIN

require('rslurm')
library(rbms)
library(tidyverse)
library(data.table)

setwd("/home/users/lhordley/tbms_bootstrapping")

rm(list = ls())
print(getwd())


  ## tBMS indices --------

 # read in data
 sindex_tbms <- readRDS("Site_indices_tBMS.rds")

  boot_fun <- function(species, iteration){

    library(rbms)
    library(tidyverse)
    library(data.table)

   print(paste("Species", species, sep=" "))
   print(paste("Iteration", iteration, sep=" "))

    # tBMS data
    sindex_temp <- sindex_tbms[sindex_tbms$SPECIES==species,]
    sindex_temp <- as.data.table(sindex_temp)
    set.seed(218795)
    bootsample <- rbms::boot_sample(sindex_temp, boot_n = 1000)

    co_index <- list()

    for(i in 1:50){
      print(i)
      co_index[[i+((iteration-1)*50)]] <- rbms::collated_index(data = sindex_temp, s_sp = species, sindex_value = "SINDEX",
                                              bootID=i+((iteration-1)*50), boot_ind= bootsample, glm_weights=TRUE, rm_zero=TRUE)
    }

    # save output
    co_index <- rbindlist(lapply(co_index, FUN = "[[","col_index"))
    co_index <- co_index[COL_INDEX > 0.0001 & COL_INDEX < 100000, ]
    co_index$SPECIES <- species
    
    saveRDS(co_index, paste("Annual_index_tBMS_spp_", species, "_iteration_", iteration, ".rds", sep=""))

  }

  ## Create roster
  species_list <- read.csv("Species_list_GAM.csv", header=TRUE)
  species<-species_list$SPECIES
  iteration<-rep(1:20,23)

  roster <- data.frame(species = species, iteration = iteration)

  ## submit job

  # Create the job scipt and the R script needed to run the process on
  # lotus using slurm. Note: you can edit the templates used. These are
  # found in the slurm folder in your R library (run '.Library' to find).
  # You will need to add the command to load jaspy: module add jaspy
  sjob <- slurm_apply(f = boot_fun,
                      params = roster,
                      jobname = 'tBMS_boot',
                      nodes = nrow(roster),
                      cpus_per_node = 1,
                      submit = TRUE,
                      global_objects = 'sindex_tbms',
                      slurm_options = list(account = 'no-project',
                                           time = '23:59:00',
                                           mem = 60000,
                                           partition = 'standard',
                                           qos = 'standard',
                                           error = '%a.err'))




## WCBS indices --------

setwd("/home/users/lhordley/wcbs_bootstrapping")

rm(list = ls())
print(getwd())

sindex_wcbs <- readRDS("Site_indices_WCBS.rds")

boot_fun <- function(species, iteration){

  library(rbms)
  library(tidyverse)
  library(data.table)

  print(paste("Species", species, sep=" "))
  print(paste("Iteration", iteration, sep=" "))

  # WCBS data
  sindex_temp <- sindex_wcbs[sindex_wcbs$SPECIES==species,]
  sindex_temp <- as.data.table(sindex_temp)
  set.seed(218795)
  bootsample <- rbms::boot_sample(sindex_temp, boot_n = 1000)


  co_index <- list()

  for(i in 1:100){
    print(i)
    co_index[[i+((iteration-1)*100)]] <- rbms::collated_index(data = sindex_temp, s_sp = species, sindex_value = "SINDEX",
                                                              bootID=i+((iteration-1)*100), boot_ind= bootsample, glm_weights=TRUE, rm_zero=TRUE)
  }

  # save output
  co_index <- rbindlist(lapply(co_index, FUN = "[[","col_index"))
  co_index <- co_index[COL_INDEX > 0.0001 & COL_INDEX < 100000, ]
  co_index$SPECIES <- species

  saveRDS(co_index, paste("Annual_index_WCBS_spp_", species, "_iteration_", iteration, ".rds", sep=""))


}

# Create roster
species_list <- read.csv("Species_list_GAM.csv", header=TRUE)
species<-species_list$SPECIES
iteration<-rep(1:10,23)

roster <- data.frame(species = species, iteration = iteration)

# test <- lapply(roster$species, boot_fun)

## submit job

# Create the job scipt and the R script needed to run the process on
# lotus using slurm. Note: you can edit the templates used. These are
# found in the slurm folder in your R library (run '.Library' to find).
# You will need to add the command to load jaspy: module add jaspy
sjob <- slurm_apply(f = boot_fun,
                    params = roster,
                    jobname = 'WCBS_boot',
                    nodes = nrow(roster),
                    cpus_per_node = 1,
                    submit = TRUE,
                    global_objects = 'sindex_wcbs',
                    slurm_options = list(account = 'no-project',
                                         time = '23:59:00',
                                         mem = 60000,
                                         partition = 'standard',
                                         qos = 'standard',
                                         error = '%a.err'))
