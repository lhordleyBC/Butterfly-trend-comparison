# Growth rate analysis for UKBMS data, limited to wcbs months (July and August) each year

library(speedglm)
library(Matrix)
library(dplyr)
# library(plyr)
library(data.table)
library(plotrix)
library(tidyverse)

# Likelihood function - growth rate model as in Roy et al 2015 paper, but using a concentrated likelihood
ll_func2 <- function(parm, Counts, Visits, l.type="conc"){
  
  # Number of years
  nyears <- ncol(Counts)
  
  # Growth rate parameters
  gr.est <- c(0, parm[1:(nyears-1)])
  par.index <- nyears-1
  
  # Site effects for full likelihood
  if(l.type=="full"){
    Si <- parm[(par.index+1):(par.index+nsites)]
  }
  
  # Take sums of the growth rates over years
  gr.estS <- matrix(cumsum(gr.est), nrow = nsites, ncol = nyears, byrow = TRUE)
  
  
  # Form site effects for concentrated likelihood
  if(l.type=="conc"){
    Si <- log(apply(Counts, 1, sum, na.rm = TRUE)/apply(Visits*exp(gr.estS), 1, sum, na.rm = TRUE)) # growth rate for each site?
  }
  
  # Poisson expectaction
  mu.est <- exp(gr.estS + Si + log(Visits))
  
  # Likelihood formulation
  llik <- dpois(Counts, lambda = mu.est, log = TRUE)
  
  -1*sum(llik, na.rm = TRUE)
}

num <- function(x){as.numeric(as.character(x))}

# Read in UKBMS data
bmsdata <- readRDS("Data/UKBMS/sp_weeklycount_1976-2024_zerofilled.rds")

## Create two files:
# 1. BMS transect data
tbmsdata <- bmsdata[SITENO < 50000]
# 2. WCBS data
wcbsdata <- bmsdata[SITENO >= 50000]
# Filter to wcbs years, July/August only 
tbmsdata <- tbmsdata[YEAR >= 2009 & MONTH %in% 7:8]
wcbsdata <- wcbsdata[YEAR >= 2009 & MONTH %in% 7:8]

# Restrict analyses to species which occur at more than 25 sites in each survey each year
sp_site_count <- wcbsdata[COUNT>0, .(no_sites=uniqueN(SITENO)), by=c("COMMON_NAME","SPECIES", "YEAR")]
sp_site_count2 <- sp_site_count %>% group_by(COMMON_NAME) %>% filter(no_sites > 25) %>% filter(all(c(2009:2024) %in% YEAR))
length(unique(sp_site_count2$COMMON_NAME)) # 22 species (26 in Roy 2015 paper and 24 in Emily's 2021 analysis)
# Extra species compared to Emily's analysis: Painted Lady (2012 recorded at 24 sites) and Orange-tip
# Extra species in Roy's analysis: Painted Lady, Dark-green fritillary, Scotch Argus and Brown argus (orange-tip isn't in Roy paper)
# Roy's filter can't have been a minimum of 25 squares each year as the above species don't fit that for 2009-2013
splist <- unique(sp_site_count2[c("COMMON_NAME", "SPECIES")])
write.csv(splist, file="Data/UKBMS/Species_list_growthrate.csv", row.names=FALSE)
# Filter to 22 species 
tbmsdata <- tbmsdata[SPECIES %in% splist$SPECIES]
length(unique(tbmsdata$COMMON_NAME)) # 22
length(unique(tbmsdata$SITENO)) # 3504
wcbsdata <- wcbsdata[SPECIES %in% splist$SPECIES]
length(unique(tbmsdata$COMMON_NAME)) # 22
length(unique(wcbsdata$SITENO)) # 2145

## 1. Growth rate model for tBMS only

out_tbms_net <- out_tbms_annual <- vcov_tbms <- glm_tbms <- NULL
for(spp in splist){
  cat(spp,"\n")
  # tBMS data for species spp
  tbmsdata_sp <- tbmsdata[SPECIES == spp]
  if(nrow(tbmsdata_sp) == 0) next() # probably not needed here 
  
  # Sum counts over multiple visits
  tbmsdata_sp_sum <- tbmsdata_sp[, .(eCount = sum(COUNT), nV = length(COUNT)), by = .(SITENO, YEAR)]
  # eCount = total counts over visits at each site and year
  # nV = total number of visits at each site and year
  
  colnames(tbmsdata_sp_sum)[1:2] <- c("Site","Year")
  
  # Set up for growth-rate GLM model as in Roy et al. (2015) WCBS paper
  tbmsdata_sp_sum <- tbmsdata_sp_sum[order(tbmsdata_sp_sum$Site, tbmsdata_sp_sum$Year, decreasing=FALSE),]
  
  nsites <- length(unique(tbmsdata_sp_sum$Site))
  nyears <- length(unique(tbmsdata_sp_sum$Year))
  years <- sort(unique(tbmsdata_sp_sum$Year))
  Counts <- as.matrix(reshape(tbmsdata_sp_sum[,c("Site","Year","eCount")],
                              direction = "wide", timevar = "Year", idvar = "Site")[,paste("eCount", years, sep = "."), with = FALSE], ncol = nyears)
  Visits <- as.matrix(reshape(tbmsdata_sp_sum[,c("Site","Year","nV")],
                              direction = "wide", timevar = "Year", idvar = "Site")[,paste("nV", years, sep = "."), with = FALSE], ncol = nyears)
  
  parm <- rep(0, nyears-1)
  cat("GLM for tBMS for",spp,"starting at",date(),"\n")
  stc <- Sys.time()
  # Fit model with concentrated likelihood
  f1_tbms <- try(optim(par = parm,
                      fn = ll_func2,
                      l.type = "conc",
                      Counts = Counts,
                      Visits = Visits,
                      method = "BFGS",
                      hessian = TRUE,
                      control = list(maxit = 2000)),
                silent = TRUE)
  etc <- Sys.time()
  
  cat("GLM for tBMS for", spp, "done at", date(), "\n")
  # Par = annual growth rates 
  # Hessian = square matrix of second-order partial derivatives of a function - 
  # provides info on the curvature of the likelihood function. Inverse of the 
  # hessian often used to esimate variance of maximum likelihood estimators 
  # for remaining parameters 
  
  # Net change estimates up to each year, in matrix form for later calculations
  gr.estS <- matrix(cumsum(c(0, f1_tbms$par)), nrow = nsites, ncol = nyears, byrow = TRUE) 
  # Estimate site indices
  Si.out <- log(apply(Counts, 1, sum, na.rm = TRUE)/apply(Visits * exp(gr.estS), 1, sum, na.rm = TRUE))
  # Estimated total counts
  Fitted <- exp(gr.estS + Si.out + log(Visits))
  # Pearson residuals
  pearsonR <- (Counts - Fitted)/sqrt(Fitted)
  
  # Calculate net change estimate
  ngr <- sum(f1_tbms$par)
  # Pearson chi-square statistic/df
  chat <- sum(pearsonR^2, na.rm = TRUE)/(nrow(tbmsdata_sp_sum) - length(f1_tbms$par) - length(Si.out))
  # Variance-covariance matrix
  vcov <- solve(f1_tbms$hessian)
  
  # Data frame with net growth rates with ci (scaled by chat), chat, ci width and significance
  out_tbms_net_spp <- data.table(SPECIES = spp,
                                COMMON_NAME = tbmsdata_sp$COMMON_NAME[1],
                                ngr = ngr,
                                ngr_lower = ngr - qnorm(.975)*sqrt(sum(vcov))*sqrt(chat),
                                ngr_upper = ngr + qnorm(.975)*sqrt(sum(vcov))*sqrt(chat),
                                chat = chat,
                                vcov_sum = sum(vcov))
  out_tbms_net_spp[, ngr_ci_width := ngr_upper - ngr_lower]
  out_tbms_net_spp[, ngr_sig := ifelse(ngr_upper*ngr_lower > 0, "*", "")]
  out_tbms_net_spp$nsites <- nsites # total number of sites (not unique)
  out_tbms_net_spp$nvisits <- sum(tbmsdata_sp_sum$nV) # total number of visits
  out_tbms_net_spp$mean_count <- mean(tbmsdata_sp_sum$eCount) # mean count
  out_tbms_net_spp$std_count <- std.error(tbmsdata_sp_sum$eCount) # standard error of mean count
  out_tbms_net <- rbind(out_tbms_net, out_tbms_net_spp)
  
  
  # Save annual growth rate outputs to a data frame
  out_tbms_annual <- rbind(out_tbms_annual,
                          data.table(SPECIES = spp,
                                     COMMON_NAME = tbmsdata_sp$COMMON_NAME[1],
                                     YEAR = head(years, -1),
                                     agr = f1_tbms$par,
                                     agr_se = sqrt(diag(solve(f1_tbms$hessian))),
                                     chat = chat))
  # Save variance-covariance matrix
  vcov_tbms[[tbmsdata_sp$COMMON_NAME[1]]] <- solve(f1_tbms$hessian)
  # Save all outputs
  glm_tbms[[tbmsdata_sp$COMMON_NAME[1]]] <- list(SPECIES = spp,
                                               COMMON_NAME = tbmsdata_sp$COMMON_NAME[1],
                                               f1_tbms = f1_tbms,
                                               tbmsdata_sp_sum = tbmsdata_sp_sum,
                                               Counts = Counts,
                                               Visits = Visits,
                                               Si.out = Si.out,
                                               Fitted = Fitted,
                                               pearsonR = pearsonR)
}




saveRDS(out_tbms_net, paste0("Output/Growth rate analysis/tBMS/tBMS_gra_net_", min(years), "_", max(years), ".rds"))
saveRDS(out_tbms_annual, paste0("Output/Growth rate analysis/tBMS/tBMS_gra_annual_", min(years), "_", max(years), ".rds"))
saveRDS(vcov_tbms, paste0("Output/Growth rate analysis/tBMS/tBMS_gra_vcov_", min(years), "_", max(years), ".rds"))
saveRDS(glm_tbms, paste0("Output/Growth rate analysis/tBMS/tBMS_gra_output_", min(years), "_", max(years), ".rds"))




## 2. Growth rate model for WCBS only

out_wcbs_net <- out_wcbs_annual <- vcov_wcbs <- glm_wcbs <- NULL
for(spp in splist){
  cat(spp,"\n")
  # wcbs data for species spp
  wcbsdata_sp <- wcbsdata[SPECIES == spp]
  if(nrow(wcbsdata_sp) == 0) next() # probably not needed here 
  
  # Sum counts over multiple visits
  wcbsdata_sp_sum <- wcbsdata_sp[, .(eCount = sum(COUNT), nV = length(COUNT)), by = .(SITENO, YEAR)]
  # eCount = total counts over visits at each site and year
  # nV = total number of visits at each site and year

  colnames(wcbsdata_sp_sum)[1:2] <- c("Site","Year")
  
  # Set up for growth-rate GLM model as in Roy et al. (2015) WCBS paper
  wcbsdata_sp_sum <- wcbsdata_sp_sum[order(wcbsdata_sp_sum$Site, wcbsdata_sp_sum$Year, decreasing=FALSE),]
  
  nsites <- length(unique(wcbsdata_sp_sum$Site))
  nyears <- length(unique(wcbsdata_sp_sum$Year))
  years <- sort(unique(wcbsdata_sp_sum$Year))
  Counts <- as.matrix(reshape(wcbsdata_sp_sum[,c("Site","Year","eCount")],
                              direction = "wide", timevar = "Year", idvar = "Site")[,paste("eCount", years, sep = "."), with = FALSE], ncol = nyears)
  Visits <- as.matrix(reshape(wcbsdata_sp_sum[,c("Site","Year","nV")],
                              direction = "wide", timevar = "Year", idvar = "Site")[,paste("nV", years, sep = "."), with = FALSE], ncol = nyears)
  
  parm <- rep(0, nyears-1)
  cat("GLM for wcbs for",spp,"starting at",date(),"\n")
  stc <- Sys.time()
  # Fit model with concentrated likelihood
  f1_wcbs <- try(optim(par = parm,
                       fn = ll_func2,
                       l.type = "conc",
                       Counts = Counts,
                       Visits = Visits,
                       method = "BFGS",
                       hessian = TRUE,
                       control = list(maxit = 2000)),
                 silent = TRUE)
  etc <- Sys.time()
  
  cat("GLM for wcbs for", spp, "done at", date(), "\n")
  
  # Net change estimates up to each year, in matrix form for later calculations
  gr.estS <- matrix(cumsum(c(0, f1_wcbs$par)), nrow = nsites, ncol = nyears, byrow = TRUE) 
  # Estimate site indices
  Si.out <- log(apply(Counts, 1, sum, na.rm = TRUE)/apply(Visits * exp(gr.estS), 1, sum, na.rm = TRUE))
  # Estimated total counts
  Fitted <- exp(gr.estS + Si.out + log(Visits))
  # Pearson residuals
  pearsonR <- (Counts - Fitted)/sqrt(Fitted)
  
  # Calculate net change estimate
  ngr <- sum(f1_wcbs$par)
  # Pearson chi-square statistic/df
  chat <- sum(pearsonR^2, na.rm = TRUE)/(nrow(wcbsdata_sp_sum) - length(f1_wcbs$par) - length(Si.out))
  # Variance-covariance matrix
  vcov <- solve(f1_wcbs$hessian)
  
  # Data frame with net growth rates with ci (scaled by chat), chat, ci width and significance
  out_wcbs_net_spp <- data.table(SPECIES = spp,
                                 COMMON_NAME = wcbsdata_sp$COMMON_NAME[1],
                                 ngr = ngr,
                                 ngr_lower = ngr - qnorm(.975)*sqrt(sum(vcov))*sqrt(chat),
                                 ngr_upper = ngr + qnorm(.975)*sqrt(sum(vcov))*sqrt(chat),
                                 chat = chat,
                                 vcov_sum = sum(vcov))
  out_wcbs_net_spp[, ngr_ci_width := ngr_upper - ngr_lower]
  out_wcbs_net_spp[, ngr_sig := ifelse(ngr_upper*ngr_lower > 0, "*", "")]
  out_wcbs_net_spp$nsites <- nsites # total number of sites (not unique)
  out_wcbs_net_spp$nvisits <- sum(wcbsdata_sp_sum$nV[wcbsdata_sp_sum$eCount>0 & wcbsdata_sp_sum$Year<=2013]) # total number of visits
  out_wcbs_net_spp$mean_count <- mean(wcbsdata_sp_sum$eCount[wcbsdata_sp_sum$Year<=2013]) # mean count
  out_wcbs_net_spp$std_count <- std.error(wcbsdata_sp_sum$eCount) # standard error of mean count
  
  out_wcbs_net <- rbind(out_wcbs_net, out_wcbs_net_spp)
  
  
  # Save annual growth rate outputs to a data frame
  out_wcbs_annual <- rbind(out_wcbs_annual,
                           data.table(SPECIES = spp,
                                      COMMON_NAME = wcbsdata_sp$COMMON_NAME[1],
                                      YEAR = head(years, -1),
                                      agr = f1_wcbs$par,
                                      agr_se = sqrt(diag(solve(f1_wcbs$hessian))),
                                      chat = chat))
  # Save variance-covariance matrix
  vcov_wcbs[[wcbsdata_sp$COMMON_NAME[1]]] <- solve(f1_wcbs$hessian)
  # Save all outputs
  glm_wcbs[[wcbsdata_sp$COMMON_NAME[1]]] <- list(SPECIES = spp,
                                                 COMMON_NAME = wcbsdata_sp$COMMON_NAME[1],
                                                 f1_wcbs = f1_wcbs,
                                                 wcbsdata_sp_sum = wcbsdata_sp_sum,
                                                 Counts = Counts,
                                                 Visits = Visits,
                                                 Si.out = Si.out,
                                                 Fitted = Fitted,
                                                 pearsonR = pearsonR)
}




saveRDS(out_wcbs_net, paste0("Output/Growth rate analysis/WCBS/WCBS_gra_net_", min(years), "_", max(years), ".rds"))
saveRDS(out_wcbs_annual, paste0("Output/Growth rate analysis/WCBS/WCBS_gra_annual_", min(years), "_", max(years), ".rds"))
saveRDS(vcov_wcbs, paste0("Output/Growth rate analysis/WCBS/WCBS_gra_vcov_", min(years), "_", max(years), ".rds"))
saveRDS(glm_wcbs, paste0("Output/Growth rate analysis/WCBS/WCBS_gra_output_", min(years), "_", max(years), ".rds"))



####################################

## 10-year rolling growth rate model

## 1. Growth rate model for tBMS only

year_list <- 2009:2015

for(year in year_list){
  out_tbms_net <- NULL
  print(year)
  start_year <- year
  end_year <- as.integer(year+9)
  tbmsdata_yr <- tbmsdata[tbmsdata$YEAR>=start_year & tbmsdata$YEAR<=end_year,]
  
for(spp in splist$COMMON_NAME){
  cat(spp,"\n")
  # tBMS data for species spp
  tbmsdata_sp <- tbmsdata_yr[COMMON_NAME == spp]
  
  if(nrow(tbmsdata_sp) == 0) next() # probably not needed here 
  
  # Sum counts over multiple visits
  tbmsdata_sp_sum <- tbmsdata_sp[, .(eCount = sum(COUNT), nV = length(COUNT)), by = .(SITENO, YEAR)]
  # eCount = total counts over visits at each site and year
  # nV = total number of visits at each site and year
  
  colnames(tbmsdata_sp_sum)[1:2] <- c("Site","Year")
  
  # Set up for growth-rate GLM model as in Roy et al. (2015) WCBS paper
  tbmsdata_sp_sum <- tbmsdata_sp_sum[order(tbmsdata_sp_sum$Site, tbmsdata_sp_sum$Year, decreasing=FALSE),]
  
  nsites <- length(unique(tbmsdata_sp_sum$Site))
  nyears <- length(unique(tbmsdata_sp_sum$Year))
  years <- sort(unique(tbmsdata_sp_sum$Year))
  Counts <- as.matrix(reshape(tbmsdata_sp_sum[,c("Site","Year","eCount")],
                              direction = "wide", timevar = "Year", idvar = "Site")[,paste("eCount", years, sep = "."), with = FALSE], ncol = nyears)
  Visits <- as.matrix(reshape(tbmsdata_sp_sum[,c("Site","Year","nV")],
                              direction = "wide", timevar = "Year", idvar = "Site")[,paste("nV", years, sep = "."), with = FALSE], ncol = nyears)
  
  parm <- rep(0, nyears-1)
  cat("GLM for tBMS for",spp,"starting at",date(),"\n")
  stc <- Sys.time()
  # Fit model with concentrated likelihood
  f1_tbms <- try(optim(par = parm,
                       fn = ll_func2,
                       l.type = "conc",
                       Counts = Counts,
                       Visits = Visits,
                       method = "BFGS",
                       hessian = TRUE,
                       control = list(maxit = 2000)),
                 silent = TRUE)
  etc <- Sys.time()
  
  cat("GLM for tBMS for", spp, "done at", date(), "\n")
  # Par = annual growth rates 
  # Hessian = square matrix of second-order partial derivatives of a function - 
  # provides info on the curvature of the likelihood function. Inverse of the 
  # hessian often used to esimate variance of maximum likelihood estimators 
  # for remaining parameters 
  
  # Net change estimates up to each year, in matrix form for later calculations
  gr.estS <- matrix(cumsum(c(0, f1_tbms$par)), nrow = nsites, ncol = nyears, byrow = TRUE) 
  # Estimate site indices
  Si.out <- log(apply(Counts, 1, sum, na.rm = TRUE)/apply(Visits * exp(gr.estS), 1, sum, na.rm = TRUE))
  # Estimated total counts
  Fitted <- exp(gr.estS + Si.out + log(Visits))
  # Pearson residuals
  pearsonR <- (Counts - Fitted)/sqrt(Fitted)
  
  # Calculate net change estimate
  ngr <- sum(f1_tbms$par)
  # Pearson chi-square statistic/df
  chat <- sum(pearsonR^2, na.rm = TRUE)/(nrow(tbmsdata_sp_sum) - length(f1_tbms$par) - length(Si.out))
  # Variance-covariance matrix
  vcov <- solve(f1_tbms$hessian)
  
  # Data frame with net growth rates with ci (scaled by chat), chat, ci width and significance
  out_tbms_net_spp <- data.table(SPECIES = spp,
                                 COMMON_NAME = tbmsdata_sp$COMMON_NAME[1],
                                 ngr = ngr,
                                 ngr_lower = ngr - qnorm(.975)*sqrt(sum(vcov))*sqrt(chat),
                                 ngr_upper = ngr + qnorm(.975)*sqrt(sum(vcov))*sqrt(chat),
                                 chat = chat,
                                 vcov_sum = sum(vcov))
  out_tbms_net_spp[, ngr_ci_width := ngr_upper - ngr_lower]
  out_tbms_net_spp[, ngr_sig := ifelse(ngr_upper*ngr_lower > 0, "*", "")]
  out_tbms_net_spp$nsites <- nsites # total number of sites (not unique)
  out_tbms_net_spp$nvisits <- sum(tbmsdata_sp_sum$nV) # total number of visits
  out_tbms_net_spp$mean_count <- mean(tbmsdata_sp_sum$eCount) # mean count
  out_tbms_net_spp$std_count <- std.error(tbmsdata_sp_sum$eCount) # standard error of mean count
  out_tbms_net <- rbind(out_tbms_net, out_tbms_net_spp)
  
  
}
  saveRDS(out_tbms_net, paste0("Output/Growth rate analysis/tBMS/tBMS_gra_net_", start_year, "_", end_year, ".rds"))
  
}





## 2. Growth rate model for WCBS only

year_list <- 2009:2015

for(year in year_list){
  out_wcbs_net <- NULL
  print(year)
  start_year <- year
  end_year <- as.integer(year+9)
  wcbsdata_yr <- wcbsdata[wcbsdata$YEAR>=start_year & wcbsdata$YEAR<=end_year,]
  
for(spp in splist$COMMON_NAME){
  cat(spp,"\n")
  # wcbs data for species spp
  wcbsdata_sp <- wcbsdata_yr[COMMON_NAME == spp]
  if(nrow(wcbsdata_sp) == 0) next() # probably not needed here 
  
  # Sum counts over multiple visits
  wcbsdata_sp_sum <- wcbsdata_sp[, .(eCount = sum(COUNT), nV = length(COUNT)), by = .(SITENO, YEAR)]
  # eCount = total counts over visits at each site and year
  # nV = total number of visits at each site and year
  
  colnames(wcbsdata_sp_sum)[1:2] <- c("Site","Year")
  
  # Set up for growth-rate GLM model as in Roy et al. (2015) WCBS paper
  wcbsdata_sp_sum <- wcbsdata_sp_sum[order(wcbsdata_sp_sum$Site, wcbsdata_sp_sum$Year, decreasing=FALSE),]
  
  nsites <- length(unique(wcbsdata_sp_sum$Site))
  nyears <- length(unique(wcbsdata_sp_sum$Year))
  years <- sort(unique(wcbsdata_sp_sum$Year))
  Counts <- as.matrix(reshape(wcbsdata_sp_sum[,c("Site","Year","eCount")],
                              direction = "wide", timevar = "Year", idvar = "Site")[,paste("eCount", years, sep = "."), with = FALSE], ncol = nyears)
  Visits <- as.matrix(reshape(wcbsdata_sp_sum[,c("Site","Year","nV")],
                              direction = "wide", timevar = "Year", idvar = "Site")[,paste("nV", years, sep = "."), with = FALSE], ncol = nyears)
  
  parm <- rep(0, nyears-1)
  cat("GLM for wcbs for",spp,"starting at",date(),"\n")
  stc <- Sys.time()
  # Fit model with concentrated likelihood
  f1_wcbs <- try(optim(par = parm,
                        fn = ll_func2,
                        l.type = "conc",
                        Counts = Counts,
                        Visits = Visits,
                        method = "BFGS",
                        hessian = TRUE,
                        control = list(maxit = 2000)),
                  silent = TRUE)
  etc <- Sys.time()
  
  cat("GLM for wcbs for", spp, "done at", date(), "\n")
  
  # Net change estimates up to each year, in matrix form for later calculations
  gr.estS <- matrix(cumsum(c(0, f1_wcbs$par)), nrow = nsites, ncol = nyears, byrow = TRUE) 
  # Estimate site indices
  Si.out <- log(apply(Counts, 1, sum, na.rm = TRUE)/apply(Visits * exp(gr.estS), 1, sum, na.rm = TRUE))
  # Estimated total counts
  Fitted <- exp(gr.estS + Si.out + log(Visits))
  # Pearson residuals
  pearsonR <- (Counts - Fitted)/sqrt(Fitted)
  
  # Calculate net change estimate
  ngr <- sum(f1_wcbs$par)
  # Pearson chi-square statistic/df
  chat <- sum(pearsonR^2, na.rm = TRUE)/(nrow(wcbsdata_sp_sum) - length(f1_wcbs$par) - length(Si.out))
  # Variance-covariance matrix
  vcov <- solve(f1_wcbs$hessian)
  
  # Data frame with net growth rates with ci (scaled by chat), chat, ci width and significance
  out_wcbs_net_spp <- data.table(SPECIES = spp,
                                  COMMON_NAME = wcbsdata_sp$COMMON_NAME[1],
                                  ngr = ngr,
                                  ngr_lower = ngr - qnorm(.975)*sqrt(sum(vcov))*sqrt(chat),
                                  ngr_upper = ngr + qnorm(.975)*sqrt(sum(vcov))*sqrt(chat),
                                  chat = chat,
                                  vcov_sum = sum(vcov))
  out_wcbs_net_spp[, ngr_ci_width := ngr_upper - ngr_lower]
  out_wcbs_net_spp[, ngr_sig := ifelse(ngr_upper*ngr_lower > 0, "*", "")]
  out_wcbs_net_spp$nsites <- nsites # total number of sites (not unique)
  out_wcbs_net_spp$nvisits <- sum(wcbsdata_sp_sum$nV[wcbsdata_sp_sum$eCount>0 & wcbsdata_sp_sum$Year<=2013]) # total number of visits
  out_wcbs_net_spp$mean_count <- mean(wcbsdata_sp_sum$eCount[wcbsdata_sp_sum$Year<=2013]) # mean count
  out_wcbs_net_spp$std_count <- std.error(wcbsdata_sp_sum$eCount) # standard error of mean count
  
  out_wcbs_net <- rbind(out_wcbs_net, out_wcbs_net_spp)

}

saveRDS(out_wcbs_net, paste0("Output/Growth rate analysis/WCBS/WCBS_gra_net_", min(years), "_", max(years), ".rds"))

}



