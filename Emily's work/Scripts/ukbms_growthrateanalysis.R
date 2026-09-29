
# Growth rate analysis for UKBMS data, limited to BBC official date period each year
# THOUGHT - could do for July and August as a whole too??
# WHAT TO DO ABOUT WCBS?

library(speedglm)
library(Matrix)
library(plyr)
library(data.table)


# Likelihood function - growth rate model as in Roy et al WCBS paper, but using a concentrated likelihood
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
    Si <- log(apply(Counts, 1, sum, na.rm = TRUE)/apply(Visits*exp(gr.estS), 1, sum, na.rm = TRUE))
  }

  # Poisson expectaction
  mu.est <- exp(gr.estS + Si + log(Visits))

  # Likelihood formulation
  llik <- dpois(Counts, lambda = mu.est, log = TRUE)

  -1*sum(llik, na.rm = TRUE)
}

num <- function(x){as.numeric(as.character(x))}

# List of 18 BBC butterfly species
specieslist <- read.table("Data/BBC species list.txt", header = TRUE)

# BBC official date periods per year
bbcdates <- read.csv("Data/BBC_official_date_periods.csv")
colnames(bbcdates)[colnames(bbcdates) == "Year"] <- "YEAR"

# Read in UKBMS data
bmsdata <- readRDS("../../Data/UKBMS_data/Data to 2021/Formatted/sp_weeklycount_1976-2021_zerofilled.rds")
# Filter to BBC years and July/August only
bmsdata <- bmsdata[YEAR >= 2011 &
                MONTH %in% 7:8]
# Filter out WCBS or not
excludeWCBS <- FALSE
if(excludeWCBS == TRUE)
  bmsdata <- bmsdata[SITENO < 50000]

excludeUKBMS <- TRUE
if(excludeUKBMS == TRUE)
  bmsdata <- bmsdata[SITENO >= 50000]

# Use data from official BBC period (TRUE) or all of July/August (FALSE)
BBCdateperiods <- TRUE

out_bms_net <- out_bms_annual <- vcov_bms <- glm_bms <- NULL
for(spp in specieslist$SPECIES){
  cat(spp,"\n")
  # UKBMS data for species spp
  bmsdata_sp <- bmsdata[SPECIES == spp]
  if(nrow(bmsdata_sp) == 0) next()
  # Limit to BBC date periods
  bmsdata_sp <- merge(bmsdata_sp, bbcdates, by = "YEAR")
  if(BBCdateperiods)
    bmsdata_sp <- bmsdata_sp[(DAY >= startday & MONTH == 7) | (DAY <= endday & MONTH == 8)]

  # Sum counts over multiple visits
  bmsdata_sp_sum <- bmsdata_sp[, .(eCount = sum(COUNT), nV = length(COUNT)), by = .(SITENO, YEAR)]

  colnames(bmsdata_sp_sum)[1:2] <- c("Site","Year")

  # Set up for growth-rate GLM model as in Roy et al. (2015) WCBS paper
    bmsdata_sp_sum <- bmsdata_sp_sum[order(bmsdata_sp_sum$Site, bmsdata_sp_sum$Year, decreasing=FALSE),]
    nsites <- length(unique(bmsdata_sp_sum$Site))
    nyears <- length(unique(bmsdata_sp_sum$Year))
    years <- sort(unique(bmsdata_sp_sum$Year))
    Counts <- as.matrix(reshape(bmsdata_sp_sum[,c("Site","Year","eCount")],
                                direction = "wide", timevar = "Year", idvar = "Site")[,paste("eCount", years, sep = "."), with = FALSE], ncol = nyears)
    Visits <- as.matrix(reshape(bmsdata_sp_sum[,c("Site","Year","nV")],
                                direction = "wide", timevar = "Year", idvar = "Site")[,paste("nV", years, sep = "."), with = FALSE], ncol = nyears)

    parm <- rep(0, nyears-1)
    cat("GLM for UKBMS for",spp,"starting at",date(),"\n")
    stc <- Sys.time()
    # Fit model with concentrated likelihood
    f1_bms <- try(optim(par = parm,
                        fn = ll_func2,
                        l.type = "conc",
                        Counts = Counts,
                        Visits = Visits,
                        method = "BFGS",
                        hessian = TRUE,
                        control = list(maxit = 2000)),
                  silent = TRUE)
    etc <- Sys.time()

    cat("GLM for UKBMS for", spp, "done at", date(), "\n")
    # Net change estimates up to each year, in matrix form for later calculations
    gr.estS <- matrix(cumsum(c(0, f1_bms$par)), nrow = nsites, ncol = nyears, byrow = TRUE)
    # Estimate site indices
    Si.out <- log(apply(Counts, 1, sum, na.rm = TRUE)/apply(Visits * exp(gr.estS), 1, sum, na.rm = TRUE))
    # Estimated total counts
    Fitted <- exp(gr.estS + Si.out + log(Visits))
    # Pearson residuals
    pearsonR <- (Counts - Fitted)/sqrt(Fitted)

    # Calculate net change estimate
    ngr <- sum(f1_bms$par)
    # Pearson chi-square statistic/df
    chat <- sum(pearsonR^2, na.rm = TRUE)/(nrow(bmsdata_sp_sum) - length(f1_bms$par) - length(Si.out))
    # Variance-covariance matrix
    vcov <- solve(f1_bms$hessian)

    # Data frame with net growth rates with ci (scaled by chat), chat, ci width and significance
    out_bms_net_spp <- data.table(SPECIES = spp,
                                  COMMON_NAME = bmsdata_sp$COMMON_NAME[1],
                                  ngr = ngr,
                                  ngr_lower = ngr - qnorm(.975)*sqrt(sum(vcov))*sqrt(chat),
                                  ngr_upper = ngr + qnorm(.975)*sqrt(sum(vcov))*sqrt(chat),
                                  chat = chat,
                                  vcov_sum = sum(vcov))
    out_bms_net_spp[, ngr_ci_width := ngr_upper - ngr_lower]
    out_bms_net_spp[, ngr_sig := ifelse(ngr_upper*ngr_lower > 0, "*", "")]

    out_bms_net <- rbind(out_bms_net, out_bms_net_spp)


     # Save annual growth rate outputs to a data frame
     out_bms_annual <- rbind(out_bms_annual,
                           data.table(SPECIES = spp,
                                      COMMON_NAME = bmsdata_sp$COMMON_NAME[1],
                                      YEAR = head(years, -1),
                                      agr = f1_bms$par,
                                      agr_se = sqrt(diag(solve(f1_bms$hessian))),
                                      chat = chat))
    # Save variance-covariance matrix
    vcov_bms[[bmsdata_sp$COMMON_NAME[1]]] <- solve(f1_bms$hessian)
    # Save all outputs
    glm_bms[[bmsdata_sp$COMMON_NAME[1]]] <- list(SPECIES = spp,
                                                 COMMON_NAME = bmsdata_sp$COMMON_NAME[1],
                                                 f1_bms = f1_bms,
                                                 bmsdata_sp_sum = bmsdata_sp_sum,
                                                 Counts = Counts,
                                                 Visits = Visits,
                                                 Si.out = Si.out,
                                                 Fitted = Fitted,
                                                 pearsonR = pearsonR)
  }




saveRDS(out_bms_net, paste0("Output/growthrateanalysis/UKBMS/ukbms_gra_net_", min(years), "_", max(years),  if(!BBCdateperiods) "_JULAUG", if(excludeWCBS == FALSE) "_incWCBS", if(excludeUKBMS == TRUE)"_WCBSonly" , ".rds"))
saveRDS(out_bms_annual, paste0("Output/growthrateanalysis/UKBMS/ukbms_gra_annual_", min(years), "_", max(years), if(!BBCdateperiods) "_JULAUG", if(excludeWCBS == FALSE) "_incWCBS", if(excludeUKBMS == TRUE)"_WCBSonly" , ".rds"))
saveRDS(vcov_bms, paste0("Output/growthrateanalysis/UKBMS/ukbms_gra_vcov_", min(years), "_", max(years), if(!BBCdateperiods) "_JULAUG", if(excludeWCBS == FALSE) "_incWCBS",if(excludeUKBMS == TRUE)"_WCBSonly" ,  ".rds"))
saveRDS(glm_bms, paste0("Output/growthrateanalysis/UKBMS/ukbms_gra_output_", min(years), "_", max(years), if(!BBCdateperiods) "_JULAUG", if(excludeWCBS == FALSE) "_incWCBS",if(excludeUKBMS == TRUE)"_WCBSonly" ,  ".rds"))









