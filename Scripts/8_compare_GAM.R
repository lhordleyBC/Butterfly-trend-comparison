## Use the tBMS and WCBS annual indices for each of the 23 species to compare between the surveys
## Includes calculating confidence intervals from bootstrap iterations
## Comparing overall trends

# load packages
library(rbms)
library(dplyr)
library(ggplot2)
library(data.table)
library(raster)
library(ggpubr)
library(ggrepel)

# Read in data
wcbs_annual_indices <- readRDS("Output/GAM abundance analysis/Annual_indices_WCBS.rds")
tbms_annual_indices <- readRDS("Output/GAM abundance analysis/Annual_indices_tBMS.rds")

sp_names <- read.csv("Data/UKBMS/Species_list_GAM.csv", header=TRUE)
tbms_annual_indices <- merge(tbms_annual_indices, sp_names, by.x="sp", by.y="SPECIES")
wcbs_annual_indices <- merge(wcbs_annual_indices, sp_names, by.x="sp", by.y="SPECIES")

# calculate confidence intervals around each annual index for each species using bootstrapped values
# tBMS
boot_CIs <- tbms_annual_indices %>% group_by(M_YEAR, COMMON_NAME) %>% filter(BOOTi>0) %>% summarise(lowerCI=quantile(LCI, 0.025),upperCI=quantile(LCI, 0.975))
true_values <- tbms_annual_indices[tbms_annual_indices$BOOTi==0,] # filter to true value
tbms_annual_indices2 <-  merge(true_values, boot_CIs, by=c("M_YEAR", "COMMON_NAME"))
# WCBS
boot_CIs <- wcbs_annual_indices %>% group_by(M_YEAR, COMMON_NAME) %>% filter(BOOTi>0) %>% summarise(lowerCI=quantile(LCI, 0.025),upperCI=quantile(LCI, 0.975))
true_values <- wcbs_annual_indices[wcbs_annual_indices$BOOTi==0,] # filter to true value
wcbs_annual_indices2 <-  merge(true_values, boot_CIs, by=c("M_YEAR", "COMMON_NAME"))

tbms_annual_indices2$SURVEY <- "tBMS"
wcbs_annual_indices2$SURVEY <- "WCBS"

annual_indices2 <- rbind(tbms_annual_indices2, wcbs_annual_indices2)
# save file
write.csv(annual_indices2, file="Output/GAM abundance analysis/Annual_indices_tBMS_WCBS.csv", row.names=FALSE)

annual_indices2 <- read.csv("Output/GAM abundance analysis/Annual_indices_tBMS_WCBS.csv", header=TRUE)

# plot annual indices for each species
annual_indices_p <- ggplot(data=annual_indices2, aes(x=M_YEAR, y=LCI, group=SURVEY))+
  geom_point(size=0.5, aes(colour=SURVEY), position=position_dodge(width=0.5))+
  geom_line(lwd=0.3, aes(colour=SURVEY))+
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI, colour=SURVEY), lwd=0.3, position=position_dodge(width=0.5), width=0.3)+
  facet_wrap(~COMMON_NAME, ncol=3, scales="free")+
  theme_minimal()
annual_indices_p
ggsave(annual_indices_p, file="Output/Figures/GAM_annual_indices_23spp.png", height=8, width=7)

# run linear models and back-transform linear model coefficients from log10 scale to growth rates

# loop through each bootstrap iteration to extract regression coefficient
# back-transform all coefficients

### tBMS ###
boots <- unique(tbms_annual_indices$BOOTi)
sp <- unique(tbms_annual_indices$sp)
tbms_growth_rates <- NULL
for(s in sp){
  print(s)
  annual_index_sp <- tbms_annual_indices[tbms_annual_indices$sp==s,]
  for(i in boots){
    
    annual_index_sp_temp <- annual_index_sp[annual_index_sp$BOOTi==i,]
    
    lm_fit <- lm(LCI ~ M_YEAR, data = annual_index_sp_temp)
    # extract and back-transform annual growth rate
    growth_rate <- exp(coef(lm_fit)[2]*2.303)
    perc_change <- 100*(exp(coef(lm_fit)[2]*2.303)-1)
    total_perc_change <- 100*(exp(coef(lm_fit)[2]*2.303)^(max(annual_index_sp_temp$M_YEAR)-min(annual_index_sp_temp$M_YEAR))-1)
    # build data frame
    tbms_growth_rates_t <- data.frame(growth_rate=growth_rate, perc_change=perc_change, total_perc_change=total_perc_change, BOOTi=i, sp=s)
    tbms_growth_rates <- rbind(tbms_growth_rates, tbms_growth_rates_t)
  }
}

# calculate upper and lower CIs from the bootstrapped coefficients (excluding the 'true' one)
bootstrap_CIs <- tbms_growth_rates %>% group_by(sp) %>% filter(BOOTi>0) %>% 
  dplyr::summarise(lowerCI_gr=quantile(growth_rate, 0.025), upperCI_gr=quantile(growth_rate, 0.975),
            lowerCI_perc=quantile(perc_change, 0.025), upperCI_perc=quantile(perc_change, 0.975),
            lowerCI_perc_tot=quantile(total_perc_change, 0.025), upperCI_perc_tot=quantile(total_perc_change, 0.975))
tbms_growth_rates_true <- tbms_growth_rates[tbms_growth_rates$BOOTi==0,]
tbms_growth_rates2 <-  merge(tbms_growth_rates_true, bootstrap_CIs, by="sp")
tbms_growth_rates2$SURVEY <- "tBMS"


### WCBS ###
boots <- unique(wcbs_annual_indices$BOOTi)
sp <- unique(wcbs_annual_indices$sp)
wcbs_growth_rates <- NULL
for(s in sp){
  print(s)
  annual_index_sp <- wcbs_annual_indices[wcbs_annual_indices$sp==s,]
  for(i in boots){
    
    annual_index_sp_temp <- annual_index_sp[annual_index_sp$BOOTi==i,]
    
    lm_fit <- lm(LCI ~ M_YEAR, data = annual_index_sp_temp)
    # extract and back-transform annual growth rate
    growth_rate <- exp(coef(lm_fit)[2]*2.303)
    perc_change <- 100*(exp(coef(lm_fit)[2]*2.303)-1)
    total_perc_change <- 100*(exp(coef(lm_fit)[2]*2.303)^(max(annual_index_sp_temp$M_YEAR)-min(annual_index_sp_temp$M_YEAR))-1)
    # build data frame
    wcbs_growth_rates_t <- data.frame(growth_rate=growth_rate, perc_change=perc_change, total_perc_change=total_perc_change, BOOTi=i, sp=s)
    wcbs_growth_rates <- rbind(wcbs_growth_rates, wcbs_growth_rates_t)
  }
}

# calculate upper and lower CIs from the bootstrapped coefficients (excluding the 'true' one)
bootstrap_CIs <- wcbs_growth_rates %>% group_by(sp) %>% filter(BOOTi>0) %>% 
  summarise(lowerCI_gr=quantile(growth_rate, 0.025),upperCI_gr=quantile(growth_rate, 0.975),
            lowerCI_perc=quantile(perc_change, 0.025), upperCI_perc=quantile(perc_change, 0.975),
            lowerCI_perc_tot=quantile(total_perc_change, 0.025), upperCI_perc_tot=quantile(total_perc_change, 0.975))
wcbs_growth_rates_true <- wcbs_growth_rates[wcbs_growth_rates$BOOTi==0,]
wcbs_growth_rates2 <-  merge(wcbs_growth_rates_true, bootstrap_CIs, by="sp")
wcbs_growth_rates2$SURVEY <- "WCBS"


growth_rates <- rbind(tbms_growth_rates2,wcbs_growth_rates2)
growth_rates <- merge(growth_rates, sp_names, by.x="sp", by.y="SPECIES")

# change from long to wide to plot
growth_rates2 <- reshape(growth_rates, idvar = c("COMMON_NAME", "sp", "BOOTi"), timevar = "SURVEY", direction = "wide")
growth_rates2$growth_rate_diff <- growth_rates2$growth_rate.tBMS - growth_rates2$growth_rate.WCBS

# save file
write.csv(growth_rates2, file="Output/GAM abundance analysis/Abundance_trends_tBMS_WCBS_GAM.csv", row.names=FALSE)
growth_rates2 <- read.csv("Output/GAM abundance analysis/Abundance_trends_tBMS_WCBS_GAM.csv", header=TRUE)

growth_rate_p <- ggplot(growth_rates2, aes(growth_rate.WCBS, growth_rate.tBMS))+
  geom_hline(yintercept = 1, color = "black", linetype="dotted")+
  geom_vline(xintercept = 1, color = "black", linetype="dotted")+
  geom_abline(linetype="dashed", color = "black")+
  geom_point(size = 2)+
  geom_text_repel(aes(label = COMMON_NAME), color="grey", size = 5)+
  geom_errorbar(aes(ymin = lowerCI_gr.tBMS, ymax = upperCI_gr.tBMS), width=0.005)+
  geom_errorbar(aes(xmin = lowerCI_gr.WCBS, xmax = upperCI_gr.WCBS), width=0.005)+
  theme_bw()+
  #ylim(-3,1.2)+
  #xlim(-3,1.2)+
  theme(text = element_text(size = 20)) +
  xlab("WCBS abundance trend 2009-2024")+
  ylab("tBMS abundance trend 2009-2024")
growth_rate_p
ggsave(growth_rate_p, file="Output/Figures/GAM_results_gr2.png", height=20, width=20, units="cm")

cor.test(growth_rates2$growth_rate.tBMS, growth_rates2$growth_rate.WCBS) # r=0.95, p<0.001

# compare width of confidence intervals between surveys 
growth_rates2$ciwidth_tBMS <- growth_rates2$upperCI_gr.tBMS - growth_rates2$lowerCI_gr.tBMS
growth_rates2$ciwidth_WCBS <- growth_rates2$upperCI_gr.WCBS - growth_rates2$lowerCI_gr.WCBS

mean(growth_rates2$ciwidth_tBMS)
mean(growth_rates2$ciwidth_WCBS)
sd(growth_rates2$ciwidth_tBMS)
sd(growth_rates2$ciwidth_WCBS)

t.test(growth_rates2$ciwidth_tBMS, growth_rates2$ciwidth_WCBS, paired = TRUE)


#

## calculate difference in GR between surveys + CIs for each species
tbms_growth_rates$SURVEY <- "tBMS"
wcbs_growth_rates$SURVEY <- "WCBS"
growth_rates3 <- rbind(tbms_growth_rates,wcbs_growth_rates)
# long to wide
growth_rates_diff <- reshape(growth_rates3, idvar = c("sp", "BOOTi"), timevar = "SURVEY", direction = "wide")
growth_rates_diff <- growth_rates_diff %>% group_by(sp, BOOTi) %>% summarise(gr_diff=growth_rate.tBMS-growth_rate.WCBS, 
                     perc_diff=perc_change.tBMS-perc_change.WCBS, tot_perc_diff=total_perc_change.tBMS-total_perc_change.WCBS)
growth_rates_diff2 <- growth_rates_diff %>% group_by(sp) %>% filter(BOOTi>0) %>% summarise(lowerCI_gr_diff=quantile(gr_diff, 0.025),upperCI_gr_diff=quantile(gr_diff, 0.975),
                      lowerCI_perc_diff=quantile(perc_diff, 0.025), upperCI_perc_diff=quantile(perc_diff, 0.975),
                      lowerCI_perc_tot_diff=quantile(tot_perc_diff, 0.025), upperCI_perc_tot_diff=quantile(tot_perc_diff, 0.975))
growth_rates_diff_true <- growth_rates_diff[growth_rates_diff$BOOTi==0,]
growth_rates_diff3 <-  merge(growth_rates_diff_true, growth_rates_diff2, by="sp")
growth_rates_diff3 <- merge(growth_rates_diff3, sp_names, by.x="sp", by.y="SPECIES")

growth_rates_diff3$significance <- ifelse(growth_rates_diff3$lowerCI_gr_diff<0 & growth_rates_diff3$upperCI_gr_diff<0, "yes", "no")

## plot of difference between schemes for growth rate
growth_rates_diff3$COMMON_NAME = with(growth_rates_diff3, reorder(COMMON_NAME, gr_diff))

gr_diff_p <- ggplot(growth_rates_diff3, aes(x=gr_diff, y=COMMON_NAME))+
  geom_point(aes(colour=significance))+
  geom_errorbar(aes(xmin = lowerCI_gr_diff, xmax = upperCI_gr_diff, colour=significance))+
  scale_colour_manual(values = c("darkgrey", "black")) +
  theme_bw()+
  geom_vline(xintercept = 0, linetype="dashed", color = "black")+
  theme(text = element_text(size = 20)) +
  ylab("")+
  xlab("Abundance trend \ndifference (tBMS-WCBS)")+
  guides(colour="none")
# positive = tBMS more positive trend than WCBS
# negative = WCBS more positive trend than tBMS
gr_diff_p
ggsave(gr_diff_p, file="Output/Figures/GAM_results_gr1.png", height=25, width=15, units="cm")

# species not overlapping 0:
# small white
# comma
# meadow brown
# holly blue
# all negative => WCBS trends are more positive than tBMS // tBMS trends are more negative than WCBS

## put graphs together for growth rate

gam_results2 <- ggarrange(growth_rate_p, gr_diff_p, labels=c("(a)", "(b)"), 
                         font.label=list(color="black",size=20))
gam_results2

ggsave(gam_results2, file="Output/Figures/Figure4.png", height=18, width=32, units="cm")



# ########################################################
# #### calculate correlation and plot graphs for each 10-year rolling window from 2009-2015
# 
# # Read in data
# wcbs_annual_indices <- readRDS("Output/GAM abundance analysis/Annual_indices_WCBS.rds")
# tbms_annual_indices <- readRDS("Output/GAM abundance analysis/Annual_indices_tBMS.rds")
# 
# sp_names <- read.csv("Data/UKBMS/Species_list_GAM.csv", header=TRUE)
# tbms_annual_indices <- merge(tbms_annual_indices, sp_names, by.x="sp", by.y="SPECIES")
# wcbs_annual_indices <- merge(wcbs_annual_indices, sp_names, by.x="sp", by.y="SPECIES")
# 
# # run linear models and back-transform linear model coefficients from log10 scale to growth rates
# 
# # loop through each bootstrap iteration to extract regression coefficient
# # back-transform all coefficients
# 
# ### tBMS ###
# boots <- unique(tbms_annual_indices$BOOTi)
# sp <- unique(tbms_annual_indices$sp)
# year_list <- 2009:2015
# tbms_growth_rates <- NULL
# 
# for(year in year_list){
#   print(year)
#   start_year <- year
#   end_year <- as.integer(year+9)
#   tbms_annual_indices_yr <- tbms_annual_indices[tbms_annual_indices$M_YEAR>=start_year & tbms_annual_indices$M_YEAR<=end_year,]
#   
# for(s in sp){
#   print(s)
#   annual_index_sp <- tbms_annual_indices_yr[tbms_annual_indices_yr$sp==s,]
#   for(i in boots){
#     
#     annual_index_sp_temp <- annual_index_sp[annual_index_sp$BOOTi==i,]
#     
#     lm_fit <- lm(LCI ~ M_YEAR, data = annual_index_sp_temp)
#     # extract and back-transform annual growth rate
#     growth_rate <- exp(coef(lm_fit)[2]*2.303)
#     perc_change <- 100*(exp(coef(lm_fit)[2]*2.303)-1)
#     total_perc_change <- 100*(exp(coef(lm_fit)[2]*2.303)^(max(annual_index_sp_temp$M_YEAR)-min(annual_index_sp_temp$M_YEAR))-1)
#     # build data frame
#     tbms_growth_rates_t <- data.frame(growth_rate=growth_rate, perc_change=perc_change, total_perc_change=total_perc_change, BOOTi=i, sp=s, years=paste(start_year, "-", end_year, sep=""))
#     tbms_growth_rates <- rbind(tbms_growth_rates, tbms_growth_rates_t)
#     }
#   }
# }
# 
# # calculate upper and lower CIs from the bootstrapped coefficients (excluding the 'true' one)
# bootstrap_CIs <- tbms_growth_rates %>% group_by(sp, years) %>% filter(BOOTi>0) %>% 
#   dplyr::summarise(lowerCI_gr=quantile(growth_rate, 0.025), upperCI_gr=quantile(growth_rate, 0.975),
#                    lowerCI_perc=quantile(perc_change, 0.025), upperCI_perc=quantile(perc_change, 0.975),
#                    lowerCI_perc_tot=quantile(total_perc_change, 0.025), upperCI_perc_tot=quantile(total_perc_change, 0.975))
# tbms_growth_rates_true <- tbms_growth_rates[tbms_growth_rates$BOOTi==0,]
# tbms_growth_rates2 <-  merge(tbms_growth_rates_true, bootstrap_CIs, by=c("sp", "years"))
# tbms_growth_rates2$SURVEY <- "tBMS"
# 
# 
# ### WCBS ###
# boots <- unique(wcbs_annual_indices$BOOTi)
# sp <- unique(wcbs_annual_indices$sp)
# year_list <- 2009:2015
# wcbs_growth_rates <- NULL
# 
# for(year in year_list){
#   print(year)
#   start_year <- year
#   end_year <- as.integer(year+9)
#   wcbs_annual_indices_yr <- wcbs_annual_indices[wcbs_annual_indices$M_YEAR>=start_year & wcbs_annual_indices$M_YEAR<=end_year,]
#   
# for(s in sp){
#   print(s)
#   annual_index_sp <- wcbs_annual_indices_yr[wcbs_annual_indices_yr$sp==s,]
#   for(i in boots){
#     
#     annual_index_sp_temp <- annual_index_sp[annual_index_sp$BOOTi==i,]
#     
#     lm_fit <- lm(LCI ~ M_YEAR, data = annual_index_sp_temp)
#     # extract and back-transform annual growth rate
#     growth_rate <- exp(coef(lm_fit)[2]*2.303)
#     perc_change <- 100*(exp(coef(lm_fit)[2]*2.303)-1)
#     total_perc_change <- 100*(exp(coef(lm_fit)[2]*2.303)^(max(annual_index_sp_temp$M_YEAR)-min(annual_index_sp_temp$M_YEAR))-1)
#     # build data frame
#     wcbs_growth_rates_t <- data.frame(growth_rate=growth_rate, perc_change=perc_change, total_perc_change=total_perc_change, BOOTi=i, sp=s, years=paste(start_year, "-", end_year, sep=""))
#     wcbs_growth_rates <- rbind(wcbs_growth_rates, wcbs_growth_rates_t)
#     }
#   }
# }
# 
# # calculate upper and lower CIs from the bootstrapped coefficients (excluding the 'true' one)
# bootstrap_CIs <- wcbs_growth_rates %>% group_by(sp, years) %>% filter(BOOTi>0) %>% 
#   summarise(lowerCI_gr=quantile(growth_rate, 0.025),upperCI_gr=quantile(growth_rate, 0.975),
#             lowerCI_perc=quantile(perc_change, 0.025), upperCI_perc=quantile(perc_change, 0.975),
#             lowerCI_perc_tot=quantile(total_perc_change, 0.025), upperCI_perc_tot=quantile(total_perc_change, 0.975))
# wcbs_growth_rates_true <- wcbs_growth_rates[wcbs_growth_rates$BOOTi==0,]
# wcbs_growth_rates2 <-  merge(wcbs_growth_rates_true, bootstrap_CIs, by=c("sp", "years"))
# wcbs_growth_rates2$SURVEY <- "WCBS"
# 
# 
# growth_rates <- rbind(tbms_growth_rates2,wcbs_growth_rates2)
# growth_rates <- merge(growth_rates, sp_names, by.x="sp", by.y="SPECIES")
# 
# # change from long to wide to plot
# growth_rates2 <- reshape(growth_rates, idvar = c("COMMON_NAME", "sp", "BOOTi", "years"), timevar = "SURVEY", direction = "wide")
# 
# # save file
# write.csv(growth_rates2, file="Output/GAM abundance analysis/Growth_rates_tBMS_WCBS_GAM_10year.csv", row.names=FALSE)
# growth_rates2 <- read.csv("Output/GAM abundance analysis/Growth_rates_tBMS_WCBS_GAM_10year.csv", header=TRUE)
# 
# # correlation for each 10-year window
# correlate <- growth_rates2 %>%
#   group_by(years) %>% 
#   summarise(r = cor(growth_rate.WCBS, growth_rate.tBMS))
# 
# min_max_values <- growth_rates2 %>% group_by(years) %>% summarise(max_tBMS=max(upperCI_gr.tBMS), min_WCBS=min(lowerCI_gr.WCBS))
# min_max_values <- as.data.frame(min_max_values)
# 
# dat_text <- data.frame(
#   label = c(paste("r =",round(correlate[1,2], digits=2)), paste("r =",round(correlate[2,2], digits=2)), 
#             paste("r =",round(correlate[3,2], digits=2)), paste("r =",round(correlate[4,2], digits=2)),
#             paste("r =",round(correlate[5,2], digits=2)), paste("r =",round(correlate[6,2], digits=2)), 
#             paste("r =",round(correlate[7,2], digits=2))),
#   years   = c("2009-2018", "2010-2019", "2011-2020", "2012-2021", "2013-2022", "2014-2023", "2015-2024"),
#   y = c(min_max_values[1,2]*0.98, min_max_values[2,2]*0.98, min_max_values[3,2]*0.98, min_max_values[4,2]*0.98,
#         min_max_values[5,2]*0.98, min_max_values[6,2]*0.97, min_max_values[7,2]*0.97),
#   x = c(min_max_values[1,3]*1.03, min_max_values[2,3]*1.03, min_max_values[3,3]*1.02, min_max_values[4,3]*1.02,
#         min_max_values[5,3]*1.03, min_max_values[6,3]*1.04, min_max_values[7,3]*1.03)
# )
# 
# growth_rate_p <- ggplot(growth_rates2, aes(growth_rate.WCBS, growth_rate.tBMS))+
#   geom_hline(yintercept = 1, color = "black", linetype="dotted")+
#   geom_vline(xintercept = 1, color = "black", linetype="dotted")+
#   geom_abline(linetype="dashed", color = "black")+
#   geom_point(size = 2)+
#   geom_errorbar(aes(ymin = lowerCI_gr.tBMS, ymax = upperCI_gr.tBMS))+
#   geom_errorbar(aes(xmin = lowerCI_gr.WCBS, xmax = upperCI_gr.WCBS))+
#   theme_bw()+
#   #ylim(-3,1.2)+
#   #xlim(-3,1.2)+
#   theme(text = element_text(size = 20)) +
#   xlab("WCBS abundance trend")+
#   ylab("tBMS abundance trend")+
#   facet_wrap(~years, scales="free")+
#   geom_text(
#     data    = dat_text, size=5,
#     mapping = aes(x = x, y = y, label = label))
# growth_rate_p
# ggsave(growth_rate_p, file="Output/Figures/GAM_results_gr_10years.png", height=30, width=30, units="cm")
