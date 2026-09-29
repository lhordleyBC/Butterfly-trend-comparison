# Compare net growth rates for wcbs and tBMS
library(data.table)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Read in net growth rate estimates for 2009-2023

tbms_netgr <- readRDS("Output/Growth rate analysis/tBMS/tbms_gra_net_2009_2024.rds")
wcbs_netgr <- readRDS("Output/Growth rate analysis/WCBS/WCBS_gra_net_2009_2024.rds")

setnames(tbms_netgr, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
         paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),"_tBMS"))
setnames(wcbs_netgr, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
         paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"), "_wcbs"))
both_ngr <- merge(tbms_netgr[, .(COMMON_NAME, ngr_tBMS, ngr_lower_tBMS, ngr_upper_tBMS, vcov_sum_tBMS, chat_tBMS)],
                  wcbs_netgr[, .(COMMON_NAME, ngr_wcbs, ngr_lower_wcbs, ngr_upper_wcbs, vcov_sum_wcbs, chat_wcbs)],
                  by = "COMMON_NAME")
write.csv(both_ngr, file="Output/Growth rate analysis/Species_growth_rates.csv", row.names=FALSE)

growth_rate_p <- ggplot(both_ngr, aes(ngr_wcbs, ngr_tBMS))+
  geom_hline(yintercept = 0, color = "black", linetype="dotted")+
  geom_vline(xintercept = 0, color = "black", linetype="dotted")+
  geom_abline(linetype="dashed", color = "black")+
  geom_point(size = 2)+
  geom_text_repel(aes(label = COMMON_NAME), color="grey", size = 5)+
  geom_errorbar(aes(ymin = ngr_lower_tBMS, ymax = ngr_upper_tBMS), width=0.1)+
  geom_errorbar(aes(xmin = ngr_lower_wcbs, xmax = ngr_upper_wcbs), width=0.1)+
  theme_bw()+
  ylim(-3,1.2)+
  xlim(-3,1.2)+
  theme(text = element_text(size = 20)) +
  xlab("Log WCBS growth rate 2009-2024")+
  ylab("Log tBMS growth rate 2009-2024")
growth_rate_p
ggsave(growth_rate_p, file="Output/Figures/Growth_rate_results_gr2.png", height=20, width=20, units="cm")

# Convert growth rates to pc?
both_ngr[, pc_tBMS := (exp(ngr_tBMS)-1)*100]
both_ngr[, pc_wcbs := (exp(ngr_wcbs)-1)*100]

both_ngr_w <- melt(both_ngr, id.vars = "COMMON_NAME", measure.vars = c("ngr_tBMS","ngr_wcbs"))
both_trends_w <- melt(both_ngr, id.vars = "COMMON_NAME", measure.vars = c("pc_tBMS","pc_wcbs"))

# Percentage change
ggplot(both_ngr, aes(pc_tBMS, pc_wcbs))+
  geom_hline(yintercept = 0, color = "grey")+
  geom_vline(xintercept = 0, color = "grey")+
  geom_abline(linetype="dashed", color = "grey")+
  geom_point(size = 2)+
  geom_text_repel(aes(label = COMMON_NAME), color="grey", size = 3)+
  theme_classic()+
  theme(text = element_text(size = 18)) +
  geom_smooth(method = "lm", se = FALSE)+
  xlab("tBMS percentage change 2009-2024")+
  ylab("wcbs percentage change 2009-2024")


# Percentage change as boxplot
ggplot(both_trends_w, aes(variable, value, fill = variable))+
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey")+
  theme_classic()+
  geom_boxplot(alpha=.5)+
  geom_line(aes(group=COMMON_NAME), position = position_dodge(0.2), col = "grey") +
  geom_point(aes(fill=variable, group=COMMON_NAME),
             size=2, shape=21, position = position_dodge(0.2)) +
  ylab("Percentage change 2009-2024")+
  scale_x_discrete(name="Scheme", labels=c("tBMS","wcbs"))+
  theme(legend.position = "none")

# Growth rates as boxplot
ggplot(both_ngr_w, aes(variable, value, fill = variable))+
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey")+
  theme_classic()+
  geom_boxplot(alpha=.5)+
  geom_line(aes(group=COMMON_NAME), position = position_dodge(0.2), col = "grey") +
  geom_point(aes(fill=variable, group=COMMON_NAME),
             size=2, shape=21, position = position_dodge(0.2)) +
  ylab("Log growth rate 2009-2024")+
  theme(legend.position = "none")


# Note these don't account for phylogeny
cor.test(both_ngr$ngr_tBMS, both_ngr$ngr_wcbs)

mean(both_ngr$ngr_tBMS)
mean(both_ngr$ngr_wcbs)

both_ngr[,ciwidth_tBMS := ngr_upper_tBMS - ngr_lower_tBMS]
both_ngr[,ciwidth_wcbs := ngr_upper_wcbs - ngr_lower_wcbs]

mean(both_ngr$ciwidth_tBMS)
mean(both_ngr$ciwidth_wcbs)
sd(both_ngr$ciwidth_tBMS)
sd(both_ngr$ciwidth_wcbs)

t.test(both_ngr$ciwidth_tBMS, both_ngr$ciwidth_wcbs, paired = TRUE)

mean(both_ngr$pc_tBMS)
mean(both_ngr$pc_wcbs)
std.error(both_ngr$pc_tBMS)
std.error(both_ngr$pc_wcbs)

t.test(both_ngr$pc_tBMS, both_ngr$pc_wcbs, paired = TRUE)

# Calculate differences between tBMS and WCBS and whether they are significant
both_ngr[, ngr_diff := ngr_tBMS - ngr_wcbs]
both_ngr[, ngr_diffSE := sqrt(vcov_sum_tBMS*chat_tBMS + vcov_sum_wcbs*chat_wcbs)]
both_ngr[, ngr_diff_lower := ngr_diff - 1.96*ngr_diffSE]
both_ngr[, ngr_diff_upper := ngr_diff + 1.96*ngr_diffSE]
both_ngr[, ngr_diff_sig := ifelse(ngr_diff_lower*ngr_diff_upper > 0, TRUE, FALSE)]
# save file
write.csv(both_ngr, file="Output/Growth rate analysis/Growth_rate_tBMS_wcbs.csv", row.names=FALSE)

both_ngr$ngr_diff_sig <- factor(both_ngr$ngr_diff_sig, levels = c("TRUE", "FALSE"))
gr_sig_diff <- ggplot(both_ngr, aes(ngr_wcbs, ngr_tBMS))+
  geom_hline(yintercept = 0, color = "black")+
  geom_vline(xintercept = 0, color = "black")+
  geom_abline(linetype="dashed", color = "black")+
  geom_errorbar(aes(ymin = ngr_lower_tBMS, ymax = ngr_upper_tBMS, color = ngr_diff_sig))+
  geom_errorbar(aes(xmin = ngr_lower_wcbs, xmax = ngr_upper_wcbs, color = ngr_diff_sig))+
  geom_text_repel(aes(label = COMMON_NAME), color="grey", size = 4)+
  geom_point(size = 2, aes(color = ngr_diff_sig))+
  theme_bw()+
  theme(text = element_text(size = 20)) +
  scale_color_manual(values=c("#339999","#FF6666"))+
  ylim(-3,1.2)+
  xlim(-3,1.2)+
  xlab("Log WCBS growth rate 2009-2024")+
  ylab("Log tBMS growth rate 2009-2024")+
  theme(legend.position = "bottom")+
  labs(color = "Significant difference in growth rate")
gr_sig_diff
ggsave(gr_sig_diff, file="Output/Figures/Growth_rate_sig_diff_gr.png", height=20, width=20, units="cm")

# also produce plot with difference in trend as x-axis in % 
# convert CIs to % too
both_ngr[, pc_lower_tBMS := (exp(ngr_lower_tBMS)-1)*100]
both_ngr[, pc_upper_tBMS := (exp(ngr_upper_tBMS)-1)*100]
both_ngr[, pc_lower_wcbs := (exp(ngr_lower_wcbs)-1)*100]
both_ngr[, pc_upper_wcbs := (exp(ngr_upper_wcbs)-1)*100]
# then calculate difference in the % CIs
both_ngr[, pc_diff_lower := pc_lower_tBMS - pc_lower_wcbs]
both_ngr[, pc_diff_upper := pc_upper_tBMS - pc_upper_wcbs]

both_ngr$COMMON_NAME = with(both_ngr, reorder(COMMON_NAME, pc_diff))
ggplot(both_ngr, aes(pc_diff, COMMON_NAME))+
  geom_point()+
  geom_errorbar(aes(xmin = pc_diff_lower, xmax = pc_diff_upper))+
  geom_vline(xintercept = 0, linetype="dashed", color = "grey")+
  theme_minimal()
# not convinced this is right

# same plot but for growth rate diff
both_ngr$significance <- ifelse(both_ngr$ngr_diff_lower<0 & both_ngr$ngr_diff_upper<0 | 
                                  both_ngr$ngr_diff_lower>0 & both_ngr$ngr_diff_upper>0, "yes", "no")

both_ngr$COMMON_NAME = with(both_ngr, reorder(COMMON_NAME, ngr_diff))
growth_rate_diff_p <- ggplot(both_ngr, aes(ngr_diff, COMMON_NAME))+
  geom_point(aes(colour=significance))+
  geom_errorbar(aes(xmin = ngr_diff_lower, xmax = ngr_diff_upper, colour=significance))+
  scale_colour_manual(values = c("darkgrey", "black")) +
  geom_vline(xintercept = 0, linetype="dashed", color = "black")+
  ylab("")+
  xlab("Growth rate difference \n(tBMS-WCBS)")+
  theme_bw()+
  theme(text = element_text(size = 20))+
  guides(colour="none")
growth_rate_diff_p
ggsave(growth_rate_diff_p, file="Output/Figures/Growth_rate_results_gr1.png", height=25, width=15, units="cm")

gr_results <- ggarrange(growth_rate_p, growth_rate_diff_p, labels=c("(a)", "(b)"), 
                          font.label=list(color="black",size=20))
gr_results
ggsave(gr_results, file="Output/Figures/Growthrate_results_gr.png", height=20, width=30, units="cm")

# Look at annual growth rate values
tbms_annual <- readRDS("Output/Growth rate analysis/tBMS/tBMS_gra_annual_2009_2024.rds")
wcbs_annual <- readRDS("Output/Growth rate analysis/wcbs/wcbs_gra_annual_2009_2024.rds")


tbms_annual[, agr_lower_tBMS := agr - qnorm(.975)*agr_se*sqrt(chat)]
tbms_annual[, agr_upper_tBMS := agr + qnorm(.975)*agr_se*sqrt(chat)]
tbms_annual[, agr_ciwidth_tBMS := agr_upper_tBMS - agr_lower_tBMS]

wcbs_annual[, agr_lower_wcbs := agr - qnorm(.975)*agr_se*sqrt(chat)]
wcbs_annual[, agr_upper_wcbs := agr + qnorm(.975)*agr_se*sqrt(chat)]
wcbs_annual[, agr_ciwidth_wcbs := agr_upper_wcbs - agr_lower_wcbs]

setnames(tbms_annual, "agr", "agr_tBMS")
setnames(tbms_annual, "agr_se", "agr_se_tBMS")
setnames(tbms_annual, "chat", "chat_tBMS")

setnames(wcbs_annual, "agr", "agr_wcbs")
setnames(wcbs_annual, "agr_se", "agr_se_wcbs")
setnames(wcbs_annual, "chat", "chat_wcbs")

both_annual <- merge(tbms_annual[, .(COMMON_NAME, YEAR, agr_tBMS, agr_lower_tBMS, agr_upper_tBMS, agr_ciwidth_tBMS, chat_tBMS)],
                     wcbs_annual[, .(COMMON_NAME, YEAR, agr_wcbs, agr_lower_wcbs, agr_upper_wcbs, agr_ciwidth_wcbs, chat_wcbs)],
                     by = c("COMMON_NAME", "YEAR"))

ggplot(both_annual, aes(agr_tBMS, agr_wcbs))+
  geom_hline(yintercept = 0, color = "grey")+
  geom_vline(xintercept = 0, color = "grey")+
  geom_abline(linetype="dashed", color = "grey")+
  geom_point(size = 2, aes(color = YEAR))+
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic()+
  theme(text = element_text(size = 18)) +
  xlab("Log tBMS growth rate")+
  ylab("Log wcbs growth rate")

years <- c(
  `2009` = "2009 - 2010",
  `2010` = "2010 - 2011",
  `2011` = "2011 - 2012",
  `2012` = "2012 - 2013",
  `2013` = "2013 - 2014",
  `2014` = "2014 - 2015",
  `2015` = "2015 - 2016",
  `2016` = "2016 - 2017",
  `2017` = "2017 - 2018",
  `2018` = "2018 - 2019",
  `2019` = "2019 - 2020",
  `2020` = "2020 - 2021",
  `2021` = "2021 - 2022",
  `2022` = "2022 - 2023",
  `2023` = "2023 - 2024")


# correlation for each year-year change 
correlate <- both_annual %>%
  group_by(YEAR) %>% 
  summarise(r = cor(agr_wcbs, agr_tBMS))

min_max_values <- both_annual %>% group_by(YEAR) %>% summarise(max_tBMS=max(agr_upper_tBMS), min_wcbs=min(agr_lower_wcbs))
min_max_values <- as.data.frame(min_max_values)

dat_text <- data.frame(
  label = c(paste("r =",round(correlate[1,2], digits=2)), paste("r =",round(correlate[2,2], digits=2)), 
            paste("r =",round(correlate[3,2], digits=2)), paste("r =",round(correlate[4,2], digits=2)),
            paste("r =",round(correlate[5,2], digits=2)), paste("r =",round(correlate[6,2], digits=2)), 
            paste("r =",round(correlate[7,2], digits=2)), paste("r =",round(correlate[8,2], digits=2)),
            paste("r =",round(correlate[9,2], digits=2)), paste("r =",round(correlate[10,2], digits=2)),
            paste("r =",round(correlate[11,2], digits=2)), paste("r =",round(correlate[12,2], digits=2)),
            paste("r =",round(correlate[13,2], digits=2)), paste("r =",round(correlate[14,2], digits=2)),
            paste("r =",round(correlate[15,2], digits=2))),
  YEAR   = c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023),
  y = c(min_max_values[1,2]*0.8, min_max_values[2,2]*0.8, min_max_values[3,2]*0.8, min_max_values[4,2]*0.8,
        min_max_values[5,2]*0.8, min_max_values[6,2]*0.8, min_max_values[7,2]*0.8, min_max_values[8,2]*0.8,
        min_max_values[9,2]*0.8, min_max_values[10,2]*0.8, min_max_values[11,2]*0.8, min_max_values[12,2]*0.8,
        min_max_values[13,2]*0.8, min_max_values[14,2]*0.8, min_max_values[15,2]-0.5),
  x = c(min_max_values[1,3]*0.7, min_max_values[2,3]*0.7, min_max_values[3,3]*0.7, min_max_values[4,3]*0.5,
        min_max_values[5,3]*0.7, min_max_values[6,3]*0.7, min_max_values[7,3]*0.7, min_max_values[8,3]*0.7,
        min_max_values[9,3]*0.7, min_max_values[10,3]*0.7, min_max_values[11,3]*0.7, min_max_values[12,3]*0.7,
        min_max_values[13,3]*0.8, min_max_values[14,3]*0.7, min_max_values[15,3]*0.85)
)

p <- ggplot(both_annual, aes(agr_wcbs, agr_tBMS))+
  facet_wrap(~YEAR ~ ., labeller = as_labeller(years), scales = "free", ncol=3)+
  geom_hline(yintercept = 0, color = "grey")+
  geom_vline(xintercept = 0, color = "grey")+
  geom_abline(linetype="dashed", color = "grey")+
  geom_point(size = 1.5)+
  theme_classic()+
  theme(text = element_text(size = 18)) +
  xlab("Log wcbs growth rate 2009-2024")+
  ylab("Log tBMS growth rate 2009-2024")+
  geom_errorbar(aes(ymin = agr_lower_tBMS, ymax = agr_upper_tBMS))+
  geom_errorbar(aes(xmin = agr_lower_wcbs, xmax = agr_upper_wcbs))+
  geom_text(
    data    = dat_text, size=4,
    mapping = aes(x = x, y = y, label = label))
p

ggsave(p, file="Output/Figures/Growthrate_annual_correlation_gr.png", height=12, width=9)

cor.test(both_annual$agr_tBMS, both_annual$agr_wcbs)
mean(both_annual$agr_tBMS)
mean(both_annual$agr_wcbs)
t.test(both_annual$agr_tBMS, both_annual$agr_wcbs, paired = TRUE)

mean(both_annual$agr_ciwidth_BBC)
mean(both_annual$agr_ciwidth_UKBMS)
t.test(both_annual$agr_ciwidth_BBC, both_annual$agr_ciwidth_UKBMS, paired = TRUE)

both_annual[, mean(agr_ciwidth_tBMS), by = YEAR]

ggplot(both_annual, aes(YEAR, agr_ciwidth_tBMS))+
  geom_point()
ggplot(both_annual, aes(YEAR, agr_ciwidth_wcbs))+
  geom_point()




#### 10-year rolling correlations and graphs

# read in tBMS data
tbms_netgr_2009_2018 <- readRDS("Output/Growth rate analysis/tBMS/tBMS_gra_net_2009_2018.rds")
setnames(tbms_netgr_2009_2018, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
         paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),"_tBMS"))
tbms_netgr_2010_2019 <- readRDS("Output/Growth rate analysis/tBMS/tBMS_gra_net_2010_2019.rds")
setnames(tbms_netgr_2010_2019, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
         paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),"_tBMS"))
tbms_netgr_2011_2020 <- readRDS("Output/Growth rate analysis/tBMS/tBMS_gra_net_2011_2020.rds")
setnames(tbms_netgr_2011_2020, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
         paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),"_tBMS"))
tbms_netgr_2012_2021 <- readRDS("Output/Growth rate analysis/tBMS/tBMS_gra_net_2012_2021.rds")
setnames(tbms_netgr_2012_2021, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
         paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),"_tBMS"))
tbms_netgr_2013_2022 <- readRDS("Output/Growth rate analysis/tBMS/tBMS_gra_net_2013_2022.rds")
setnames(tbms_netgr_2013_2022, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
         paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),"_tBMS"))
tbms_netgr_2014_2023 <- readRDS("Output/Growth rate analysis/tBMS/tBMS_gra_net_2014_2023.rds")
setnames(tbms_netgr_2014_2023, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
         paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),"_tBMS"))
tbms_netgr_2015_2024 <- readRDS("Output/Growth rate analysis/tBMS/tBMS_gra_net_2015_2024.rds")
setnames(tbms_netgr_2015_2024, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
         paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),"_tBMS"))

# read in WCBS data
wcbs_netgr_2009_2018 <- readRDS("Output/Growth rate analysis/WCBS/WCBS_gra_net_2009_2018.rds")
setnames(wcbs_netgr_2009_2018, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
         paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"), "_WCBS"))
wcbs_netgr_2010_2019 <- readRDS("Output/Growth rate analysis/WCBS/WCBS_gra_net_2010_2019.rds")
setnames(wcbs_netgr_2010_2019, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
         paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"), "_WCBS"))
wcbs_netgr_2011_2020 <- readRDS("Output/Growth rate analysis/WCBS/WCBS_gra_net_2011_2020.rds")
setnames(wcbs_netgr_2011_2020, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
         paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"), "_WCBS"))
wcbs_netgr_2012_2021 <- readRDS("Output/Growth rate analysis/WCBS/WCBS_gra_net_2012_2021.rds")
setnames(wcbs_netgr_2012_2021, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
         paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"), "_WCBS"))
wcbs_netgr_2013_2022 <- readRDS("Output/Growth rate analysis/WCBS/WCBS_gra_net_2013_2022.rds")
setnames(wcbs_netgr_2013_2022, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
         paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"), "_WCBS"))
wcbs_netgr_2014_2023 <- readRDS("Output/Growth rate analysis/WCBS/WCBS_gra_net_2014_2023.rds")
setnames(wcbs_netgr_2014_2023, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
         paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"), "_WCBS"))
wcbs_netgr_2015_2024 <- readRDS("Output/Growth rate analysis/WCBS/WCBS_gra_net_2015_2024.rds")
setnames(wcbs_netgr_2015_2024, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
         paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"), "_WCBS"))

## merge each 10-year growth rates and calculate correlation
# 2009-2018
both_ngr_2009_2018 <- merge(tbms_netgr_2009_2018[, .(COMMON_NAME, ngr_tBMS, ngr_lower_tBMS, ngr_upper_tBMS, vcov_sum_tBMS, chat_tBMS)],
                            wcbs_netgr_2009_2018[, .(COMMON_NAME, ngr_WCBS, ngr_lower_WCBS, ngr_upper_WCBS, vcov_sum_WCBS, chat_WCBS)],
                                        by = "COMMON_NAME")
cor.test(both_ngr_2009_2018$ngr_tBMS, both_ngr_2009_2018$ngr_WCBS) # 0.85 for 2009-2018
# 2010-2019
both_ngr_2010_2019 <- merge(tbms_netgr_2010_2019[, .(COMMON_NAME, ngr_tBMS, ngr_lower_tBMS, ngr_upper_tBMS, vcov_sum_tBMS, chat_tBMS)],
                            wcbs_netgr_2010_2019[, .(COMMON_NAME, ngr_WCBS, ngr_lower_WCBS, ngr_upper_WCBS, vcov_sum_WCBS, chat_WCBS)],
                            by = "COMMON_NAME")
cor.test(both_ngr_2010_2019$ngr_tBMS, both_ngr_2010_2019$ngr_WCBS) # 0.94 for 2010-2019
# 2011-2020
both_ngr_2011_2020 <- merge(tbms_netgr_2011_2020[, .(COMMON_NAME, ngr_tBMS, ngr_lower_tBMS, ngr_upper_tBMS, vcov_sum_tBMS, chat_tBMS)],
                            wcbs_netgr_2011_2020[, .(COMMON_NAME, ngr_WCBS, ngr_lower_WCBS, ngr_upper_WCBS, vcov_sum_WCBS, chat_WCBS)],
                            by = "COMMON_NAME")
cor.test(both_ngr_2011_2020$ngr_tBMS, both_ngr_2011_2020$ngr_WCBS) # 0.78 for 2011-2020
# 2012-2021
both_ngr_2012_2021 <- merge(tbms_netgr_2012_2021[, .(COMMON_NAME, ngr_tBMS, ngr_lower_tBMS, ngr_upper_tBMS, vcov_sum_tBMS, chat_tBMS)],
                            wcbs_netgr_2012_2021[, .(COMMON_NAME, ngr_WCBS, ngr_lower_WCBS, ngr_upper_WCBS, vcov_sum_WCBS, chat_WCBS)],
                            by = "COMMON_NAME")
cor.test(both_ngr_2012_2021$ngr_tBMS, both_ngr_2012_2021$ngr_WCBS) # 0.78 for 2012-2021
# 2013-2022
both_ngr_2013_2022 <- merge(tbms_netgr_2013_2022[, .(COMMON_NAME, ngr_tBMS, ngr_lower_tBMS, ngr_upper_tBMS, vcov_sum_tBMS, chat_tBMS)],
                            wcbs_netgr_2013_2022[, .(COMMON_NAME, ngr_WCBS, ngr_lower_WCBS, ngr_upper_WCBS, vcov_sum_WCBS, chat_WCBS)],
                            by = "COMMON_NAME")
cor.test(both_ngr_2013_2022$ngr_tBMS, both_ngr_2013_2022$ngr_WCBS) # 0.95 for 2013-2022
# 2014-2023
both_ngr_2014_2023 <- merge(tbms_netgr_2014_2023[, .(COMMON_NAME, ngr_tBMS, ngr_lower_tBMS, ngr_upper_tBMS, vcov_sum_tBMS, chat_tBMS)],
                            wcbs_netgr_2014_2023[, .(COMMON_NAME, ngr_WCBS, ngr_lower_WCBS, ngr_upper_WCBS, vcov_sum_WCBS, chat_WCBS)],
                            by = "COMMON_NAME")
cor.test(both_ngr_2014_2023$ngr_tBMS, both_ngr_2014_2023$ngr_WCBS) # 0.95 for 2014-2023
# 2015-2024
both_ngr_2015_2024 <- merge(tbms_netgr_2015_2024[, .(COMMON_NAME, ngr_tBMS, ngr_lower_tBMS, ngr_upper_tBMS, vcov_sum_tBMS, chat_tBMS)],
                            wcbs_netgr_2015_2024[, .(COMMON_NAME, ngr_WCBS, ngr_lower_WCBS, ngr_upper_WCBS, vcov_sum_WCBS, chat_WCBS)],
                            by = "COMMON_NAME")
cor.test(both_ngr_2015_2024$ngr_tBMS, both_ngr_2015_2024$ngr_WCBS) # 0.92 for 2009-2018

# merge all data together to plot
both_ngr_2009_2018$Years <- "2009-2018"
both_ngr_2010_2019$Years <- "2010-2019"
both_ngr_2011_2020$Years <- "2011-2020"
both_ngr_2012_2021$Years <- "2012-2021"
both_ngr_2013_2022$Years <- "2013-2022"
both_ngr_2014_2023$Years <- "2014-2023"
both_ngr_2015_2024$Years <- "2015-2024"

both_ngr_all <- rbind(both_ngr_2009_2018, both_ngr_2010_2019, both_ngr_2011_2020, both_ngr_2012_2021, both_ngr_2013_2022, both_ngr_2014_2023,
                      both_ngr_2015_2024)

# correlation for each 10-year change 
correlate <- both_ngr_all %>%
  group_by(Years) %>% 
  summarise(r = cor(ngr_WCBS, ngr_tBMS))

min_max_values <- both_ngr_all %>% group_by(Years) %>% summarise(max_tBMS=max(ngr_upper_tBMS), min_WCBS=min(ngr_lower_WCBS))
min_max_values <- as.data.frame(min_max_values)

dat_text <- data.frame(
  label = c(paste("r =",round(correlate[1,2], digits=2)), paste("r =",round(correlate[2,2], digits=2)), 
            paste("r =",round(correlate[3,2], digits=2)), paste("r =",round(correlate[4,2], digits=2)),
            paste("r =",round(correlate[5,2], digits=2)), paste("r =",round(correlate[6,2], digits=2)), 
            paste("r =",round(correlate[7,2], digits=2))),
  Years   = c("2009-2018", "2010-2019", "2011-2020", "2012-2021", "2013-2022", "2014-2023", "2015-2024"),
  y = c(min_max_values[1,2]*0.8, min_max_values[2,2]*0.8, min_max_values[3,2]*0.8, min_max_values[4,2]*0.8,
        min_max_values[5,2]*0.8, min_max_values[6,2]*0.8, min_max_values[7,2]*0.8),
  x = c(min_max_values[1,3]*0.7, min_max_values[2,3]*0.7, min_max_values[3,3]*0.7, min_max_values[4,3]*0.5,
        min_max_values[5,3]*0.8, min_max_values[6,3]*0.7, min_max_values[7,3]*0.8)
)

growth_rate_10years <- ggplot(both_ngr_all, aes(ngr_WCBS, ngr_tBMS))+
  geom_hline(yintercept = 0, color = "black", linetype="dotted")+
  geom_vline(xintercept = 0, color = "black", linetype="dotted")+
  geom_abline(linetype="dashed", color = "black")+
  geom_point(size = 2)+
  # geom_text_repel(aes(label = COMMON_NAME), color="grey", size = 3, max.overlaps = Inf)+
  geom_errorbar(aes(ymin = ngr_lower_tBMS, ymax = ngr_upper_tBMS))+
  geom_errorbar(aes(xmin = ngr_lower_WCBS, xmax = ngr_upper_WCBS))+
  theme_bw()+
  #ylim(-3,1.2)+
  #xlim(-3,1.2)+
  theme(text = element_text(size = 20)) +
  xlab("Log WCBS growth rate")+
  ylab("Log tBMS growth rate")+
  facet_wrap(~Years, scales="free")+
  geom_text(
    data    = dat_text, size=5,
    mapping = aes(x = x, y = y, label = label))
growth_rate_10years
ggsave(growth_rate_10years, file="Output/Figures/Growth_rate_results_10years.png", height=30, width=30, units="cm")

