# Compare net growth rates for BBC and UKBMS
library(data.table)
library(ggplot2)
library(ggrepel)

# Read in net growth rate estimates for 2011-2021
# Exclude Large Skipper as removed from BBC 2018 onwards
bbc_netgr <- readRDS("Output/growthrateanalysis/BBC/bbc_gra_net_2011_2021.rds")[COMMON_NAME != "Large Skipper"]
bms_netgr <- readRDS("Output/growthrateanalysis/UKBMS/ukbms_gra_net_2011_2021_incWCBS.rds")

setnames(bbc_netgr, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
                    paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),"_BBC"))
setnames(bms_netgr, c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"),
                    paste0(c("ngr", "ngr_lower", "ngr_upper","vcov_sum","chat"), "_UKBMS"))

both_ngr <- merge(bbc_netgr[, .(COMMON_NAME, ngr_BBC, ngr_lower_BBC, ngr_upper_BBC, vcov_sum_BBC, chat_BBC)],
                  bms_netgr[, .(COMMON_NAME, ngr_UKBMS, ngr_lower_UKBMS, ngr_upper_UKBMS, vcov_sum_UKBMS, chat_UKBMS)],
                  by = "COMMON_NAME")

ggplot(both_ngr, aes(ngr_BBC, ngr_UKBMS))+
  geom_hline(yintercept = 0, color = "grey")+
  geom_vline(xintercept = 0, color = "grey")+
  geom_abline(linetype="dashed", color = "grey")+
  geom_point(size = 2)+
  geom_text_repel(aes(label = COMMON_NAME), color="grey", size = 3)+
  geom_errorbar(aes(ymin = ngr_lower_UKBMS, ymax = ngr_upper_UKBMS))+
  geom_errorbar(aes(xmin = ngr_lower_BBC, xmax = ngr_upper_BBC))+
  theme_classic()+
  theme(text = element_text(size = 18)) +
  geom_smooth(method = "lm", se = FALSE)+
  xlab("Log BBC growth rate 2011-2021")+
  ylab("Log UKBMS growth rate 2011-2021")

ggplot(both_ngr[!COMMON_NAME %in% c("Painted Lady","Marbled White")], aes(ngr_BBC, ngr_UKBMS))+
  geom_text_repel(aes(label = COMMON_NAME), color="grey", size = 3)+
  geom_hline(yintercept = 0, color = "grey")+
  geom_vline(xintercept = 0, color = "grey")+
  geom_abline(linetype="dashed", color = "grey")+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = ngr_lower_UKBMS, ymax = ngr_upper_UKBMS))+
  geom_errorbar(aes(xmin = ngr_lower_BBC, xmax = ngr_upper_BBC))+
  theme_classic()+
  theme(text = element_text(size = 18)) +
  geom_smooth(method = "lm", se = FALSE)+
  xlab("Log BBC growth rate 2011-2021")+
  ylab("Log UKBMS growth rate 2011-2021")



# Convert growth rates to pc?
both_ngr[, pc_BBC := (exp(ngr_BBC)-1)*100]
both_ngr[, pc_UKBMS := (exp(ngr_UKBMS)-1)*100]



both_ngr_w <- melt(both_ngr, id.vars = "COMMON_NAME", measure.vars = c("ngr_BBC","ngr_UKBMS"))
both_trends_w <- melt(both_ngr, id.vars = "COMMON_NAME", measure.vars = c("pc_BBC","pc_UKBMS"))



ggplot(both_ngr_w, aes(variable, value, fill = variable))+
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey")+
  theme_classic()+
  geom_boxplot(alpha=.5)+
  geom_line(aes(group=COMMON_NAME), position = position_dodge(0.2), col = "grey") +
  geom_point(aes(fill=variable, group=COMMON_NAME),
             size=2, shape=21, position = position_dodge(0.2)) +
  theme(legend.position = "none")


ggplot(both_trends_w[!(COMMON_NAME %in% c("Marbled White","Painted Lady"))], aes(variable, value, fill = variable))+
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey")+
  theme_classic()+
  geom_boxplot(alpha=.5)+
  geom_line(aes(group=COMMON_NAME), position = position_dodge(0.2), col = "grey") +
  geom_point(aes(fill=variable, group=COMMON_NAME),
             size=2, shape=21, position = position_dodge(0.2)) +
  theme(legend.position = "none")



ggplot(both_ngr, aes(pc_BBC, pc_UKBMS))+
  geom_hline(yintercept = 0, color = "grey")+
  geom_vline(xintercept = 0, color = "grey")+
  geom_abline(linetype="dashed", color = "grey")+
  geom_point(size = 2)+
  geom_text_repel(aes(label = COMMON_NAME), color="grey", size = 3)+
  theme_classic()+
  theme(text = element_text(size = 18)) +
  geom_smooth(method = "lm", se = FALSE)+
  xlab("Log BBC growth rate 2011-2021")+
  ylab("Log UKBMS growth rate 2011-2021")


# Note these don't account for phylogeny
cor.test(both_ngr$ngr_BBC, both_ngr$ngr_UKBMS)
cor.test(both_ngr[!COMMON_NAME %in% c("Painted Lady","Marbled White")]$ngr_BBC,
         both_ngr[!COMMON_NAME %in% c("Painted Lady","Marbled White")]$ngr_UKBMS)

t.test(both_ngr$ngr_BBC, both_ngr$ngr_UKBMS, paired = TRUE)

mean(both_ngr$ngr_BBC)
mean(both_ngr$ngr_UKBMS)

mean(both_ngr[!COMMON_NAME %in% c("Painted Lady","Marbled White")]$ngr_BBC)
mean(both_ngr[!COMMON_NAME %in% c("Painted Lady","Marbled White")]$ngr_UKBMS)
t.test(both_ngr[!COMMON_NAME %in% c("Painted Lady","Marbled White")]$ngr_BBC,
       both_ngr[!COMMON_NAME %in% c("Painted Lady","Marbled White")]$ngr_UKBMS, paired = TRUE)



both_ngr[,ciwidth_BBC := ngr_upper_BBC - ngr_lower_BBC]
both_ngr[,ciwidth_UKBMS := ngr_upper_UKBMS - ngr_lower_UKBMS]

mean(both_ngr$ciwidth_BBC)
mean(both_ngr$ciwidth_UKBMS)
t.test(both_ngr$ciwidth_BBC, both_ngr$ciwidth_UKBMS, paired = TRUE)


# Calculate differences between BBC and UKBMS and whether they are significant
both_ngr[, ngr_diff := ngr_BBC - ngr_UKBMS]
both_ngr[, ngr_diffSE := sqrt(vcov_sum_BBC*chat_BBC + vcov_sum_UKBMS*chat_UKBMS)]
both_ngr[, ngr_diff_lower := ngr_diff - 1.96*ngr_diffSE]
both_ngr[, ngr_diff_upper := ngr_diff + 1.96*ngr_diffSE]
both_ngr[, ngr_diff_sig := ifelse(ngr_diff_lower*ngr_diff_upper > 0, TRUE, FALSE)]


ggplot(both_ngr, aes(ngr_BBC, ngr_UKBMS))+
  geom_hline(yintercept = 0, color = "grey")+
  geom_vline(xintercept = 0, color = "grey")+
  geom_smooth(method = "lm", se = FALSE)+
  geom_abline(linetype="dashed", color = "grey")+
  geom_errorbar(aes(ymin = ngr_lower_UKBMS, ymax = ngr_upper_UKBMS))+
  geom_errorbar(aes(xmin = ngr_lower_BBC, xmax = ngr_upper_BBC))+
  geom_text_repel(aes(label = COMMON_NAME), color="grey", size = 3)+
  geom_point(size = 2, aes(color = ngr_diff_sig))+
  theme_classic()+
  theme(text = element_text(size = 18)) +
  xlab("Log BBC growth rate 2011-2019")+
  ylab("Log UKBMS growth rate 2011-2019")+
  theme(legend.position = "bottom")+
  labs(color = "Sig. difference in GR")




# Look at annual growth rate values
bbc_annual <- readRDS("Output/growthrateanalysis/BBC/bbc_gra_annual_2011_2021.rds")[COMMON_NAME != "Large Skipper"]
bms_annual <- readRDS("Output/growthrateanalysis/UKBMS/ukbms_gra_annual_2011_2021.rds")

bbc_annual[, agr_lower_BBC := agr - qnorm(.975)*agr_se*sqrt(chat)]
bbc_annual[, agr_upper_BBC := agr + qnorm(.975)*agr_se*sqrt(chat)]
bbc_annual[, agr_ciwidth_BBC := agr_upper_BBC - agr_lower_BBC]

bms_annual[, agr_lower_UKBMS := agr - qnorm(.975)*agr_se*sqrt(chat)]
bms_annual[, agr_upper_UKBMS := agr + qnorm(.975)*agr_se*sqrt(chat)]
bms_annual[, agr_ciwidth_UKBMS := agr_upper_UKBMS - agr_lower_UKBMS]

setnames(bbc_annual, "agr", "agr_BBC")
setnames(bms_annual, "agr", "agr_UKBMS")


both_annual <- merge(bbc_annual[, .(COMMON_NAME, YEAR, agr_BBC, agr_lower_BBC, agr_upper_BBC, agr_ciwidth_BBC)],
                  bms_annual[, .(COMMON_NAME, YEAR, agr_UKBMS, agr_lower_UKBMS, agr_upper_UKBMS, agr_ciwidth_UKBMS)],
                  by = c("COMMON_NAME", "YEAR"))

ggplot(both_annual, aes(agr_BBC, agr_UKBMS))+
  geom_hline(yintercept = 0, color = "grey")+
  geom_vline(xintercept = 0, color = "grey")+
  geom_abline(linetype="dashed", color = "grey")+
  geom_point(size = 2, aes(color = YEAR))+
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic()+
  theme(text = element_text(size = 18)) +
  xlab("Log BBC growth rate")+
  ylab("Log UKBMS growth rate")

ggplot(both_annual, aes(agr_BBC, agr_UKBMS))+
  facet_wrap(~YEAR)+
  geom_hline(yintercept = 0, color = "grey")+
  geom_vline(xintercept = 0, color = "grey")+
  geom_abline(linetype="dashed", color = "grey")+
  geom_point(size = 1.5)+
  theme_classic()+
  theme(text = element_text(size = 14)) +
  geom_smooth(method = "lm", se = FALSE)+
  xlab("Log BBC growth rate")+
  ylab("Log UKBMS growth rate")+
  geom_errorbar(aes(ymin = agr_lower_UKBMS, ymax = agr_upper_UKBMS))+
  geom_errorbar(aes(xmin = agr_lower_BBC, xmax = agr_upper_BBC))


ggplot(both_annual[!(COMMON_NAME %in% "Painted Lady")], aes(agr_BBC, agr_UKBMS))+
  facet_wrap(~YEAR, scales = "free")+
  geom_hline(yintercept = 0, color = "grey")+
  geom_vline(xintercept = 0, color = "grey")+
  geom_abline(linetype="dashed", color = "grey")+
  geom_point(size = 1.5)+
  theme_classic()+
  theme(text = element_text(size = 14)) +
  geom_smooth(method = "lm", se = FALSE)+
  xlab("Log BBC growth rate")+
  ylab("Log UKBMS growth rate")+
  geom_errorbar(aes(ymin = agr_lower_UKBMS, ymax = agr_upper_UKBMS))+
  geom_errorbar(aes(xmin = agr_lower_BBC, xmax = agr_upper_BBC))


cor.test(both_annual$agr_BBC, both_annual$agr_UKBMS)
mean(both_annual$agr_BBC)
mean(both_annual$agr_UKBMS)
t.test(both_annual$agr_BBC, both_annual$agr_UKBMS, paired = TRUE)


mean(both_annual$agr_ciwidth_BBC)
mean(both_annual$agr_ciwidth_UKBMS)
t.test(both_annual$agr_ciwidth_BBC, both_annual$agr_ciwidth_UKBMS, paired = TRUE)

both_annual[, mean(agr_ciwidth_BBC), by = YEAR]

ggplot(both_annual, aes(YEAR, agr_ciwidth_BBC))+
  geom_point()
ggplot(both_annual, aes(YEAR, agr_ciwidth_UKBMS))+
  geom_point()





# Compare BBC growth rates and UKBMS GAI trends
# Exclude Large Skipper as removed from BBC 2018 onwards
bms_gai <- fread(paste0("Output/GAI_UKBMS/UKBMS_GAI_trends_2011_2021_JULAUG.csv"))
setnames(bms_gai, "pc","pc_gai")

bbc_vs_gai <-  merge(bbc_netgr, bms_gai, by = "COMMON_NAME")
bbc_vs_gai[, pc_bbc := (exp(ngr_BBC)-1)*100]


ggplot(bbc_vs_gai[COMMON_NAME != "Painted Lady"], aes(pc_gai, pc_bbc))+
  geom_hline(yintercept = 0, color = "grey")+
  geom_vline(xintercept = 0, color = "grey")+
  geom_abline(linetype="dashed", color = "grey")+
  geom_point(size = 2)+
  geom_text_repel(aes(label = COMMON_NAME), color="grey", size = 3)+
  theme_classic()+
  theme(text = element_text(size = 18)) +
  xlab("% change 2011-2021 (UKBMS GAI)")+
  ylab("% change 2011-2021 (BBC)")



bms_vs_gai <-  merge(bms_netgr[,.(COMMON_NAME, ngr_UKBMS)], bms_gai, by = "COMMON_NAME")
bms_vs_gai[, pc_bms := (exp(ngr_UKBMS)-1)*100]

ggplot(bms_vs_gai[!(COMMON_NAME %in% c("Painted Lady", "Marbled White"))], aes(pc_bms, pc_gai))+
  geom_hline(yintercept = 0, color = "grey")+
  geom_vline(xintercept = 0, color = "grey")+
  geom_abline(linetype="dashed", color = "grey")+
  geom_point(size = 2)+
  geom_text_repel(aes(label = COMMON_NAME), color="grey", size = 3)+
  theme_classic()+
  theme(text = element_text(size = 18)) +
  #geom_smooth(method = "lm", se = FALSE)+
  xlab("% change 2011-2019 (UKBMS gr)")+
  ylab("% change 2011-2019 (UKBMS GAI)")


bms_vs_gai_w <- melt(bms_vs_gai, id.vars = "COMMON_NAME", measure.vars = c("pc_bms","pc_gai"))

ggplot(bms_vs_gai_w[!(COMMON_NAME %in% c("Painted Lady", "Marbled White"))],
       aes(variable, value, fill = variable))+
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey")+
  theme_classic()+
  geom_boxplot(alpha=.5)+
  geom_line(aes(group=COMMON_NAME), position = position_dodge(0.2), col = "grey") +
  geom_point(aes(fill=variable, group=COMMON_NAME),
             size=2, shape=21, position = position_dodge(0.2)) +
  theme(legend.position = "none")


