### Traits analysis

library(dplyr)
library(ggplot2)
library(data.table)
library(ggeffects)
library(lme4)
library(rlang)
library(lmerTest)

# Load in files
# Growth rate data
growth_rates <- read.csv("Output/Growth rate analysis/Growth_rate_tBMS_WCBS.csv", header=TRUE)
# GAM data
growth_rates_gam <- read.csv("Output/GAM abundance analysis/Abundance_trends_tBMS_WCBS_GAM.csv", header=TRUE)
# Trait data
traits <- read.csv("Data/ecological_traits.csv", header=TRUE)
traits2 <- read.csv("Data/mobility_spec_traits.csv", header=TRUE)

# take only columns of interest (species and growth rate difference)
growth_rates <- growth_rates[,c("COMMON_NAME", "ngr_diff")]
growth_rates_gam <- growth_rates_gam[,c("COMMON_NAME", "growth_rate_diff")]
colnames(growth_rates)[2] <- "growth_rate_diff"
colnames(growth_rates_gam)[2] <- "GAM_diff"
growth_rates2 <- merge(growth_rates, growth_rates_gam, by="COMMON_NAME", all=TRUE)

# reduce columns in traits data
traits <- traits[,c("common_name", "family", "red_list_gb", "forewing_minimum", "forewing_maximum", "estimated_dry_mass",
                    "specificity", "number_hostplants", "nitrogen")]
growth_rates_traits <- merge(growth_rates2, traits, by.x="COMMON_NAME", by.y="common_name", all.x=TRUE)
growth_rates_traits <- merge(growth_rates_traits, traits2, by="COMMON_NAME", all.x=TRUE)
# no body mass data for butterflies
growth_rates_traits$specificity <- gsub("\\s*\\([^\\)]+\\)","",as.character(growth_rates_traits$specificity))

# Models 

## 1a. Red List: growth rate
model_redlist <- lmer(growth_rate_diff ~ red_list_gb + (1|family), data=growth_rates_traits)
# singular fit warning
summary(model_redlist) # significant - higher difference for threatened butterfiles (wall and small heath)
# higher growth rates on tBMS transects compared to WCBS squares for threatened butterflies
# but only 2 vs 20 so very uneven sample size 
ggpredict(model_redlist) %>% plot()
red_list <- data.frame(ggpredict(model_redlist))
red_list_p <- ggplot(red_list, aes(x=red_list_gb.x, y=red_list_gb.predicted))+
  geom_point()+
  geom_errorbar(aes(ymin=red_list_gb.conf.low, ymax=red_list_gb.conf.high), width=0)+
  geom_hline(yintercept=0, linetype="dashed")+
  labs(x="Red List status", y="Predicted trend difference \nbetween schemes")+
  theme_bw()
red_list_p

## 1b. Red list: GAM
model_redlist2 <- lmer(GAM_diff ~ red_list_gb + (1|family), data=growth_rates_traits)
# singular fit warning
summary(model_redlist2) # non-significant for GAM result


## 2a. Forewing length: growth rate
model_forewing <- lmer(growth_rate_diff ~ forewing_maximum + (1|family), data=growth_rates_traits)
# singular fit warning
summary(model_forewing) # non-sig for minimum and maximum

## 2b. Forewing length: GAM
model_forewing2 <- lmer(GAM_diff ~ forewing_maximum + (1|family), data=growth_rates_traits)
# singular fit warning
summary(model_forewing2) # non-sig for minimum and maximum
# p = 0.07, plot result to look at direction
pred <- ggpredict(model_forewing2, terms = "forewing_maximum")
forewing_plot <- ggplot(data = pred, aes(x = x, y = predicted)) +
  geom_line(linetype="dashed") +
  geom_point(data = growth_rates_traits, aes(x = forewing_maximum, y = GAM_diff),
             size = 3, shape = 21, fill = "skyblue", alpha = 0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  labs(x="Maximum forewing length", y="Abundance trend \ndifference")+
  theme_bw()
forewing_plot

## 3a. Hostplant specificity: growth rate
model_hostplant <- lmer(growth_rate_diff ~ specificity + (1|family), data=growth_rates_traits)
# singular fit warning
summary(model_hostplant) # non-significant for growth rate 

## 3b. Hostplant specificity: GAM
model_hostplant2 <- lmer(GAM_diff ~ specificity + (1|family), data=growth_rates_traits)
# singular fit warning
summary(model_hostplant2) # non-significant for GAM result


## 4a. Mobility score: growth rate
model_mobility <- lmer(growth_rate_diff ~ mobility_score + (1|family), data=growth_rates_traits)
# singular fit warning
summary(model_mobility) # non-significant for growth rate

## 4b. Mobility score: GAM
model_mobility2 <- lmer(GAM_diff ~ mobility_score + (1|family), data=growth_rates_traits)
# singular fit warning
summary(model_mobility2) # non-significant for GAM
ggpredict(model_mobility2) %>% plot() # higher mobility score = lower growth rate at tBMS sites compared to WCBS sites
# p = 0.08, plot result to look at direction
pred2 <- ggpredict(model_mobility2, terms = "mobility_score")
mobility_plot <- ggplot(data = pred2, aes(x = x, y = predicted)) +
  geom_line(linetype="dashed") +
  geom_point(data = growth_rates_traits, aes(x = mobility_score, y = GAM_diff),
             size = 3, shape = 21, fill = "skyblue", alpha = 0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  labs(x="Mobility score", y="Abundance trend \ndifference")+
  theme_bw()
mobility_plot

# put plots together for manuscript in case we put them in 
library(ggpubr)
trait_plots <- ggarrange(forewing_plot, mobility_plot)
ggsave(trait_plots, file="Output/Figures/Trait_plots.png", height=4, width=8)

## 5a. Biotype specialism: growth rate
model_spec <- lmer(growth_rate_diff ~ biotype_specialism + (1|family), data=growth_rates_traits)
# singular fit warning
summary(model_spec) # significant for growth rate - specialist species have lower growth rates on tBMS compared to WCBS
# BUT only one species is a habitat specialist (silver washed fritillary)
ggpredict(model_spec) %>% plot()

## 5b. Biotype specialism: GAM 
model_spec2 <- lmer(GAM_diff ~ biotype_specialism + (1|family), data=growth_rates_traits)
# singular fit warning
summary(model_spec2) # significant for GAM result too - same direction 
ggpredict(model_spec2) %>% plot()
gen_spec <- data.frame(ggpredict(model_spec2))
gen_spec_p <- ggplot(gen_spec, aes(x=biotype_specialism.x, y=biotype_specialism.predicted))+
  geom_point()+
  geom_errorbar(aes(ymin=biotype_specialism.conf.low, ymax=biotype_specialism.conf.high), width=0)+
  geom_hline(yintercept=0, linetype="dashed")+
  labs(x="Biotype specialism", y="Predicted trend difference \nbetween schemes")+
  theme_bw()
gen_spec_p

## 6a. Ellenberg nitrogen: growth rate



# Summary: 
  ## Red list: significant for GR only (more threatened butterflies have higher GR on tBMS sites compared to WCBS sites)
  ## Biotype specialism: significant for both (more specialist butterflies have higher GR on WCBS sites compared to tBMS sites) - odd result
  ## because only one species is specialist (silver-washed fritillary), and you would assume doing better on tBMS sites



