# Resource Quality Experiment Analysis

# Code Written by: Michelle L Fearon
# Last updated: 4/24/2024
# Last update by: Michelle L Fearon


# Terminology Used

      # Bacterium = Pasteruia ramosa
      # Fungus = Metschnikowia bicuspidata







# Libraries
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(emmeans)
library(ggeffects)
library(epiR)
library(car)
library(lubridate)
library(survival)
library(rms)
library(rsq)
library(survminer)
library(ggplot2)
library(viridis)
library(EnvStats)
library(olsrr)
library(DHARMa)
library(RColorBrewer)
library(here)

# set the path to the script relative to the project root directory
here::i_am("scripts/ResourceQuality_Analysis_cleanMLF.R")
#setwd("/Users/duffymeg/Dropbox (University of Michigan)/Manuscripts/FearonResourceQuality/2024/resource-quality-expt-main")
setwd("C:/Users/mlfearon/OneDrive - Umich/PROJECTS/MHMP Daphnia Duffy/Resource Quality/Data and Code/resource-quality-expt")

# load Pasteuria data
p_data <- read.csv(here("data/ResourceQuality_Pasteuria_Full.csv"), stringsAsFactors = F)
p_longdata <- read.csv(here("data/ResourceQuality_Pasteuria_Full_ByExptWeek.csv"), stringsAsFactors = F)
p_lifetable <- read.csv(here("data/ResourceQuality_Pasteuria_LifeTableResults.csv"), stringsAsFactors = F)
p_lifetable_clone_inf <- read.csv(here("data/ResourceQuality_Pasteuria_LifeTable_clone_x_infstatus_Results.csv"), stringsAsFactors = F)

# load Metsch data
m_data <- read.csv(here("data/ResourceQuality_Metsch_Full.csv"), stringsAsFactors = F)
m_longdata <- read.csv(here("data/ResourceQuality_Metsch_Full_ByExptWeek.csv"), stringsAsFactors = F)
m_lifetable <- read.csv(here("data/ResourceQuality_Metsch_LifeTableResults.csv"), stringsAsFactors = F)
m_lifetable_clone_inf <- read.csv(here("data/ResourceQuality_Metsch_LifeTable_clone_x_infstatus_Results.csv"), stringsAsFactors = F)


# Reorder factors
p_data$Diet <- factor(p_data$Diet, levels = c("S", "SM", "M", "M+"))
p_data$Infection <- factor(p_data$Infection, levels = c("Uninfected", "Exposed", "Infected"))
p_data$Parasites <- factor(p_data$Parasites, levels = c("Uninfected", "Pasteuria"))
p_longdata$Diet <- factor(p_longdata$Diet, levels = c("S", "SM", "M", "M+"))
p_longdata$Infection <- factor(p_longdata$Infection, levels = c("Uninfected", "Exposed", "Infected"))
p_lifetable$Diet <- factor(p_lifetable$Diet, levels = c("S", "SM", "M", "M+"))
p_lifetable$Infection <- factor(p_lifetable$Infection, levels = c("Uninfected", "Exposed", "Infected"))
p_lifetable_clone_inf$Infection <- factor(p_lifetable_clone_inf$Infection, levels = c("Uninfected", "Exposed", "Infected"))

m_data$Diet <- factor(m_data$Diet, levels = c("S", "SM", "M", "M+"))
m_data$Infection <- factor(m_data$Infection, levels = c("Uninfected", "Exposed", "Infected"))
m_data$Parasites <- factor(m_data$Parasites, levels = c("Uninfected", "Metsch"))
m_longdata$Diet <- factor(m_longdata$Diet, levels = c("S", "SM", "M", "M+"))
m_longdata$Infection <- factor(m_longdata$Infection, levels = c("Uninfected", "Exposed", "Infected"))
m_lifetable$Diet <- factor(m_lifetable$Diet, levels = c("S", "SM", "M", "M+"))
m_lifetable$Infection <- factor(m_lifetable$Infection, levels = c("Uninfected", "Exposed", "Infected"))
m_lifetable_clone_inf$Infection <- factor(m_lifetable_clone_inf$Infection, levels = c("Uninfected", "Exposed", "Infected"))


# Overdispersion parameter estimation function from Ben Bolker: http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#testing-for-overdispersioncomputing-overdispersion-factor
# Comment additions from http://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}



# Parasite infectivity and host susceptibility ----------------------------------

### Comparison of Metsch infection prevalence based on diet, and host clone.

# remove animals that died during exposure and remove uninfected treatments
m_data_prev <- m_data %>% filter(!is.na(InfectionStatus) & Parasites == "Metsch") 

m_data_prev$Block <- as.factor(m_data_prev$Block)

mprevmod <- glmer(InfectionStatus ~ Diet * Clone + (1|Block), family = binomial, data = m_data_prev)
# Note: left off the block random effect here since it seemed to cause a problem with the model... 
summary(mprevmod)
Anova(mprevmod)
overdisp_fun(mprevmod)
AIC(mprevmod)
testDispersion(mprevmod)
testZeroInflation(mprevmod)
mprev_simResid <- simulateResiduals(fittedModel = mprevmod)
plot(mprev_simResid)  # model looks good

mprevmod_contrasts <- emmeans(mprevmod, specs = pairwise ~ Diet | Clone, type = "response")
mprevmod_contrasts

mprevmod_contrasts2 <- emmeans(mprevmod, specs = pairwise ~ Clone | Diet, type = "response")
mprevmod_contrasts2



## Fungus infectivity results in manuscript
# test a model with only S and SM diets for Mid37 to remove treatments with no infection (except Std M+ did have a little infection)
m_data_prev_std <- filter(m_data_prev, Clone == "Standard")
m_data_prev_midS_SM <- filter(m_data_prev, Clone == "Mid37", Diet != "M", Diet != "M+")
m_data_prev2 <- rbind(m_data_prev_midS_SM, m_data_prev_std)

mprevmod2 <- glmer(InfectionStatus ~ Diet + Clone + (1|Block), family = binomial, data = m_data_prev2)
summary(mprevmod2)  

# Table 1: Fungus prevalence
Anova(mprevmod2) # reported model results

# block random effect does work in this model, warnings because there are NAs for some 
# groups when the interaction is present(the M and M+ treatments for Mid37s that I took out). 
# Since the interaction is not sig, I removed it and this model looks much better with 
# lower Std error values. --> USE THIS MODEL in manuscript

# checking model assumptions and fit
overdisp_fun(mprevmod2)
AIC(mprevmod2)
testDispersion(mprevmod2)
testZeroInflation(mprevmod2)
mprev2_simResid <- simulateResiduals(fittedModel = mprevmod2)
plot(mprev2_simResid)  # model looks good

mprevmod_contrasts3 <- emmeans(mprevmod2, specs = pairwise ~ Diet | Clone, type = "response")
mprevmod_contrasts3



# Marginal effects from best model with Abundance on the x axis
metsch_prev <- ggpredict(mprevmod2, c("Clone", "Diet"))


# figures made with marginal effects (not used in manuscript)
m_prev_fig1 <- plot(metsch_prev, dodge = 0.5, show.title = F) +
  scale_color_manual(values = diet_colors) +
  #ggtitle("Metsch prevalence by diet x host clone") +
  labs(x = "Host Clone", y = "Fungus Prevalence (%)") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=9, color="black"))
m_prev_fig1
ggsave("figures/MetschPrev_DietxClone.tiff", plot = m_prev_fig1, dpi = 600, width = 4, height = 3, units = "in", compression="lzw")



# calculate prevalence for each diet and clone
m_data_prev$Rep <- as.character(m_data_prev$Rep)

metsch_prev_sum <- m_data_prev %>%
  group_by(Diet, Clone) %>%
  summarise(Infected = sum(InfectionStatus), Tested = length(Rep),
            N  = length(Clearance_rel.Week1), 
            Mean_Clearance_rel = mean(Clearance_rel.Week1, na.rm = T),
            clearance.rel.sd   = sd(Clearance_rel.Week1, na.rm = T),
            clearance.rel.se   = clearance.rel.sd / sqrt(N),
            Mean_Clearance = mean(Clearance.Week1, na.rm = T),
            clearance.sd   = sd(Clearance.Week1, na.rm = T),
            clearance.se   = clearance.sd / sqrt(N),
            .groups = "rowwise") %>%
  mutate(Prevalence = unlist(epi.prev(Infected, Tested, method = "blaker", conf.level = 0.95, sp = 1, se = 0.95))[1],
         Prev_Lower = unlist(epi.prev(Infected, Tested, method = "blaker", conf.level = 0.95, sp = 1, se = 0.95))[2],
         Prev_Upper = unlist(epi.prev(Infected, Tested, method = "blaker", conf.level = 0.95, sp = 1, se = 0.95))[3])


## Figure 1A: Fungus prevalence for each clone and diet treatment (used in manuscript)
# figures  with calculated prevalence and confidence intervals
m_prev_fig <- ggplot(metsch_prev_sum, aes(x = Diet, y= Prevalence)) +
  geom_point(aes(color = Diet), size = 2, show.legend = F) +
  geom_errorbar(aes(ymin = Prev_Lower, ymax = Prev_Upper, color = Diet), width = 0.2, show.legend = F) + 
  scale_color_manual(values = diet_colors) +
  labs(x = "Diet Treatment", y = "Fungus Prevalence (%)") +
  ggtitle("Fungus Experiment") +
  ylim(0, 100) +
  facet_wrap(~Clone) +
  theme_classic() +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
m_prev_fig
ggsave("figures/MetschPrevCalcualtions_DietxClone.tiff", plot = m_prev_fig, dpi = 600, width = 4, height = 3, units = "in", compression="lzw")



# Results presented in Appendix 
# S1.6 Results for feeding rate and infection prevalence
m_prev_feed_model <- glmer(InfectionStatus ~ Clearance.Week1 * Clone + (1|Diet) + (1|Block), family = "binomial", glmerControl(optimizer = "bobyqa"), data = m_data_prev)
summary(m_prev_feed_model) 
Anova(m_prev_feed_model)
overdisp_fun(m_prev_feed_model)

mprevfeed_contrasts <- emtrends(m_prev_feed_model, specs = ~ Clone, var = "Clearance.Week1", type = "response")
mprevfeed_contrasts


## Figure S3A in Appendix
# Fungus prevalence vs Feeding rate per hr
m_prev_feed_fig <- ggplot(metsch_prev_sum, aes(x = Mean_Clearance, y= Prevalence)) +
  geom_point(aes(color = Diet, shape = Clone), size = 2) +
  geom_errorbar(aes(ymin = Prev_Lower, ymax = Prev_Upper, color = Diet), width = 0.01) + 
  geom_errorbarh(aes(xmin = Mean_Clearance-clearance.se, xmax = Mean_Clearance+clearance.se, color = Diet), height = 3) + 
  scale_color_manual(values = diet_colors) +
  ggtitle("Fungus Experiment") +
  labs(x = "Mean Feeding Rate (mL per hr)", y = "Fungus Prevalence (%)") +
  ylim(-2, 100) +
  facet_wrap(~Clone) +
  theme_classic() +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
m_prev_feed_fig
ggsave(here("figures/MetschPrevCalcualtions_FeedingRate.tiff"), plot = m_prev_feed_fig, dpi = 600, width = 4.5, height = 3, units = "in", compression="lzw")


#### Figure S3 in Appendix
FigS3 <- ggarrange(m_prev_feed_fig, p_prev_feed_fig, nrow = 1, ncol = 2, labels = "auto", common.legend = T, legend = "right")
ggsave(here("figures/manuscript/FigS3_Prevalence_FeedingRate.tiff"), plot = FigS3, dpi = 600, width = 8, height = 3.5, units = "in", compression="lzw")


# Relative feeding rate (not in manuscript)
m_prev_feed_rel_fig <- ggplot(metsch_prev_sum, aes(x = Mean_Clearance_rel, y= Prevalence)) +
  geom_point(aes(color = Diet, shape = Clone), size = 2) +
  geom_errorbar(aes(ymin = Prev_Lower, ymax = Prev_Upper, color = Diet), width = 0.01) + 
  geom_errorbarh(aes(xmin = Mean_Clearance_rel-clearance.rel.se, xmax = Mean_Clearance_rel+clearance.rel.se, color = Diet), width = 0.01) + 
  scale_color_manual(values = diet_colors) +
  labs(x = "Mean Relative Feeding Rate (mL per hr per mm)", y = "Fungus Prevalence (%)") +
  ylim(0, 100) +
  facet_wrap(~Clone) +
  theme_classic() +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
m_prev_feed_rel_fig
ggsave(here("figures/MetschPrevCalcualtions_RelFeedingRate3.tiff"), plot = m_prev_feed_rel_fig, dpi = 600, width = 4.5, height = 3, units = "in", compression="lzw")





### Comparison of Pasteuria infection prevalence based on diet and host clone.

# remove animals that died during exposure and remove uninfected treatments
p_data_prev <- p_data %>% filter(!is.na(InfectionStatus) & Parasites == "Pasteuria") 

pprevmod <- glm(InfectionStatus ~ Diet * Clone, family = binomial, data = p_data_prev)
# Note: left off the block random effect here since it seemed to cause a problem with the model.
summary(pprevmod)

# Table 1: Bacterium prevalence
Anova(pprevmod) # reported model results

# checking model assumptions and fit
overdisp_fun(pprevmod)
testDispersion(pprevmod)
testZeroInflation(pprevmod)
pprev_simResid <- simulateResiduals(fittedModel = pprevmod)
plot(pprev_simResid)  # model looks good

pprevmod_contrasts <- emmeans(pprevmod, specs = pairwise ~ Diet | Clone, type = "response")
pprevmod_contrasts



# Marginal effects from best model with Abundance on the x axis
past_prev <- ggpredict(pprevmod, c("Clone", "Diet"))

# set the color scheme for the Diet treatments
diet_colors <- c("#ADDD8E", "#41AB5D", "#006837", "#1D91C0")


# figures made with marginal effects (not used in manuscript)
p_prev_fig1 <- plot(past_prev, dodge = 0.5, show.title = F) +
  scale_color_manual(values = diet_colors) +
  #ggtitle("Pasteuria prevalence by diet x host clone") +
  labs(x = "Host Clone", y = "Bacterium Prevalence (%)") +
  theme_classic() +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=9, color="black"))
p_prev_fig1
ggsave(here("figures/PasteuriaPrev_DietxClone.tiff"), plot = p_prev_fig1, dpi = 600, width = 4, height = 3, units = "in", compression="lzw")




# calculate prevalence for each diet and clone
p_data_prev$Rep <- as.character(p_data_prev$Rep)

past_prev_sum <- p_data_prev %>%
  group_by(Diet, Clone) %>%
  summarise(Infected = sum(InfectionStatus), Tested = length(Rep), 
            N  = length(Clearance_rel.Week1), 
            Mean_Clearance_rel = mean(Clearance_rel.Week1, na.rm = T),
            clearance.rel.sd   = sd(Clearance_rel.Week1, na.rm = T),
            clearance.rel.se   = clearance.rel.sd / sqrt(N), 
            Mean_Clearance = mean(Clearance.Week1, na.rm = T),
            clearance.sd   = sd(Clearance.Week1, na.rm = T),
            clearance.se   = clearance.sd / sqrt(N),
            .groups = "rowwise") %>%
  mutate(Prevalence = unlist(epi.prev(Infected, Tested, method = "blaker", conf.level = 0.95, sp = 1, se = 0.95))[1],
         Prev_Lower = unlist(epi.prev(Infected, Tested, method = "blaker", conf.level = 0.95, sp = 1, se = 0.95))[2],
         Prev_Upper = unlist(epi.prev(Infected, Tested, method = "blaker", conf.level = 0.95, sp = 1, se = 0.95))[3])


## Figure 1B: Bacterium prevalence for each clone and diet treatment (used in manuscript)
# figures  with calculated prevalence and confidence intervals
p_prev_fig <- ggplot(past_prev_sum, aes(x = Diet, y= Prevalence)) +
  geom_point(aes(color = Diet), size = 2, show.legend = F) +
  geom_errorbar(aes(ymin = Prev_Lower, ymax = Prev_Upper, color = Diet), width = 0.2, show.legend = F) + 
  scale_color_manual(values = diet_colors) +
  labs(x = "Diet Treatment", y = "Bacterium Prevalence (%)") +
  ggtitle("Bacterium Experiment") +
  ylim(0, 100) +
  facet_wrap(~Clone) +
  theme_classic() +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
p_prev_fig
ggsave(here("figures/PasteuriaPrevCalcualtions_DietxClone.tiff"), plot = p_prev_fig, dpi = 600, width = 4, height = 3, units = "in", compression="lzw")


# Results presented in Appendix 
#S1.6 Results for feeding rate and infection prevalence
p_prev_feed_model <- glmer(InfectionStatus ~ Clearance.Week1 * Clone + (1|Diet) + (1|Block), family = "binomial", data = p_data_prev)
summary(p_prev_feed_model)
Anova(p_prev_feed_model)


## Figure S3B in Appendix
# Bacterium prevalence vs Feeding rate per hr
p_prev_feed_fig <- ggplot(past_prev_sum, aes(x = Mean_Clearance, y= Prevalence)) +
  geom_point(aes(color = Diet, shape = Clone), size = 2) +
  geom_errorbar(aes(ymin = Prev_Lower, ymax = Prev_Upper, color = Diet), width = 0.015) + 
  geom_errorbarh(aes(xmin = Mean_Clearance-clearance.se, xmax = Mean_Clearance+clearance.se, color = Diet), height = 3) + 
  scale_color_manual(values = diet_colors) +
  ggtitle("Bacterium Experiment") +
  labs(x = "Mean Feeding Rate (mL per hr)", y = "Bacterium Prevalence (%)") +
  ylim(-20, 100) +
  facet_wrap(~Clone) +
  theme_classic() +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=10, color="black"))
p_prev_feed_fig
ggsave(here("figures/PasteuriaPrevCalcualtions_FeedingRate.tiff"), plot = p_prev_feed_fig, dpi = 600, width = 4, height = 3, units = "in", compression="lzw")

# Relative feeding rate (not shown in manuscript)
p_prev_feed_rel_fig <- ggplot(past_prev_sum, aes(x = Mean_Clearance_rel, y= Prevalence)) +
  geom_point(aes(color = Diet, shape = Clone), size = 2) +
  geom_errorbar(aes(ymin = Prev_Lower, ymax = Prev_Upper, color = Diet), width = 0.02) + 
  geom_errorbarh(aes(xmin = Mean_Clearance_rel-clearance.rel.se, xmax = Mean_Clearance_rel+clearance.rel.se, color = Diet), height = 3) + 
  scale_color_manual(values = diet_colors) +
  labs(x = "Mean Relative Feeding Rate (mL per hr per mm)", y = "Bacterium Prevalence (%)") +
  ylim(-20, 100) +
  facet_wrap(~Clone) +
  theme_classic() +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=10, color="black"))
p_prev_feed_rel_fig
ggsave(here("figures/PasteuriaPrevCalcualtions_RelFeedingRate2.tiff"), plot = p_prev_feed_rel_fig, dpi = 600, width = 4, height = 3, units = "in", compression="lzw")






## Parasite exposure and feeding rate --------------------------------------------


### Fungus feeding rate analyses
# Feeding rate during the inoculation period (week 1 of the study)

m_data_feeding <- m_data %>% filter(!is.na(Clearance.Week1))
hist(m_data_feeding$Clearance.Week1) # wow normal data!! 
hist(m_data_feeding$Clearance_rel.Week1) # wow normal data!! 


mfeedmod <- lmer(Clearance.Week1 ~ Diet + Infection + Clone + Diet:Infection + Diet:Clone + Infection:Clone +(1|Block), data = m_data_feeding)
# removed the three-way interaction, not enough data for the model to run the 3-way interaction
summary(mfeedmod)

# Table 1: Feeding rate/Parasite exposure rate
Anova(mfeedmod)
AIC(mfeedmod)
plot(mfeedmod)
qqnorm(resid(mfeedmod))
qqline(resid(mfeedmod))
shapiro.test(resid(mfeedmod))

# Appendix Table S3: Fungus experiment feeding rate/exposure rate
# contrasts of diet effects by each infection status and clones
mfeedmod_contrasts <- emmeans(mfeedmod, specs = pairwise ~ Diet | Infection + Clone, type = "response")
mfeedmod_contrasts 

# Appendix Table S4: Fungus experiment feeding rate/exposure rate
# contrasts of infection status effects by each clone (averaged over all diets)
mfeedmod_contrasts2 <- emmeans(mfeedmod, specs = pairwise ~ Infection | Clone, type = "response")
mfeedmod_contrasts2  

# assess the nearly sig interaction between infection and clone
mfeedmod_contrasts3 <- emmeans(mfeedmod, specs = pairwise ~ Clone | Infection, type = "response")
mfeedmod_contrasts3


# Relative Feeding rate (bodysize controlled) - analysis described in Appendix Section S1.7
mfeedmod2 <- lmer(Clearance_rel.Week1 ~ Diet + Infection + Clone + Diet:Infection + Diet:Clone + Infection:Clone +(1|Block), data = m_data_feeding)
# removed the three-way interaction, not enough data for the model to run the 3-way interaction
summary(mfeedmod2)

# Appendix Table S2: Fungus Experiment Relative Feeding Rate/Exposure Rate
Anova(mfeedmod2) 
AIC(mfeedmod2)
plot(mfeedmod2)
qqnorm(resid(mfeedmod2))
qqline(resid(mfeedmod2))
shapiro.test(resid(mfeedmod2))

# contrasts of diet effects by each infection status and clones
mfeedmod_contrasts4 <- emmeans(mfeedmod2, specs = pairwise ~ Diet | Infection + Clone, type = "response")
mfeedmod_contrasts4   # Mid37 there is a dif between S - M for uninfected and SM - M+ for exposed, and no differences for infected; and for Std uninfected there are dif between S - M, S - M+, and SM - M, for exposed and infected dif between S - M for both


# Appendix Table S4: Relative body size corrected feeding rate
# contrasts of infection status effects by each clone (averaged over all diets)
mfeedmod_contrasts5 <- emmeans(mfeedmod2, specs = pairwise ~ Infection | Clone, type = "response")
mfeedmod_contrasts5      # This result tells us that there is a sig difference between uninfected and exposed for Standard but NOT for Mid37 (same for both feedingand relative feeding)




# Figures of Fungus Feeding rate
# Marginal effects for feeding rate
metsch_feed <- ggpredict(mfeedmod, c("Infection", "Diet", "Clone"))
plot(metsch_feed)

# figure with feeding rate using marginal effect - points are not positioned correctly
m_feeding_fig1 <- plot(metsch_feed, dodge = 0.6, add.data = T, show.title = F) +
  scale_color_manual(values = diet_colors) +
  #ggtitle("Relative Feeding Rate (mL per hr per mm)") +
  labs(x = "Parasite Exposure", y = "Feeding Rate (mL per hr per mm)") +
  theme_classic()
m_feeding_fig1
ggsave("figures/MetschFeedingRate_DietxClonexInfection2.tiff", plot = m_feeding_fig1, dpi = 600, width = 5, height = 4, units = "in", compression="lzw")


# Figure 1C: Fungus experiment Feeding Rate
# figure with boxplots and data points (week 1 feeding rate - NOT body size corrected)
m_feeding_fig <- ggplot(m_data_feeding, aes(x = Infection, y = Clearance.Week1)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width =0.9), show.legend = F) +
  geom_point(aes(color=Diet), size=1.5, position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.2), alpha = 0.4, show.legend = F) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninfected", "Exposed", "Infected"), labels = c("Uninf", "Exp", "Inf")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  labs(x = "Infection Status", y = bquote("Feeding Rate (mL per hr)")) +
  #ggtitle("Fungus Experiment") +
  theme_classic() + 
  theme(axis.text.x = element_text(size=9, color = "black"), axis.text.y = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
m_feeding_fig
ggsave("figures/MetschFeedingRate_DietxClonexInfection_Week1_NotBodysizeCorrected.tiff", plot = m_feeding_fig, dpi = 600, width = 5, height = 4, units = "in", compression="lzw")


# Figures of Fungus Relative Feeding rate (body size corrected)
# These figures are NOT in the manuscript
# Marginal effects for body size corrected feeding rate
metsch_feed <- ggpredict(mfeedmod2, c("Infection", "Diet", "Clone"))
plot(metsch_feed)

# figure with bodysize corrected feeding rate using marginal effects
m_rel_feeding_fig1 <- plot(metsch_feed, dodge = 0.6, add.data = T, show.title = F) +
  scale_color_manual(values = diet_colors) +
  #ggtitle("Relative Feeding Rate (mL per hr per mm)") +
  labs(x = "Parasite Exposure", y = "Relative Feeding Rate (mL per hr per mm)") +
  theme_classic()
m_rel_feeding_fig1
ggsave("figures/MetschFeedingRate_DietxClonexInfection2.tiff", plot = m_rel_feeding_fig1, dpi = 600, width = 5, height = 4, units = "in", compression="lzw")


# figure with boxplot and data points for bodysize corrected feeding rate
m_rel_feeding_fig <- ggplot(m_data_feeding, aes(x = Infection, y = Clearance_rel.Week1)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width =0.9)) +
  geom_point(aes(color=Diet), size=2, position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.2), alpha = 0.4) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninfected", "Exposed", "Infected"), labels = c("Uninf", "Exp", "Inf")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  #ggtitle("Growth Rate from Week 1 to 2 by infection status x diet x host clone for Metsch Expt") +
  labs(x = "Parasite Exposure", y = bquote("Relative Feeding Rate (mL per hr per mm)")) +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
m_rel_feeding_fig

ggsave("figures/MetschFeedingRate_DietxClonexInfection.tiff", plot = m_rel_feeding_fig, dpi = 600, width = 5, height = 4, units = "in", compression="lzw")





### Pasteuria feeding rate analyses

p_data_feeding <- p_data %>% filter(!is.na(Clearance.Week1))
hist(p_data_feeding$Clearance.Week1) # wow normal data!!
hist(p_data_feeding$Clearance_rel.Week1) # wow normal data!! 


# Feeding rate (not bodysize controlled)
pfeedmod <- lmer(Clearance.Week1 ~ Diet * Infection * Clone + (1|Block), data = p_data_feeding)
summary(pfeedmod)

# Table 1: Bacterium Experiment Feeding rate/exposure rate
Anova(pfeedmod)
AIC(pfeedmod)
plot(pfeedmod)
qqnorm(resid(pfeedmod))
qqline(resid(pfeedmod))         # looks good
shapiro.test(resid(pfeedmod))   # significant! But doesn't look too severe from the qqnorm plots (see https://stats.stackexchange.com/questions/32957/testing-normality-assumptions-for-linear-mixed-models-and-mixed-repeated-glm-a)


# Appendix Table S5
# contrasts of diet effects by each infection status (averaged over both clones)
pfeedmod_contrasts <- emmeans(pfeedmod, specs = pairwise ~ Diet | Infection + Clone, type = "response")
pfeedmod_contrasts 

# Appendix Table S4: Bacterium Experiment Feeding rate/exposure rate
# contrasts of infection status effects by each clone (averaged over all diets)
pfeedmod_contrasts2 <- emmeans(pfeedmod, specs = pairwise ~ Infection | Clone, type = "response")
pfeedmod_contrasts2

pfeedmod_infstatus_est <- as.data.frame(pfeedmod_contrasts2)
pfeedmod_infstatus_est <- filter(pfeedmod_infstatus_est, contrast == ".")
pfeedmod_infstatus_est$Infection <- factor(pfeedmod_infstatus_est$Infection, levels = c("Uninfected", "Exposed", "Infected"))



# Relative Feeding rate (bodysize controlled) - analysis described in Appendix Section S1.7
pfeedmod2 <- lmer(Clearance_rel.Week1 ~ Diet * Infection * Clone +(1|Block), data = p_data_feeding)
summary(pfeedmod2)

# Appendix Table S2
Anova(pfeedmod2)
AIC(pfeedmod2)
plot(pfeedmod2)
qqnorm(resid(pfeedmod2))       # looks good
qqline(resid(pfeedmod2))
shapiro.test(resid(pfeedmod2)) # significant! But doesn't look too severe from the qqnorm plots (see https://stats.stackexchange.com/questions/32957/testing-normality-assumptions-for-linear-mixed-models-and-mixed-repeated-glm-a)


pfeedmod3 <- lmer(Clearance_rel.Week1 ~ Diet + Infection + Clone + Diet:Clone +(1|Block), data = p_data_feeding)
summary(pfeedmod3)
Anova(pfeedmod3)
AIC(pfeedmod3)
plot(pfeedmod3)
qqnorm(resid(pfeedmod3))     # looks good
qqline(resid(pfeedmod3))
shapiro.test(resid(pfeedmod3))  # not significant
# though this simplifies the model, it makes predicted effects across clone x infection and diet x infection the same
# prefer to use the pfeedmod2 above over this one


# tests contrasts for each factor holding the other two constant in separate tests
# contrasts of diet effects by each infection status and clone. 
pfeedmod2_contrasts <- emmeans(pfeedmod2, specs = pairwise ~ Diet | Infection + Clone, type = "response")
pfeedmod2_contrasts# Mid37 there is a dif between S - M for uninfected and SM - M+ for exposed, and no differences for infected; and for Std uninfected there are dif between S - M, S - M+, and SM - M, for exposed and infected dif between S - M for both

# Appendix Table S4
# contrasts of infection status effects by each clone (averaged over all diets)
pfeedmod2_contrasts2 <- emmeans(pfeedmod2, specs = pairwise ~ Infection | Clone, type = "response")
pfeedmod2_contrasts2   # This result tells us that there is a sig difference between uninfected and exposed for Standard but NOT for Mid37 (same for both feeding and relative feeding)



##### figure with feeding rate (NOT body size corrected) - NOT shown in the manuscript
past_feed <- ggpredict(pfeedmod, c("Infection", "Diet", "Clone"))
plot(past_feed)

# figure with feeding rate (not body size corrected)
p_feeding_fig1 <- plot(past_feed, dodge = 0.6, add.data = T, show.title = F) +
  scale_color_manual(values = diet_colors) +
  #ggtitle("Relative Feeding Rate (mL per hr per mm)") +
  labs(x = "Parasite Exposure", y = "Feeding Rate (mL per hr)") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
p_feeding_fig1
ggsave("figures/PastFeedingRate_DietxClonexInfection.tiff", plot = p_feeding_fig1, dpi = 600, width = 5, height = 4, units = "in", compression="lzw")



# figure for Meg (week 1 feeding rate - NOT body size corrected)
p_feeding_fig <- ggplot(p_data_feeding, aes(x = Infection, y = Clearance.Week1)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width =0.9)) +
  geom_point(aes(color=Diet), size=1.5, position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.2), alpha = 0.4, show.legend = F) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninfected", "Exposed", "Infected"), labels = c("Uninf", "Exp", "Inf")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  labs(x = "Infection Status", y = bquote("Feeding Rate (mL per hr)")) +
  #ggtitle("Bacteria Experiment") +
  theme_classic() + 
  theme(axis.text.x = element_text(size=9, color = "black"), axis.text.y = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
p_feeding_fig

ggsave("figures/PastFeedingRate_DietxClonexInfection_Week1_NotBodysizeCorrected.tiff", plot = p_feeding_fig, dpi = 600, width = 5, height = 4, units = "in", compression="lzw")




##### # figure with feeding rate (body size corrected)
past_feed2 <- ggpredict(pfeedmod3, c("Infection", "Diet", "Clone"))
plot(past_feed)

# figure with bodysize corrected feeding rate using marginal effects
p_feeding_fig2 <- plot(past_feed2, dodge = 0.6, add.data = T, show.title = F) +
  scale_color_manual(values = diet_colors) +
  #ggtitle("Relative Feeding Rate (mL per hr per mm)") +
  labs(x = "Parasite Exposure", y = "Relative Feeding Rate (mL per hr per mm)") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
p_feeding_fig2
ggsave("figures/PastRelFeedingRate_DietxClonexInfection.tiff", plot = p_feeding_fig2, dpi = 600, width = 5, height = 4, units = "in", compression="lzw")


# figure with bodysize corrected feeding rate using boxplot and data points
p_feeding2 <- ggplot(p_data_feeding, aes(x = Infection, y = Clearance_rel.Week1)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width =0.9)) +
  geom_point(aes(color=Diet), size=2, position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.2), alpha = 0.4) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninfected", "Exposed", "Infected")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  #ggtitle("Growth Rate from Week 1 to 2 by infection status x diet x host clone for Metsch Expt") +
  labs(x = "Infection Status", y = bquote("Relative Feeding Rate (mL per hr per mm)")) +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
p_feeding2

ggsave("figures/PastFeedingRate_DietxClonexInfection.tiff", plot = p_feeding2, dpi = 600, width = 5, height = 4, units = "in", compression="lzw")








## Mature parasite spore yield and size --------------------------------------------


### Comparison of Mature Bacterium spore yield based on diet and host clone.
# Fungus Mature spore yield 

# remove animals that died during exposure and remove uninfected treatments and remove animals with no spores detected
m_data_spores <- m_data %>% filter(!is.na(InfectionStatus) & Parasites == "Metsch" & total_spores > 0)

# filter by mature spores only
m_data_maturespores <- m_data_spores %>% filter(mature_spores > 0) 

mmaturesporemod <- glmmTMB(mature_spores ~ Diet + Clone + (1|Block), family = nbinom2(), data = m_data_maturespores)  # lower AIC than nbinom1, otherwise fit is similar
summary(mmaturesporemod) # not enough data to run the Diet x Clone interaction

# Table 1: Fungus Mature spore yield
Anova(mmaturesporemod)
AIC(mmaturesporemod)
testDispersion(mmaturesporemod)
testZeroInflation(mmaturesporemod)
mmaturespore_simResid <- simulateResiduals(fittedModel = mmaturesporemod)
plot(mmaturespore_simResid)

mmaturesporemod_contrasts <- emmeans(mmaturesporemod, specs = pairwise ~ Diet | Clone, type = "response")
mmaturesporemod_contrasts

range(m_data_spores$mature_spores)



# Marginal effects from best model with Abundance on the x axis
metsch_maturespores <- ggpredict(mmaturesporemod, c("Clone", "Diet"))


# figure with marginal effects (not in manuscript)
m_maturespore_fig1 <- plot(metsch_maturespores, dodge = 0.5, add.data = T, show.title = F) +
  scale_color_manual(values = diet_colors) +
  labs(x = "Host Clone", y = "Mature Fungus Spore Yield (log scale)") +
  scale_y_log10(labels = scales::comma) +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
m_maturespore_fig1
ggsave("figures/MetschMatureSpores_DietxClone.tiff", plot = m_maturespore_fig1, dpi = 600, width = 4, height = 3, units = "in", compression="lzw")


# Figure 1E: Fungus Mature spore yield
# figure with box plots and points
m_maturespore_fig <- ggplot(m_data_maturespores, aes(x = Diet, y = mature_spores)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width=0.7)) +
  geom_point(aes(color=Diet), size=2, position=position_jitterdodge(dodge.width=0.7, jitter.width = 0.7), alpha = 0.4) +
  facet_wrap(~Clone) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  ylim(0,300000) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Diet Treatment", y = "Fungus Spore Yield (log scaled)") +
  theme_classic()  +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
m_maturespore_fig
ggsave("figures/MetschMatureSpores_DietxClone2.tiff", plot = m_maturespore_fig, dpi = 600, width = 4, height = 3, units = "in", compression="lzw")


# # Results presented in Appendix 
# S1.5 Results for mature spore yield and microcystin-LR concentration
#mmaturesporemod2 <- glmmTMB(mature_spores ~ microcystin.avg + Clone + Diet + (1|Block), family = nbinom1(), data = m_data_maturespores)  # Model can't run with clone included, probably because there are some diets missing for Mid37
mmaturesporemod2 <- glmmTMB(mature_spores ~ microcystin.avg + Diet + (1|Block), family = nbinom1(), data = m_data_maturespores)
overdisp_fun(mmaturesporemod2)
summary(mmaturesporemod2)

# Table S1 in Appendix Section S1.5
Anova(mmaturesporemod2)  
AIC(mmaturesporemod2)
qqnorm(resid(mmaturesporemod2))
qqline(resid(mmaturesporemod2))      #### NOTE FROM MAD: This isn't running for me (This worked for Michelle 4/24)
testDispersion(mmaturesporemod2)
testZeroInflation(mmaturesporemod2)
mmaturespore2_simResid <- simulateResiduals(fittedModel = mmaturesporemod2)
plot(mmaturespore2_simResid)

# try model without M+ samples (outliers b/c so much higher than all other samples, see if there is a pattern )
m_data_spores_noM_tox <- filter(m_data_spores, Diet != "M+")
mmaturesporemod3 <- glmmTMB(mature_spores ~ microcystin.avg + Diet + Clone + (1|Block), family = nbinom2(), data = m_data_spores_noM_tox)
summary(mmaturesporemod3)
Anova(mmaturesporemod3)







### Comparison of Mature Bacterium spore yield based on diet and host clone.

# remove animals that died during exposure and remove uninfected treatments and remove animals with no spores detected
p_data_spores <- p_data %>% filter(!is.na(InfectionStatus) & Parasites == "Pasteuria" & total_spores > 0) 

# remove the one sample without mature spores
p_data_spores_mature <- filter(p_data_spores, mature_spores > 0) 


pmaturesporemod <- glmmTMB(mature_spores ~ Diet * Clone + (1|Block), family = nbinom1(), data = p_data_spores_mature)  # best fit of the residuals
summary(pmaturesporemod)

#Table 1: Bacterium mature spore yield
Anova(pmaturesporemod)
AIC(pmaturesporemod)
testDispersion(pmaturesporemod)
testZeroInflation(pmaturesporemod)
pmaturespore_simResid <- simulateResiduals(fittedModel = pmaturesporemod)
plot(pmaturespore_simResid)  

pmaturesporemod_contrasts <- emmeans(pmaturesporemod, specs = pairwise ~ Diet | Clone, type = "response")
pmaturesporemod_contrasts

range(p_data_spores$mature_spores)


# Marginal effects from best model with Abundance on the x axis
past_maturespores <- ggpredict(pmaturesporemod, c("Clone", "Diet"))


# figure with marginal effects (not in manuscript)
p_maturespore_fig3 <- plot(past_maturespores, dodge = 0.5, add.data = T, show.title = F) +
  scale_color_manual(values = diet_colors) +
  #ggtitle("Mature Pasteuria Spore Yield by diet x host clone") +
  labs(x = "Host Clone", y = "Mature Bacterium Spore Yield (log scale)") +
  scale_y_log10(labels = scales::comma, breaks = c(100, 1000, 10000,100000,1000000)) +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=9, color="black"))
p_maturespore_fig3

ggsave("figures/PastMatureSpores_DietxClone.tiff", plot = p_maturespore_fig3, dpi = 600, width = 4, height = 3, units = "in", compression="lzw")


# Figure 1F: Bacterium mature spore yield 
# figure with box plots and points of spore yield
p_maturespore_fig <- ggplot(p_data_spores_mature, aes(x = Diet, y = mature_spores)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width=0.7), show.legend = F) +
  geom_point(aes(color=Diet), size=2, position=position_jitterdodge(dodge.width=0.7, jitter.width = 0.7), alpha = 0.4, show.legend = F) +
  facet_wrap(~Clone) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  ylim(0,300000) +
  scale_y_log10(labels = scales::comma) +
  #scale_y_log10(labels = scales::comma, breaks = c(100, 1000, 10000,100000,1000000)) +
  #ggtitle("Mature Pasteuria Spore Yield by diet x host clone") +
  labs(x = "Diet Treatment", y = "Bacterium spore yield (log scaled)") +
  theme_classic() +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
p_maturespore_fig
ggsave("figures/PastMatureSpores_DietxClone2.tiff", plot = p_maturespore_fig, dpi = 600, width = 4, height = 3, units = "in", compression="lzw")




# # Results presented in Appendix 
# S1.5 Results for mature spore yield and microcystin-LR concentration
pmaturesporemod2 <- glmmTMB(mature_spores ~ microcystin.avg + Diet + Clone + (1|Block), family = nbinom1(), data = p_data_spores_mature)  # lower AIC model
summary(pmaturesporemod2)

# Table S1 in appendix in Section S1.5
Anova(pmaturesporemod2)   
AIC(pmaturesporemod2)
testDispersion(pmaturesporemod2)
testZeroInflation(pmaturesporemod2)
pmaturespore2_simResid <- simulateResiduals(fittedModel = pmaturesporemod2)
plot(pmaturespore2_simResid)  

# Additional versions of the analysis below (not in the manuscript)
# try model without M+ samples (outliers b/c so much higher than all other samples, see if there is a pattern )
p_data_spores_noM_tox <- filter(p_data_spores_mature, Diet != "M+")
pmaturesporemod3 <- glmmTMB(mature_spores ~ microcystin.avg * Diet + (1|Block), family = nbinom1(), data = p_data_spores_noM_tox)
summary(pmaturesporemod3)
Anova(pmaturesporemod3)
testDispersion(pmaturesporemod3)
testZeroInflation(pmaturesporemod3)
pmaturespore3_simResid <- simulateResiduals(fittedModel = pmaturesporemod3)
plot(pmaturespore3_simResid)  

# Remove other diets to only compare M and M+ diets  (only 8 data points...)
p_data_spores_M_tox <- filter(p_data_spores_mature, Diet == "M" | Diet == "M+")
pmaturesporemod4 <- glmmTMB(mature_spores ~ log(microcystin.avg+0.01) * Diet + (1|Block), family = nbinom1(), data = p_data_spores_M_tox)
summary(pmaturesporemod4)
Anova(pmaturesporemod4)
testDispersion(pmaturesporemod4)
testZeroInflation(pmaturesporemod4)
pmaturespore4_simResid <- simulateResiduals(fittedModel = pmaturesporemod4)
plot(pmaturespore4_simResid)  



### Bacterium spore area (size)
p_data_spores_size <- filter(p_data_spores, !is.na(AvgSporeSize)) # remove animals where we did not measure spore size

# initial model of spore size based on diet x clone
sporesize_mod <- lm(AvgSporeSize ~ Diet * Clone, data = p_data_spores_size)
summary(sporesize_mod)
Anova(sporesize_mod)
AIC(sporesize_mod)
qqnorm(resid(sporesize_mod))
qqline(resid(sporesize_mod)) # looks ok, not the best   #### NOTE FROM MAD: This isn't running for me (it worked for Michelle 4/24)
ols_test_normality(sporesize_mod)  # shapiro-wilk is not significant, normal distribution should be fine
ols_test_breusch_pagan(sporesize_mod)
ols_plot_resid_fit(sporesize_mod)

# figure of spores size vs diet x clone
sporesize_plot <- ggplot(p_data_spores_size, aes(x = Diet, y = AvgSporeSize, fill = Diet)) +
  geom_boxplot() +
  geom_jitter(size = 2, alpha = 0.5, width = 0.2) +
  facet_wrap(~Clone) +
  scale_fill_manual(values = diet_colors) +
  #ggtitle("Pasteuria spore size based on diet") +
  labs(x = "Diet", y = bquote("Mean Bacterium Spore Area " ~ (mu *m^2))) +
  theme_classic()
sporesize_plot
ggsave(here("figures/Pasteuria_sporesize_dietxclone.tiff"), plot = sporesize_plot, dpi = 300, width = 5, height = 4, units = "in", compression="lzw")

# very few replicates among the standard diet treatments, not a robust test including host clone
# remove clone and run model with just diet


# Spore size model with just diet (analysis in Manuscript)
sporesize_mod2 <- lm(AvgSporeSize ~ Diet, data = p_data_spores_size)
summary(sporesize_mod2)

# Table 1: Bacterium spore size (surface area)
Anova(sporesize_mod2) # reported results
AIC(sporesize_mod2)
qqnorm(resid(sporesize_mod2))
qqline(resid(sporesize_mod2)) # looks ok, not the best   #### NOTE FROM MAD: This isn't running for me (it worked for Michelle 4/24)
ols_test_normality(sporesize_mod2)  # shapiro-wilk is not significant, normal distribution should be fine
ols_test_breusch_pagan(sporesize_mod2)
ols_plot_resid_fit(sporesize_mod2)

psporesize_contrasts <- emmeans(sporesize_mod2, specs = pairwise ~ Diet, type = "response")
psporesize_contrasts


# Figure 3B: Bacterium spore Area
# figure of spores size vs diet
sporesize_plot2 <- ggplot(p_data_spores_size, aes(x = Diet, y = AvgSporeSize, fill = Diet)) +
  geom_boxplot(show.legend = F) +
  geom_jitter(size = 2, alpha = 0.5, width = 0.2, show.legend = F) +
  scale_fill_manual(values = diet_colors) +
  labs(x = "Diet Treatment", y = bquote("Mean Bacterium Spore Area " ~ (mu *m^2))) +
  #annotate(geom = "text", x = 3.5, y = 23, label = "p = 0.039") +
  #annotate(geom = "text", x = c(1,2,3,4), y = c(24.5,21,21,21), label = c("a", "ab", "ab", "b")) +
  ggtitle("Bacterium Experiment") +
  theme_classic()  +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
sporesize_plot2
ggsave(here("figures/Pasteuria_sporesize_diet.tiff"), plot = sporesize_plot2, dpi = 300, width = 5, height = 4, units = "in", compression="lzw")



### Full Figure 1 in Manuscript

Figure1 <- ggarrange(m_prev_fig, p_prev_fig, m_feeding_fig, p_feeding_fig, m_maturespore_fig, p_maturespore_fig, labels = "auto",
                     ncol = 2, nrow = 3, legend = "right", common.legend = T)
Figure1
ggsave("figures/manuscript/Fig1_Prev_Exposure_Yield.tiff", plot = Figure1, dpi = 600, width = 7.5, height = 9, units = "in", compression="lzw")







## Gut spore data for Metsch --------------------------------------------------

m_data_gutspore <- m_data %>% filter(Block == 5 & Infection != "Uninfected") 
# only include animals in block 5 that were checked for gut spores and remove 3 animals that did not get exposed to parasites (STD SM, M, and M+ rep 20 for all, b/c ran out of spores for exposure)
hist(m_data_gutspore$TotalGutSpores)


### Analysis of Fungus spores attacking the gut
mgutsporemod <- glmmTMB(TotalGutSpores ~ Diet * Clone, family = poisson, data = m_data_gutspore)  # lowest AIC, best model use this one
# Note: left off the block random effect because all these animals were in the same block 
summary(mgutsporemod)
Anova(mgutsporemod)      ## Use this model for results in the manuscript
overdisp_fun(mgutsporemod)
AIC(mgutsporemod)
testDispersion(mgutsporemod)
testZeroInflation(mgutsporemod)
mgutspore_simResid <- simulateResiduals(fittedModel = mgutsporemod)
plot(mgutspore_simResid)

mgutsporemod_contrasts <- emmeans(mgutsporemod, specs = pairwise ~ Diet | Clone, type = "response")
mgutsporemod_contrasts


# Figure 3A: Fungus spores attacking gut
# figure with boxplot and data points
m_totgutspores_fig <- ggplot(m_data_gutspore, aes(x = Diet, y = TotalGutSpores)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width=0.7)) +
  geom_point(aes(color=Diet), size=2, position=position_jitterdodge(dodge.width=0.7, jitter.width = 0.7), alpha = 0.4) +
  facet_wrap(~Clone) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  labs(x = "Diet Treatment", y = "Number of Fungus Spores Attacking Gut") +
  ggtitle("Fungus Experiment") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
m_totgutspores_fig
ggsave("figures/MetschGutSpores_DietxClone.tiff", plot = m_totgutspores_fig, dpi = 600, width = 5, height = 4, units = "in", compression="lzw")





# Additional version of the model that is not included in the manuscript
## Check if body size impacts the likelihood of spores attacking the gut
mgutsporemod2 <- glm(TotalGutSpores ~ Diet * Clone * length.um, family = poisson, data = m_data_gutspore)
# Note: left off the block random effect because all these animals were in the same block 
# glmmTMB model didn't run but the lme4 glm did
summary(mgutsporemod2)
Anova(mgutsporemod2)   # host body size does not impact likelihood that spores attack gut, while accounting for diet and clone
overdisp_fun(mgutsporemod2)
AIC(mgutsporemod2)     # AIC is much higher with body size is included in the model compared to the initial model above
testDispersion(mgutsporemod2)
testZeroInflation(mgutsporemod2)
mgutspore2_simResid <- simulateResiduals(fittedModel = mgutsporemod)
plot(mgutspore2_simResid)

# model has much higher AIC than previous models
mgutsporemod3 <- glmmTMB(TotalGutSpores ~ Clone * length.um, family = poisson, data = m_data_gutspore)  # much higher AIC than the initial model
summary(mgutsporemod3)
Anova(mgutsporemod3)







### Analysis of the number of Hemocytes (Not presented in manuscript)
m_data_hemocytes <- m_data_gutspore %>% filter(TotalPunctSpores > 0)  # we only counted hemocytes where there were penetrated spores
dim(m_data_hemocytes)
m_data_hemocytes$Unique.code <- as.factor(m_data_hemocytes$Unique.code)

hist(m_data_hemocytes$Hemocytes)
hist(m_data_hemocytes$HemocytesPerSpore)


## Number of hemocytes (raw counts of hemocytes)
mhemocytemod <- glmmTMB(Hemocytes ~ Diet + Clone + (1|Unique.code), family = poisson(), data = m_data_hemocytes)  # lowest AIC model
#mhemocytemod <- glmmTMB(Hemocytes ~ Diet + Clone, family = nbinom2(), data = m_data_hemocytes)                    # very similar AIC
#mhemocytemod <- glmmTMB(Hemocytes ~ Diet + Clone, family = nbinom1(), data = m_data_hemocytes)                    # not as good AIC as above
summary(mhemocytemod)
Anova(mhemocytemod)
overdisp_fun(mhemocytemod) 
AIC(mhemocytemod)
testDispersion(mhemocytemod)
testZeroInflation(mhemocytemod)
mhemocyte_simResid <- simulateResiduals(fittedModel = mhemocytemod)
plot(mhemocyte_simResid)
# model is initially overdispersed, so I included a unique ID random effect to correct for that. The model did not have enough replicates 
# in each diet x clone treatment to be able to run the interaction.

# Hemocytes and including body size as a predictor (not in manuscript)
mhemocytemod2 <- glmmTMB(Hemocytes ~ Diet + Clone * length.um + Diet:length.um + (1|Unique.code), family = poisson(), data = m_data_hemocytes) 
# Note: left off the block random effect because all these animals were in the same block 
summary(mhemocytemod2)
Anova(mhemocytemod2)
AIC(mhemocytemod2)   # model AIC is better when length is included, but it is not significant
testDispersion(mhemocytemod2)
testZeroInflation(mhemocytemod2)
mhemocyte2_simResid <- simulateResiduals(fittedModel = mhemocytemod2)
plot(mhemocyte2_simResid)



# Marginal effects 
metsch_hemocytes <- ggpredict(mhemocytemod, c("Diet", "Clone"))

# Figures NOT in manuscript
# figure with marginal effects
m_hemocytes_fig <- plot(metsch_hemocytes, dodge = 0.4, add.data = T) +
  scale_color_manual(values = diet_colors) +
  scale_fill_manual(values = diet_colors) +
  labs(x = "Diet", y = "Recruited Hemocytes") +
  ylim(0,55) +
  theme_classic()
m_hemocytes_fig

## Figure with boxplot and data points
m_hemocytes_fig2 <- ggplot(m_data_hemocytes, aes(x = Diet, y = Hemocytes)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width=0.7)) +
  geom_point(aes(color=Diet), size=2, position=position_jitterdodge(dodge.width=0.7, jitter.width = 0.7), alpha = 0.4) +
  facet_wrap(~Clone) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_shape_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  #ggtitle("Hemocytes per attacking Metsch spore in the gut by diet x host clone") +
  labs(x = "Diet", y = "Recruited Hemocytes") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=9, color="black"))
m_hemocytes_fig2

ggsave("figures/Hemocytes_DietxClone.tiff", plot = m_hemocytes_fig2, dpi = 600, width = 5, height = 4, units = "in", compression="lzw")






## Since we know from above that the number of attacking spores per host also differ significantly, we modeled 
## The number of hemocytes per attacking spore


# Hemocytes per spore and NOT including body size as a predictor 
# [USE THIS ONE IN THE MANUSCRIPT]

mhemocytemod3 <- glmmTMB(Hemocytes ~ Diet + Clone + (1|Unique.code), offset = log(TotalPunctSpores), family = poisson(), data = m_data_hemocytes) 
# Note: left off the block random effect because all these animals were in the same block
# Not enough diet x clone replicates to include the interaction
# Adding body size made the predicted effects come out very strange, and I did not trust that model. Model is robust without body size
summary(mhemocytemod3)

# Table 1: Hemocytes per fungus spore
Anova(mhemocytemod3)
AIC(mhemocytemod3) # model AIC is better 
testDispersion(mhemocytemod3)
testZeroInflation(mhemocytemod3)
mhemocyte3_simResid <- simulateResiduals(fittedModel = mhemocytemod3)
plot(mhemocyte3_simResid)

mhemocytemod_contrasts <- emmeans(mhemocytemod3, specs = pairwise ~ Diet | Clone, type = "response", offset = T)
mhemocytemod_contrasts

mhemocytemod_contrasts2 <- emmeans(mhemocytemod3, specs = pairwise ~ Clone, type = "response", offset = T)
mhemocytemod_contrasts2

# Hemocytes per spore and WITH including body size as a predictor
mhemocytemod3 <- glmmTMB(Hemocytes ~ Diet + Clone + length.um + Diet:length.um + (1|Unique.code), offset = log(TotalPunctSpores), family = poisson(), data = m_data_hemocytes)  # lower AIC than above...


## Test to show that hemoctyes are responding to the intensity of exposure to attacking spores
mhemocytemod4 <- glmmTMB(Hemocytes ~ TotalPunctSpores + (1|Unique.code), family = poisson(), data = m_data_hemocytes) # This model has severe quantile deviations...needs diet and clone in the model to fix, but then total spores is not sig. 
mhemocytemod4 <- glmmTMB(Hemocytes ~ TotalPunctSpores + Diet + Clone + (1|Unique.code), family = poisson(), data = m_data_hemocytes) 
summary(mhemocytemod4)
Anova(mhemocytemod4)
testDispersion(mhemocytemod4)
testZeroInflation(mhemocytemod4)
mhemocyte4_simResid <- simulateResiduals(fittedModel = mhemocytemod4)
plot(mhemocyte4_simResid)



# Marginal effects
metsch_hemocytesperspore <- ggpredict(mhemocytemod3, c("Clone", "Diet"))
plot(metsch_hemocytesperspore) + ylim(0,15)


# Figure 3C: Hemocytes per attacking fungus spore
# figure boxoplot and data points
m_hemocytesperspore_fig <- ggplot(m_data_hemocytes, aes(x = Diet, y = HemocytesPerSpore)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width=0.7)) +
  geom_point(aes(color=Diet), size=2, position=position_jitterdodge(dodge.width=0.7, jitter.width = 0.7), alpha = 0.4) +
  facet_wrap(~Clone) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  #ggtitle("Hemocytes per attacking Metsch spore in the gut by diet x host clone") +
  labs(x = "Diet Treatment", y = "Hemocytes per Attacking Fungus Spore") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
m_hemocytesperspore_fig

ggsave("figures/HemocytesPerSpore_DietxClone.tiff", plot = m_hemocytesperspore_fig, dpi = 600, width = 5, height = 4, units = "in", compression="lzw")





# Full Figure 3 in Manuscript

Figure3 <- ggarrange(m_totgutspores_fig, sporesize_plot2, m_hemocytesperspore_fig, labels = "auto",
                     ncol = 2, nrow = 2, legend = "right", common.legend = T)
Figure3
ggsave("figures/manuscript/Fig3_GutSpores_Hemo_SporeSize.tiff", plot = Figure3, dpi = 600, width = 7, height = 7, units = "in", compression="lzw")




