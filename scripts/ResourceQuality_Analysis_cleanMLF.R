# Resource Quality Experiment Analysis

# This code runs all of the analyses included in the main text and appendix of Fearon et al.
# "Resource Quality Differentially Impacts Daphnia Interactions with Two Parasites"
# submitted to Ecological Monographs.

# Code Written by: Michelle L Fearon
# Last updated: 7/15/2024
# Last update by: Michelle L Fearon


# Terminology Used
      # Bacterium = Pasteuria ramosa
            # data sets for Pasteuria (Bacteria) are labeled with a "p"
      # Fungus = Metschnikowia bicuspidata (Metsch)
            # data sets for Metsch (Fungus) are labeled with a "m"




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


# load Pasteuria data
p_data <- read.csv(here("data/ResourceQuality_Pasteuria_Full.csv"), stringsAsFactors = F)
p_longdata <- read.csv(here("data/ResourceQuality_Pasteuria_Full_ByExptWeek.csv"), stringsAsFactors = F)


# load Metsch data
m_data <- read.csv(here("data/ResourceQuality_Metsch_Full.csv"), stringsAsFactors = F)
m_longdata <- read.csv(here("data/ResourceQuality_Metsch_Full_ByExptWeek.csv"), stringsAsFactors = F)



# Reorder factors
p_data$Diet <- factor(p_data$Diet, levels = c("S", "SM", "M", "M+"))
p_data$Infection <- factor(p_data$Infection, levels = c("Uninfected", "Exposed", "Infected"))
p_data$Parasites <- factor(p_data$Parasites, levels = c("Uninfected", "Pasteuria"))
p_longdata$Diet <- factor(p_longdata$Diet, levels = c("S", "SM", "M", "M+"))
p_longdata$Infection <- factor(p_longdata$Infection, levels = c("Uninfected", "Exposed", "Infected"))

m_data$Diet <- factor(m_data$Diet, levels = c("S", "SM", "M", "M+"))
m_data$Infection <- factor(m_data$Infection, levels = c("Uninfected", "Exposed", "Infected"))
m_data$Parasites <- factor(m_data$Parasites, levels = c("Uninfected", "Metsch"))
m_longdata$Diet <- factor(m_longdata$Diet, levels = c("S", "SM", "M", "M+"))
m_longdata$Infection <- factor(m_longdata$Infection, levels = c("Uninfected", "Exposed", "Infected"))


# source functions
source(here("scripts/functions.R"))

# set the color scheme for the Diet treatments
diet_colors <- c("#ADDD8E", "#41AB5D", "#006837", "#1D91C0")



# Parasite infectivity and host susceptibility ----------------------------------

### Comparison of Metsch infection prevalence based on diet, and host clone.

# remove animals that died during exposure and remove uninfected treatments
m_data_prev <- m_data %>% filter(!is.na(InfectionStatus) & Parasites == "Metsch") 

m_data_prev$Block <- as.factor(m_data_prev$Block)

m_N_per_treatment_lifespand <- m_data %>%
  count(Diet, Parasites, Clone, Block, lifespan.days) %>%
  arrange(Diet, Parasites, Clone, lifespan.days)

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


p_N_per_treatment_lifespand <- p_data %>%
  count(Diet, Parasites, Clone, Block, lifespan.days) %>%
  arrange(Diet, Parasites, Clone, lifespan.days)

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



#### Figure S3 in Appendix
FigS3 <- ggarrange(m_prev_feed_fig, p_prev_feed_fig, nrow = 1, ncol = 2, labels = "auto", common.legend = T, legend = "right")
FigS3
ggsave(here("figures/manuscript/FigS3_Prevalence_FeedingRate.tiff"), plot = FigS3, dpi = 600, width = 8, height = 3.5, units = "in", compression="lzw")




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


# Appendix Table S4: Bacterium Experiment Relative Feeding rate/exposure rate
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
qqline(resid(sporesize_mod)) # looks ok, not the best 
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
qqline(resid(sporesize_mod2)) # looks ok, not the best
ols_test_normality(sporesize_mod2)  # shapiro-wilk is not significant, normal distribution should be fine
ols_test_breusch_pagan(sporesize_mod2)
ols_plot_resid_fit(sporesize_mod2)

psporesize_contrasts <- emmeans(sporesize_mod2, specs = pairwise ~ Diet, type = "response")
psporesize_contrasts


# Figure 2B: Bacterium spore Area
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


# There is NO strong correlation between spore size and number of mature spores
cor.test(p_data_spores_size$AvgSporeSize, p_data_spores_size$mature_spores)
cor.test(p_data_spores_size$AvgSporeSize, log(p_data_spores_size$mature_spores))

sporesize_yield_plot <- ggplot(p_data_spores_size, aes(x = mature_spores, y = AvgSporeSize)) +
  #geom_boxplot(show.legend = F) +
  geom_jitter(aes(color = Diet), size = 2, width = 0.2, show.legend = F) +
  geom_smooth(method = lm, color = "black") +
  scale_color_manual(values = diet_colors) +
  labs(x = "Mature Bacterium Spore Yield (log scale)", y = bquote("Mean Bacterium Spore Area " ~ (mu *m^2))) +
  scale_x_log10(labels = scales::comma, breaks = c(100, 1000, 10000,100000,1000000)) +
  #annotate(geom = "text", x = 3.5, y = 23, label = "p = 0.039") +
  #annotate(geom = "text", x = c(1,2,3,4), y = c(24.5,21,21,21), label = c("a", "ab", "ab", "b")) +
  ggtitle("Bacterium Experiment") +
  theme_classic()  +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
sporesize_yield_plot



### Full Figure 1 in Manuscript

Figure1 <- ggarrange(m_prev_fig, p_prev_fig, m_feeding_fig, p_feeding_fig, m_maturespore_fig, p_maturespore_fig, labels = "auto",
                     ncol = 2, nrow = 3, legend = "right", common.legend = T)
Figure1
ggsave("figures/manuscript/Fig1_Prev_Exposure_Yield.tiff", plot = Figure1, dpi = 600, width = 7.5, height = 9, units = "in", compression="lzw")







## Gut spore data for Metsch --------------------------------------------------

m_data_gutspore <- m_data %>% filter(Block == 5 & Infection != "Uninfected") 
# only include animals in block 5 that were checked for gut spores and remove 3 animals that did not get exposed to parasites (STD SM, M, and M+ rep 20 for all, b/c ran out of spores for exposure)
hist(m_data_gutspore$TotalGutSpores)

# calculate N animals per treatment that died due to procedure to check for gut spores
M_gutspore_mortality <- m_data_gutspore %>%
  filter(Lifespan <= 9) %>%
  count(Diet, Clone)


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


# Figure 2A: Fungus spores attacking gut
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


# Figure 2C: Hemocytes per attacking fungus spore
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





# Full Figure 2 in Manuscript

Figure2 <- ggarrange(m_totgutspores_fig, sporesize_plot2, m_hemocytesperspore_fig, labels = "auto",
                     ncol = 2, nrow = 2, legend = "right", common.legend = T)
Figure2
ggsave("figures/manuscript/Fig2_GutSpores_Hemo_SporeSize.tiff", plot = Figure3, dpi = 600, width = 7, height = 7, units = "in", compression="lzw")



# Host vital rates: impacts of diet and infection --------------------------------


## Total Fecundity ----------------------------------------------------------------------


### Total Fecundity by treatment for Metsch experiment

m_data_fecund <- m_data %>% 
  filter(!is.na(Total.Babies) & lifespan.days > 11) %>%  # first babies observed on day 11
  mutate(Parasite_Treatment = factor(if_else(Parasites == "Metsch", "Inoculated", "Uninoculated"), levels = c("Uninoculated", "Inoculated")))
hist(m_data$Total.Babies) 


# calculate N animals that were excluded 
m_data_fecund_exluded <- m_data %>% 
  filter(is.na(Total.Babies) | lifespan.days <= 11) %>% # first babies observed on day 11
  count(Diet, Parasites, Clone, Block)

m_total_n_excluded <- sum(m_data_fecund_exluded$n)

#mfecundmod <- glmmTMB(Total.Babies ~ Diet + Infection + Clone + Diet:Infection + Diet:Clone + Infection:Clone + (1|Block) + (1|Unique.code), family = poisson, data = m_data_fecund)  # overdispersed unless include unique code random effect
#mfecundmod <- glmmTMB(Total.Babies ~ Diet + Infection + Clone + Diet:Infection + Diet:Clone + Infection:Clone + (1|Block), family = nbinom2(), data = m_data_fecund)  #higher AIC than quasipoisson, still had overdispersion
mfecundmod <- glmmTMB(Total.Babies ~ Diet + Infection + Clone + Diet:Infection + Diet:Clone + Infection:Clone + (1|Block), family = nbinom1(), data = m_data_fecund)   # lower AIC than neg binomial and corrects for overdispersion, use this model
summary(mfecundmod)

# Table 1: Fungus Total Fecundity 
Anova(mfecundmod) # model doesn't have enough power for the 3 way interaction
AIC(mfecundmod)
testDispersion(mfecundmod)
testZeroInflation(mfecundmod)
mfecund_simResid <- simulateResiduals(fittedModel = mfecundmod)
plot(mfecund_simResid)

# Appendix Table S6
# contrasts among diets conditional on infection status and clone
mfecundmod_contrasts <- emmeans(mfecundmod, specs = pairwise ~ Diet | Infection + Clone, type = "response")
mfecundmod_contrasts

# Appendix Table S4: Fungus Total Fecundity
# contrasts among infection statuses conditional on clone and averaged across all diets
mfecundmod_contrasts2 <- emmeans(mfecundmod, specs = pairwise ~ Infection | Clone, type = "response")
mfecundmod_contrasts2



# test of microcystin effect on offspring production
m_data_fecund_toxintest <- filter(m_data_fecund, Diet == "M" | Diet == "M+")
mfecundmod3 <- glmmTMB(Total.Babies ~ microcystin.avg + Infection + Clone + (1|Block), family = nbinom1(), data = m_data_fecund_toxintest)
summary(mfecundmod3)
Anova(mfecundmod3)   # technically rank deficient to handle interactions, so removed all interactions
overdisp_fun(mfecundmod3)
AIC(mfecundmod3)
testDispersion(mfecundmod3)
testZeroInflation(mfecundmod3)
mfecund3_simResid <- simulateResiduals(fittedModel = mfecundmod3)
plot(mfecund3_simResid)


# Figure 3A: Fungus Total Fecundity
m_fecund_tot <- ggplot(m_data_fecund, aes(x = Infection, y = Total.Babies)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge2(width =0.8), show.legend = F) +
  geom_point(aes(color=Diet), size=1.5, position=position_jitterdodge(dodge.width=0.8, jitter.width = 0.25), alpha = 0.4, show.legend = F) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninfected", "Exposed", "Infected"), labels = c("Uninf", "Exp","Inf")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  ggtitle("Fungus Experiment") +
  #scale_y_log10(labels = scales::comma) +
  #ggtitle("Total Fecundity by infection status x diet x host clone for Metsch Expt") +
  labs(x = "Infection Status", y = "Total Fecundity") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=9, color="black"))
m_fecund_tot
ggsave("figures/MetschTotalFecundity_DietxClonexInfection.tiff", plot = m_fecund_tot, dpi = 600, width = 6.5, height = 4, units = "in", compression="lzw")



# Figure 5C: Fungus Total Fecundity Summary Fig
parasite_colors <- c("orange", "blue")
m_fecund_tot_summary <- ggplot(m_data_fecund, aes(x = Diet, y = Total.Babies)) +
  geom_boxplot(aes(fill=Parasite_Treatment),position=position_dodge2(width =0.8), show.legend = F) +
  geom_point(aes(color=Parasite_Treatment, shape = Parasite_Treatment), size=1.5, position=position_jitterdodge(dodge.width=0.8, jitter.width = 0.25), alpha = 0.6, show.legend = F) +
  facet_wrap(~Clone) +
  #scale_x_discrete(limits = c("Uninfected", "Exposed", "Infected"), labels = c("Uninf", "Exp","Inf")) +
  #scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = parasite_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  #scale_y_log10(labels = scales::comma) +
  ggtitle("Fungus Experiment") +
  labs(x = "Diet", y = "Total Fecundity") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=9, color="black"))
m_fecund_tot_summary
ggsave("figures/MetschTotalFecundity_DietxClonexParasiteTrmT.tiff", plot = m_fecund_tot_summary, dpi = 600, width = 6.5, height = 4, units = "in", compression="lzw")




# Figure of Total fecundity vs avg microcystin concentration by diet, clone, infection status (not in manuscript)
mc_diets <- c("#006837", "#1D91C0")
m_fecund_tot_micro <- ggplot(m_data_fecund_toxintest, aes(x = microcystin.avg+1, y = Total.Babies+1)) +
  geom_point(aes(color=Diet), size=2, alpha = 0.4) +
  geom_smooth(method = "lm", color = "black") +
  facet_grid(Clone~Infection) +
  scale_color_discrete(limits = c("M", "M+")) +
  scale_color_manual(values = mc_diets) +
  scale_y_log10(labels = scales::comma) +
  scale_x_log10() +
  ylim(0.9, 30) + 
  labs(x = "Average Microcystin Concentration (ug/L) log+1", y = "Total Fecundity (log + 1 scaled)") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=9, color="black"))
m_fecund_tot_micro
ggsave("figures/MetschTotalFecundity_MicrocystinxClonexInfection.tiff", plot = m_fecund_tot_micro, dpi = 600, width = 6.5, height = 4, units = "in", compression="lzw")





### Total Fecundity by treatment for Pasteuria experiment

p_data_fecund <- p_data %>% 
  filter(!is.na(Total.Babies)  & lifespan.days > 11) %>%  # first babies were observed on day 11 
  mutate(Parasite_Treatment = factor(if_else(Parasites == "Pasteuria", "Inoculated", "Uninoculated"), levels = c("Uninoculated", "Inoculated")))
hist(p_data$Total.Babies) 


# calculate N animals that were excluded 
p_data_fecund_exluded <- p_data %>% 
  filter(is.na(Total.Babies) | lifespan.days <= 11) %>% # first babies observed on day 11
  count(Diet, Parasites, Clone, Block)

p_total_n_excluded <- sum(p_data_fecund_exluded$n)


#pfecundmod <- glmmTMB(Total.Babies ~ Diet * Infection * Clone + (1|Block) + (1|Unique.code), family = poisson, data = p_data_fecund)  # overdispersed unless include unique code random effect
#pfecundmod <- glmmTMB(Total.Babies ~ Diet * Infection * Clone + (1|Block), family = nbinom1(), data = p_data_fecund)   # higher AIC than neg binomial
pfecundmod <- glmmTMB(Total.Babies ~ Diet * Infection * Clone + (1|Block), family = nbinom2(), data = p_data_fecund)  #lower AIC than poisson and quasipoisson, use this model
summary(pfecundmod)

# Table 1: Bacterium Total Fecundity
Anova(pfecundmod)        
overdisp_fun(pfecundmod)
AIC(pfecundmod)
qqnorm(resid(pfecundmod))
qqline(resid(pfecundmod))
testDispersion(pfecundmod)
testZeroInflation(pfecundmod)
pfecund_simResid <- simulateResiduals(fittedModel = pfecundmod)
plot(pfecund_simResid)

# Appendix Table S7
# contrasts among diets conditional on infection status and clone
pfecundmod_contrasts <- emmeans(pfecundmod, specs = pairwise ~ Diet | Infection + Clone, type = "response")
pfecundmod_contrasts

# Appendix Table S4: Bacterium Total Fecundity
# contrasts among infection statuses conditional on clone and averaged across all diets
pfecundmod_contrasts2 <- emmeans(pfecundmod, specs = pairwise ~ Infection | Clone, type = "response")
pfecundmod_contrasts2



# test of microcystin effect on offspring production
p_data_fecund_toxintest <- filter(p_data_fecund, Diet == "M" | Diet == "M+")
pfecundmod2 <- glmmTMB(Total.Babies ~ microcystin.avg * Infection * Clone + (1|Block), family = nbinom2(), data = p_data_fecund_toxintest)
summary(pfecundmod2)
Anova(pfecundmod2)          # microcystin concentration did not impact offspring production
overdisp_fun(pfecundmod2)
AIC(pfecundmod2)
qqnorm(resid(pfecundmod2))
qqline(resid(pfecundmod2))
testDispersion(pfecundmod2)
testZeroInflation(pfecundmod2)
pfecund2_simResid <- simulateResiduals(fittedModel = pfecundmod2)
plot(pfecund2_simResid)



# Figure 3B: Bacterium Total Fecundity
p_fecund_tot <- ggplot(p_data_fecund, aes(x = Infection, y = Total.Babies)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge2(width=0.8)) +
  geom_point(aes(color=Diet), size=1.5, position=position_jitterdodge(dodge.width=0.8, jitter.width = 0.25), alpha = 0.4) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninfected", "Exposed", "Infected"), labels = c("Uninf", "Exp","Inf")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  ggtitle("Bacterium Experiment") +
  labs(x = "Infection Status", y = "Total Fecundity") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
p_fecund_tot
ggsave("figures/PastTotalFecundity_DietxClonexInfection.tiff", plot = p_fecund_tot, dpi = 600, width = 6.5, height = 4, units = "in", compression="lzw")



# Figure 5D: Bacterium Total Fecundity Summary Fig
p_fecund_tot_summary <- ggplot(p_data_fecund, aes(x = Diet, y = Total.Babies)) +
  geom_boxplot(aes(fill=Parasite_Treatment),position=position_dodge2(width =0.8), show.legend = T) +
  geom_point(aes(color=Parasite_Treatment, shape = Parasite_Treatment), size=1.5, position=position_jitterdodge(dodge.width=0.8, jitter.width = 0.25), alpha = 0.6, show.legend = T) +
  facet_wrap(~Clone) +
  scale_fill_manual(values = parasite_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  ggtitle("Bacterium Experiment") +
  labs(x = "Diet", y = "Total Fecundity") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=9, color="black"))
p_fecund_tot_summary
ggsave("figures/PastTotalFecundity_DietxClonexParasiteTrmT.tiff", plot = p_fecund_tot_summary, dpi = 600, width = 6.5, height = 4, units = "in", compression="lzw")





# Figure of Total fecundity vs avg microcystin concentration by diet, clone, infection status (not in manuscript)
mc_diets <- c("#006837", "#1D91C0")
p_fecund_tot_micro <- ggplot(p_data_fecund_toxintest, aes(x = microcystin.avg+1, y = Total.Babies+1)) +
  geom_point(aes(color=Diet), size=2, alpha = 0.4) +
  geom_smooth(method = "lm", color = "black") +
  facet_grid(Clone~Infection) +
  scale_color_discrete(limits = c("M", "M+")) +
  scale_color_manual(values = mc_diets) +
  scale_y_log10(labels = scales::comma) +
  scale_x_log10() +
  #ggtitle("Total Fecundity by infection status x diet x host clone for Past Expt") +
  labs(x = "Average Microcystin Concentration (ug/L) log+1", y = "Total Fecundity (log + 1 scaled)") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=9, color="black"))
p_fecund_tot_micro
ggsave("figures/PastTotalFecundity_MicrocystinxClonexInfection.tiff", plot = p_fecund_tot_micro, dpi = 600, width = 6.5, height = 4, units = "in", compression="lzw")

  





## Survival ----------------------------------------------------------------------


### Metsch survival analyses

m_survival_df <- m_data %>%
  filter(lifespan.days > 11)

m_survival_df$Diet2 = recode_factor(m_survival_df$Diet, 'SM' = "50:50 SM", .default = levels(m_data$Diet))
m_survival_df$Diet2 <- factor(m_survival_df$Diet2, levels = c("S", "50:50 SM", "M", "M+"))
m_survival_df <- mutate(m_survival_df, Parasite_Treatment = factor(if_else(Parasites == "Metsch", "Inoculated", "Uninoculated"), levels = c("Uninoculated", "Inoculated")))

length(m_survival_df$SurvObj)
table(m_survival_df$Block)

# Calculate the number of animals censored at end of expt and early deaths
m_censored <- m_survival_df %>%
  filter(Censored == "censored") %>%
  count(Lifespan)


# Make a figure of survival curves for each clone and infection status
m_survival_df$SurvObj <- with(m_survival_df, Surv(lifespan.days, Status))

m_km.by.diet_infstatus <- survfit(SurvObj ~ Diet + Infection, data = m_survival_df, conf.type = "log-log")
m_km.by.diet_infstatus

summary(m_km.by.diet_infstatus)

# Figure 3E: Fungus expt survival plot by diet, clone, and infection status
m_survival_fig <- ggsurvplot_facet(m_km.by.diet_infstatus, data = m_survival_df, xlab = "Days", xlim = c(0,25), facet.by = c("Infection", "Clone"), palette = diet_colors, short.panel.labs = T) +
  #labs(color = "Diet") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"), legend.position = "none")
ggsave("figures/MetschSurvival.tiff", plot = m_survival_fig, dpi = 600, width = 4, height = 4.5, units = "in", compression="lzw")

# Same figure as above with p-values
m_survival_fig2 <- ggsurvplot_facet(m_km.by.diet_infstatus, data = m_survival_df, pval = TRUE, xlab = "Days", facet.by = c("Infection", "Clone"), palette = diet_colors, short.panel.labs = T) +
  labs(color = "Diet") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))



# Appendix Table S8: Fungus experiment hazard ratios
## Analysis with both clones
m_survival_df$SurvObj <- with(m_survival_df, Surv(lifespan.days, Status))
cox.m_diet_infstatus_clone <- coxph(SurvObj ~ Diet2 + Infection + Clone, data = m_survival_df)
summary(cox.m_diet_infstatus_clone)

#Likelihood ratio test
anova(cox.m_diet_infstatus_clone)
cox.m_diet_infstatus_clone$n

# Reported Wald Chisquare in Table S8
# Wald test for each parameter (null hypothesis is that the factor is not associated with survival outcome)
linearHypothesis(cox.m_diet_infstatus_clone, c("Diet250:50 SM", "Diet2M", "Diet2M+"))
linearHypothesis(cox.m_diet_infstatus_clone, c("Diet250:50 SM"))
linearHypothesis(cox.m_diet_infstatus_clone, c("Diet2M"))
linearHypothesis(cox.m_diet_infstatus_clone, c("Diet2M+"))
linearHypothesis(cox.m_diet_infstatus_clone, c("InfectionExposed", "InfectionInfected"))
linearHypothesis(cox.m_diet_infstatus_clone, c("InfectionExposed"))
linearHypothesis(cox.m_diet_infstatus_clone, c("InfectionInfected"))
linearHypothesis(cox.m_diet_infstatus_clone, c("CloneStandard "))


# forest plot for fungus by diet inf status and clone
metsch_cox <- ggforest(cox.m_diet_infstatus_clone, data = m_survival_df, cpositions = c(0.02, 0.15, 0.35))
metsch_cox
ggsave("figures/MetschCoxHazard_DietxClonexInfection.tiff", plot = metsch_cox, dpi = 600, width = 5.5, height = 4, units = "in", compression="lzw")


## Survival model with parasite treatment instead of infection status
m_km.by.diet_parasitetrmt <- survfit(SurvObj ~ Diet + Parasite_Treatment, data = m_survival_df, conf.type = "log-log")
m_km.by.diet_parasitetrmt

summary(m_km.by.diet_parasitetrmt)

# Fungus expt survival plot by diet, clone, and infection status (not included in manuscript)
m_survival_summary <- ggsurvplot_facet(m_km.by.diet_parasitetrmt, data = m_survival_df, xlab = "Days", xlim = c(0,25), facet.by = c("Diet", "Clone"), palette = parasite_colors, short.panel.labs = T) +
  #labs(color = "Diet") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"), legend.position = "none")
m_survival_summary
ggsave("figures/MetschSurvival_ParasiteTrmt.tiff", plot = m_survival_summary, dpi = 600, width = 4, height = 4.5, units = "in", compression="lzw")




# split data by clone to conduct separate survival analyses for each clone
m_data_std <- m_survival_df %>% filter(Clone == "Standard") %>% mutate(Treatment_Combo = paste(Diet, Infection, sep = "_"))
m_data_mid <- m_survival_df %>% filter(Clone == "Mid37") %>% mutate(Treatment_Combo = paste(Diet, Infection, sep = "_"))


# Mid37 survival in the Metsch experiment

m_data_mid$SurvObj <- with(m_data_mid, Surv(lifespan.days, Status))

# Kaplan-Meier estimator with "log-log" confidence interval
km.by.diet_parasite_mid <- survfit(SurvObj ~ Diet + Parasites, data = m_data_mid, conf.type = "log-log")
km.by.diet_infstatus_mid <- survfit(SurvObj ~ Diet + Infection, data = m_data_mid, conf.type = "log-log")

ggsurvplot(km.by.diet_parasite_mid, data = m_data_mid, pval = TRUE, xlab = "Days")
ggsurvplot(km.by.diet_infstatus_mid, data = m_data_mid, pval = TRUE, xlab = "Days")


# Cox proportional hazards model
cox.m_diet_parasite_mid <- coxph(SurvObj ~ Diet + Parasites, data = m_data_mid)
cox.m_diet_parasite_mid
cox.m_diet_infstatus_mid <- coxph(SurvObj ~ Diet + Infection, data = m_data_mid)
cox.m_diet_infstatus_mid

# Appendix Table S9: Fungus experiment Mid37 hazard ratios
cox.m_diet_infstatus_mid2 <- coxph(SurvObj ~ Diet2 + Infection, data = m_data_mid, singular.ok = T)   # Use this for the mid37 specific model
summary(cox.m_diet_infstatus_mid2)
Anova(cox.m_diet_infstatus_mid2)


# Reported Wald Chisquare in Table S9 for Fungus experiment Mid37
# Wald test for each parameter (null hypothesis is that the factor is not associated with survival outcome)
linearHypothesis(cox.m_diet_infstatus_mid2, c("Diet250:50 SM"))
linearHypothesis(cox.m_diet_infstatus_mid2, c("Diet2M"))
linearHypothesis(cox.m_diet_infstatus_mid2, c("Diet2M+"))
linearHypothesis(cox.m_diet_infstatus_mid2, c("InfectionExposed"))
linearHypothesis(cox.m_diet_infstatus_mid2, c("InfectionInfected"))

# forest plot
ggforest(cox.m_diet_parasite_mid, data = m_data_mid)
ggforest(cox.m_diet_infstatus_mid, data = m_data_mid)
ggforest(cox.m_diet_infstatus_mid2, data = m_data_mid)  # use this for the mid37 specific model




# Standard survival in the Metsch experiment 

m_data_std$SurvObj <- with(m_data_std, Surv(lifespan.days, Status))

# Kaplan-Meier estimator with "log-log" confidence interval
km.by.diet_std <- survfit(SurvObj ~ Diet, data = m_data_std, conf.type = "log-log")
km.by.diet_parasite_std <- survfit(SurvObj ~ Diet + Parasites, data = m_data_std, conf.type = "log-log")
km.by.diet_infstatus_std <- survfit(SurvObj ~ Diet + Infection, data = m_data_std, conf.type = "log-log")

#pairwise survdiff
res <- pairwise_survdiff(Surv(lifespan.days, Status) ~ Diet + Infection,
                         data = m_data_std)
res

ggsurvplot(km.by.diet_std, data = m_data_std, pval = TRUE, xlab = "Days")
ggsurvplot(km.by.diet_parasite_std, data = m_data_std, pval = TRUE, xlab = "Days")
ggsurvplot(km.by.diet_infstatus_std, data = m_data_std, pval = TRUE, xlab = "Days")




# Cox proportional hazards model
cox.m_diet_parasite_std <- coxph(SurvObj ~ Diet + Parasites, data = m_data_std)
cox.m_diet_parasite_std
cox.m_diet_infstatus_std <- coxph(SurvObj ~ Diet + Infection, data = m_data_std)
summary(cox.m_diet_infstatus_std)

# Appendix Table S9: Fungus experiment Standard hazard ratios
cox.m_diet_infstatus_std2 <- coxph(SurvObj ~ Diet2 + Infection, data = m_data_std)   # Use this for the standard specific model
summary(cox.m_diet_infstatus_std2)

# Reported Wald Chisquare in Table S9 for Fungus experiment Standard
# Wald test for each parameter (null hypothesis is that the factor is not associated with survival outcome)
linearHypothesis(cox.m_diet_infstatus_std2, c("Diet250:50 SM"))
linearHypothesis(cox.m_diet_infstatus_std2, c("Diet2M"))
linearHypothesis(cox.m_diet_infstatus_std2, c("Diet2M+"))
linearHypothesis(cox.m_diet_infstatus_std2, c("InfectionExposed"))
linearHypothesis(cox.m_diet_infstatus_std2, c("InfectionInfected"))

# forest plot
ggforest(cox.m_diet_parasite_std, data = m_data_std)
ggforest(cox.m_diet_infstatus_std, data = m_data_std)
ggforest(cox.m_diet_infstatus_std2, data = m_data_std)  # Use this for the standard specific model










### Pasteuria survival analyses

p_survival_df <- p_data %>%
  filter(lifespan.days > 11)

p_survival_df$Diet2 = recode_factor(p_survival_df$Diet, 'SM' = "50:50 SM", .default = levels(p_survival_df$Diet2))
p_survival_df$Diet2 <- factor(p_survival_df$Diet2, levels = c("S", "50:50 SM", "M", "M+"))
p_survival_df <- mutate(p_survival_df, Parasite_Treatment = factor(if_else(Parasites == "Pasteuria", "Inoculated", "Uninoculated"), levels = c("Uninoculated", "Inoculated")))

table(p_survival_df$Parasites)

# Calculate the number of animals censored at end of expt and early deaths
p_censored <- p_survival_df %>%
  filter(Censored == "censored") %>%
  count(Lifespan)



# Make a figure of survial curves for each clone and infection status
p_survival_df$SurvObj <- with(p_survival_df, Surv(lifespan.days, Status))

p_km.by.diet_infstatus <- survfit(SurvObj ~ Diet + Infection + Clone, data = p_survival_df, conf.type = "log-log")
p_km.by.diet_infstatus

# Figure 3F: Bacterium expt survival plot by diet, clone, and infection status
p_survival_fig <- ggsurvplot_facet(p_km.by.diet_infstatus, data = p_survival_df, xlab = "Days", xlim = c(0,42), facet.by = c("Infection", "Clone"), palette = diet_colors, short.panel.labs = T) +
  labs(color = "Diet") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
p_survival_fig
ggsave("figures/PasteuriaSurvival.tiff", plot = p_survival_fig, dpi = 600, width = 4, height = 4.5, units = "in", compression="lzw")

# same figure as above with p-values
p_survival_fig2 <- ggsurvplot_facet(p_km.by.diet_infstatus, data = p_survival_df, pval = TRUE, xlab = "Days", facet.by = c("Infection", "Clone"), palette = diet_colors) +
  labs(color = "Diet") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))


# Appendix Table S8: Bacterium experiment hazard ratios
## Analysis with both clones
cox.diet_infstatus_clone <- coxph(SurvObj ~ Diet2 + Infection + Clone, data = p_survival_df)
summary(cox.diet_infstatus_clone)


# Wald test for each parameter (null hypothesis is that the factor is not associated with survival outcome)
linearHypothesis(cox.diet_infstatus_clone, c("Diet250:50 SM", "Diet2M", "Diet2M+"))
linearHypothesis(cox.diet_infstatus_clone, c("Diet250:50 SM"))
linearHypothesis(cox.diet_infstatus_clone, c("Diet2M"))
linearHypothesis(cox.diet_infstatus_clone, c("Diet2M+"))
linearHypothesis(cox.diet_infstatus_clone, c("InfectionExposed", "InfectionInfected"))
linearHypothesis(cox.diet_infstatus_clone, c("InfectionExposed"))
linearHypothesis(cox.diet_infstatus_clone, c("InfectionInfected"))
linearHypothesis(cox.diet_infstatus_clone, c("CloneStandard "))

# Likelihood ratio test
anova(cox.diet_infstatus_clone)

# forest plot for bacterium diet inf status and clone
past_cox <- ggforest(cox.diet_infstatus_clone, data = p_survival_df, cpositions = c(0.02, 0.15, 0.35))
past_cox
ggsave("figures/PastCoxHazard_DietxClonexInfection.tiff", plot = past_cox, dpi = 600, width = 6, height = 4, units = "in", compression="lzw")

# NOTE: the forest plot has a problem plotting the reference if it starts with the same letter as another factor in the same variable (e.g. S and SM), but the analysis runs correctly. I have fixed it by changing the name of the SM diet to 50:50 SM.




## survival model for diet by parasite treatment
p_km.by.diet_parasitetrmt <- survfit(SurvObj ~ Diet + Parasite_Treatment, data = p_survival_df, conf.type = "log-log")
p_km.by.diet_parasitetrmt

summary(m_km.by.diet_parasitetrmt)

# Bacterium expt survival plot by diet, clone, and infection status (figure not shown in manuscript)
p_survival_summary <- ggsurvplot_facet(p_km.by.diet_parasitetrmt, data = p_survival_df, xlab = "Days", xlim = c(0,40), facet.by = c("Diet", "Clone"), palette = parasite_colors, short.panel.labs = T) +
  #labs(color = "Diet") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"), legend.position = "none")
p_survival_summary
ggsave("figures/PastSurvival_ParasiteTrmt.tiff", plot = p_survival_summary, dpi = 600, width = 4, height = 4.5, units = "in", compression="lzw")







# split data by clone to run separate survival analyses for each clone
p_data_std <- p_survival_df %>% filter(Clone == "Standard") %>% mutate(Treatment_Combo = paste(Diet, Infection, sep = "_"))
p_data_mid <- p_survival_df %>% filter(Clone == "Mid37") %>% mutate(Treatment_Combo = paste(Diet, Infection, sep = "_"))



# Standard survival in the Pasteuria experiment
  
# Kaplan-Meier models estimate the survival probability and the long-rank test. 


p_data_std$SurvObj <- with(p_data_std, Surv(lifespan.days, Status))


# Kaplan-Meier estimator with "log-log" confidence interval
km.by.diet_parasite_std <- survfit(SurvObj ~ Diet + Parasites, data = p_data_std, conf.type = "log-log")
km.by.diet_parasite_std
km.by.diet_infstatus_std <- survfit(SurvObj ~ Diet + Infection, data = p_data_std, conf.type = "log-log")
km.by.diet_infstatus_std
km.by.diet_infstatus_std2 <- survfit(SurvObj ~ Diet2 + Infection, data = p_data_std, conf.type = "log-log")
km.by.diet_infstatus_std2

ggsurvplot(km.by.diet_parasite_std, data = p_data_std, pval = TRUE, xlab = "Days")
ggsurvplot(km.by.diet_infstatus_std, data = p_data_std, pval = TRUE, xlab = "Days")
ggsurvplot(km.by.diet_infstatus_std2, data = p_data_std, pval = TRUE, xlab = "Days")


# Cox proportional hazards model calculates the risk of death and respective hazard ratios (HR). HR < 1 indicates an increased risk of death, and HR < 1 indicates a decreased risk of death.



# Cox proportional hazards model
cox.diet_parasite_std <- coxph(SurvObj ~ Diet + Parasites, data = p_data_std)
cox.diet_parasite_std
cox.diet_infstatus_std <- coxph(SurvObj ~ Diet + Infection, data = p_data_std)
summary(cox.diet_infstatus_std)

# Appendix Table S9: Bacterium experiment Standard hazard ratios
cox.diet_infstatus_std2 <- coxph(SurvObj ~ Diet2 + Infection, data = p_data_std)  ## Use this one for the Standard specific results
summary(cox.diet_infstatus_std2)

# Reported Wald Chisquare in Table S9 for Bacterium experiment Standard
# Wald test for each parameter (null hypothesis is that the factor is not associated with survival outcome)
linearHypothesis(cox.diet_infstatus_std2, c("Diet250:50 SM"))
linearHypothesis(cox.diet_infstatus_std2, c("Diet2M"))
linearHypothesis(cox.diet_infstatus_std2, c("Diet2M+"))
linearHypothesis(cox.diet_infstatus_std2, c("InfectionExposed"))
linearHypothesis(cox.diet_infstatus_std2, c("InfectionInfected"))

# forest plot
ggforest(cox.diet_parasite_std, data = p_data_std)
ggforest(cox.diet_infstatus_std, data = p_data_std)
ggforest(cox.diet_infstatus_std2, data = p_data_std)  # Use this one for the Standard specific results
ggforest(cox.diet_treatments_std, data = p_data_std)  # doesn't plot result




# Mid37 survival in the Pasteuria experiment
  
p_data_mid$SurvObj <- with(p_data_mid, Surv(lifespan.days, Status))

# Kaplan-Meier estimator with "log-log" confidence interval
km.by.diet_parasite_mid <- survfit(SurvObj ~ Diet + Parasites, data = p_data_mid, conf.type = "log-log")
km.by.diet_infstatus_mid <- survfit(SurvObj ~ Diet + Infection, data = p_data_mid, conf.type = "log-log")

ggsurvplot(km.by.diet_parasite_mid, data = p_data_mid, pval = TRUE, xlab = "Days")
ggsurvplot(km.by.diet_infstatus_mid, data = p_data_mid, pval = TRUE, xlab = "Days")



# Cox proportional hazards model
cox.diet_parasite_mid <- coxph(SurvObj ~ Diet + Parasites, data = p_data_mid)
cox.diet_parasite_mid
cox.diet_infstatus_mid <- coxph(SurvObj ~ Diet + Infection, data = p_data_mid)
cox.diet_infstatus_mid

# Appendix Table S9: Bacterium experiment Mid37 hazard ratios
cox.diet_infstatus_mid2 <- coxph(SurvObj ~ Diet2 + Infection, data = p_data_mid)  # Use this one for the Mid37 specific results
summary(cox.diet_infstatus_mid2)

# Reported Wald Chisquare in Table S9 for Bacterium experiment Mid37
# Wald test for each parameter (null hypothesis is that the factor is not associated with survival outcome)
linearHypothesis(cox.diet_infstatus_mid2, c("Diet250:50 SM"))
linearHypothesis(cox.diet_infstatus_mid2, c("Diet2M"))
linearHypothesis(cox.diet_infstatus_mid2, c("Diet2M+"))
linearHypothesis(cox.diet_infstatus_mid2, c("InfectionExposed"))
linearHypothesis(cox.diet_infstatus_mid2, c("InfectionInfected"))

# forest plot
ggforest(cox.diet_parasite_mid, data = p_data_mid)
ggforest(cox.diet_infstatus_mid, data = p_data_mid)
ggforest(cox.diet_infstatus_mid2, data = p_data_mid)  # Use this one for the Mid37 specific results




# Figure 3
Figure3 <- ggarrange(m_fecund_tot, p_fecund_tot, m_survival_fig, p_survival_fig, labels = "auto",
                     ncol = 2, nrow = 2, heights = c(3,3), legend = "right", common.legend = T)
Figure3
ggsave("figures/manuscript/Fig3_Fecundity_Survival.tiff", plot = Figure3, dpi = 600, width = 8, height = 6, units = "in", compression="lzw")



# Figure 5: Summary figure
Figure5 <- ggarrange(m_fecund_tot_summary, p_fecund_tot_summary, labels = "auto",
                     ncol = 2, nrow = 1, heights = c(3.2,3,3), legend = "right", common.legend = T)
Figure5
ggsave("figures/manuscript/Fig5_Summary_Fecundity.tiff", plot = Figure5, dpi = 600, width = 8.5, height = 4, units = "in", compression="lzw")




### Lifespan [ANALYESES NOT INCLUDED IN MANUSCRIUPT DUE TO POOR MODEL FIT]

# Metsch lifespan
m_data_lifespan <- m_data %>% filter(Include_Lifespan == "Y")
hist(m_data_lifespan$lifespan.days)


mlifemod <- glmmTMB(lifespan.days ~ Diet + Infection + Clone + Diet:Infection + Diet:Clone + Infection:Clone + (1|Block) + (1|Unique.code), family = poisson(), data = m_data_lifespan)
summary(mlifemod) # model did not have enough power to include the three-way interaction
Anova(mlifemod)
AIC(mlifemod)
testDispersion(mlifemod)
testZeroInflation(mlifemod)
mlife_simResid <- simulateResiduals(fittedModel = mlifemod)
plot(mlife_simResid)  # this model has some problems. Not a good fit. Use the survival analyses instead.


mlifemod_contrasts <- emmeans(mlifemod, specs = pairwise ~ Diet | Infection + Clone, type = "response")
mlifemod_contrasts


mlifemod_contrasts2 <- emmeans(mlifemod, specs = pairwise ~ Infection | Clone, type = "response")
mlifemod_contrasts2


mlifemod_contrasts3 <- emmeans(mlifemod, specs = pairwise ~ Diet, type = "response")
mlifemod_contrasts3




m_lifespan <- ggplot(m_data_lifespan, aes(x = Infection, y = Lifespan)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width =0.9)) +
  geom_point(aes(color=Diet), size=2, position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.2), alpha = 0.4) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninfected", "Exposed", "Infected")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  labs(x = "Parasite Exposure", y = "Lifespan") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=9, color="black"))
m_lifespan

ggsave("figures/MetschLifespan_DietxClonexInfection.tiff", plot = m_lifespan, dpi = 600, width = 6.5, height = 4, units = "in", compression="lzw")



# does fecundity data explain differences in lifespan?
m_fecund_lifespan <- ggplot(m_data_lifespan, aes(x = Total.Babies+1, y = Lifespan)) +
  geom_point(aes(color=Diet), size=2, alpha = 0.6) +
  geom_smooth(method = "lm", color = "black") +
  facet_wrap(Clone~Infection) +
  scale_color_manual(values = diet_colors) +
  scale_x_log10() +
  labs(x = "Total Fecundity", y = "Lifespan") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=9, color="black"))
m_fecund_lifespan

# same plot as above but only for S  and SM treatments
m_data_lifespan_S <- filter(m_data_lifespan, Diet == "S" | Diet == "SM")
m_fecund_lifespan_S <- ggplot(m_data_lifespan_S, aes(x = Total.Babies+1, y = Lifespan)) +
  geom_point(aes(color=Diet), size=2, alpha = 0.6) +
  geom_smooth(method = "lm", color = "black") +
  facet_wrap(Clone~Infection) +
  scale_color_manual(values = diet_colors) +
  scale_x_log10() +
  labs(x = "Total Fecundity", y = "Lifespan") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=9, color="black"))
m_fecund_lifespan_S








# Pasteuria lifespan

p_data_lifespan <- p_data %>% filter(Include_Lifespan == "Y")
hist(p_data_lifespan$lifespan.days)


plifemod <- glmmTMB(lifespan.days ~ Diet * Infection * Clone +(1|Block), family = nbinom1(), data = p_data_lifespan)
summary(plifemod)
Anova(plifemod)   # use this for the manuscript
AIC(plifemod)
testDispersion(plifemod)
testZeroInflation(plifemod)
plife_simResid <- simulateResiduals(fittedModel = plifemod)
plot(plife_simResid) # this model has some problems. Not a good fit. Use the survival analyses instead.


plifemod_contrasts <- emmeans(plifemod, specs = pairwise ~ Diet | Infection + Clone, type = "response")
plifemod_contrasts


plifemod_contrasts2 <- emmeans(plifemod, specs = pairwise ~ Infection | Clone, type = "response")
plifemod_contrasts2


plifemod_contrasts3 <- emmeans(plifemod, specs = pairwise ~ Diet, type = "response")
plifemod_contrasts3



p_lifespan <- ggplot(p_data_lifespan, aes(x = Infection, y = Lifespan)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width =0.9)) +
  geom_point(aes(color=Diet), size=1.5, position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.2), alpha = 0.4) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninfected", "Exposed", "Infected")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  ggtitle("Bacteria Experiment") +
  labs(x = "Parasite Exposure", y = "Lifespan") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
p_lifespan

ggsave("figures/PastLifespan_DietxClonexInfection.tiff", plot = p_lifespan, dpi = 600, width = 6.5, height = 4, units = "in", compression="lzw")



# does fecundity data explain differences in lifespan?
p_fecund_lifespan <- ggplot(p_data_lifespan, aes(x = Total.Babies+1, y = Lifespan)) +
  geom_point(aes(color=Diet), size=2, alpha = 0.6) +
  geom_smooth(method = "lm", color = "black") +
  facet_wrap(Clone~Infection) +
  scale_color_manual(values = diet_colors) +
  scale_x_log10() +
  labs(x = "Total Fecundity", y = "Lifespan") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=9, color="black"))
p_fecund_lifespan

# same plot as above but only for S treatment
p_data_lifespan_S <- filter(p_data_lifespan, Diet == "S")
p_fecund_lifespan_S <- ggplot(p_data_lifespan_S, aes(x = Total.Babies+1, y = Lifespan)) +
  geom_point(aes(color=Diet), size=2, alpha = 0.6) +
  geom_smooth(method = "lm", color = "black") +
  facet_wrap(Clone~Infection) +
  scale_color_manual(values = diet_colors) +
  scale_x_log10() +
  labs(x = "Total Fecundity", y = "Lifespan") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=9, color="black"))
p_fecund_lifespan_S







## Body size during exposure and Growth Rate ----------------------------------------------------

### Body size during exposure ------------------------------------------------------


## Body size differences during Week 1 (parasite exposure) for Fungus expt
m_data_size <- m_data %>% filter(!is.na(Length_Week1))
hist(m_data_size$Length_Week1) # wow normal data!!

msizemod <- lmer(log(Length_Week1) ~ Diet + Infection + Clone + Diet:Infection + Diet:Clone  + Infection:Clone +(1|Block), data = m_data_size)
# Model had problems when all interactions and random effect included, but none of the interactions are significant
# so I removed the three-way interaction and kept the random effect
summary(msizemod)

# Appendix Table S2
Anova(msizemod)   # Use this model for the manuscript
AIC(msizemod)
plot(msizemod)
qqnorm(resid(msizemod))        # looks great
qqline(resid(msizemod))
shapiro.test(resid(msizemod)) 

msizemod_contrasts <- emmeans(msizemod, specs = pairwise ~ Diet | Infection + Clone, type = "response")
msizemod_contrasts

# Appendix Table S4
msizemod_contrasts2 <- emmeans(msizemod, specs = pairwise ~ Infection | Clone, type = "response")
msizemod_contrasts2

# difference in size between clones
msizemod_contrasts3 <- emmeans(msizemod, specs = pairwise ~ Clone, type = "response")
msizemod_contrasts3



# Figure S4A: Fungus Expt Body Size during inoculation period
m_size <- ggplot(m_data_size, aes(x = Infection, y = Length_Week1)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width =0.9), show.legend = F) +
  geom_point(aes(color=Diet), size=1.5, position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.2), alpha = 0.4, show.legend = F) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninfected", "Exposed", "Infected"), labels = c("Uninf", "Exp", "Inf")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  ggtitle("Fungus Experiment") +
  labs(x = "Infection Status", y = bquote("Body Size During Inoculation Period ("* mu* "m)")) +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
m_size
ggsave("figures/MetschBodySize_DietxClonexInfection.tiff", plot = m_size, dpi = 600, width = 6.5, height = 4, units = "in", compression="lzw")



## manipulated data to get body size over time by week
m_longdata_size_sum <- m_longdata %>%
  filter(!is.na(WeeklyLength)) %>%
  group_by(Diet, Infection, Clone, Week) %>%
  dplyr::summarize(N  = length(WeeklyLength),
                   size.mean = mean(WeeklyLength, na.rm = T),
                   size.sd   = sd(WeeklyLength, na.rm = T),
                   size.se   = size.sd / sqrt(N)) %>%
  mutate(Treatment = paste0(Diet, "_", Infection, "_", Clone))
m_longdata_size_sum$Week <- as.factor(m_longdata_size_sum$Week)
m_longdata_size_sum$Diet <- factor(m_longdata_size_sum$Diet, levels = c("S", "SM", "M", "M+"))


# Figure 4A: Fungus expt Body size by week
plot_byWeek_size_metsch <- ggplot(data = m_longdata_size_sum, aes(x = Week, y = size.mean, group = Treatment, shape = Diet)) +
  geom_errorbar(aes(x=Week, ymin=size.mean-size.se, ymax=size.mean+size.se,color=Diet), width=0.3, show.legend = F) + 
  geom_line(aes(color=Diet, group = Treatment), linewidth=1, show.legend = F) +
  geom_point(aes(color=Diet), size=2.5, alpha = 0.8, show.legend = F) +
  facet_grid(Clone~Infection) +
  ylim(600, 2000) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  annotate("rect", xmin = 0.9, xmax = 1.1, ymin = 600, ymax = 2000,
           alpha = .2,fill = "gray") +
  ggtitle("Fungus Experiment") +
  labs(x = "Week of Experiment", y = bquote("Body Size ("* mu* "m)")) +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
plot_byWeek_size_metsch
ggsave(here("figures/MetschBodySize_DietxCloneXInfection_byWeek.tiff"), plot = plot_byWeek_size_metsch, dpi = 300, width = 7, height = 6, units = "in", compression="lzw")






## Body size differences during Week 1 (parasite exposure) for Bacterium expt


p_data_size <- p_data %>% filter(!is.na(Length_Week1), Length_Week1 > 765)  # remove 3 low outliers that made the residual distribution not normal and were much lower than all other animals
hist(p_data_size$Length_Week1) # wow normal data!!

library(EnvStats)
outlier_test_bodysize <- rosnerTest(p_data_size$Length_Week1, k = 10)
outlier_test_bodysize$all.stats


psizemod <- lmer(log(Length_Week1) ~ Diet * Infection * Clone +(1|Block), data = p_data_size)
summary(psizemod)

# Appendix Table S2
Anova(psizemod)   # Use this model for the manuscript
AIC(psizemod)
plot(psizemod)
qqnorm(resid(psizemod))   # looks ok
qqline(resid(psizemod))
shapiro.test(resid(psizemod)) # not sig, but 0.079


psizemod_contrasts <- emmeans(psizemod, specs = pairwise ~ Diet | Infection + Clone, type = "response")
psizemod_contrasts

# Appendix Table S4
psizemod_contrasts2 <- emmeans(psizemod, specs = pairwise ~ Infection | Clone, type = "response")
psizemod_contrasts2

# difference in size between clones
psizemod_contrasts3 <- emmeans(psizemod, specs = pairwise ~ Clone, type = "response")
psizemod_contrasts3

psizemod_contrasts4 <- emmeans(psizemod, specs = pairwise ~ Diet, type = "response")
psizemod_contrasts4


# Figure S4B: Bacterium Expt Body Size during inoculation period
p_size <- ggplot(p_data_size, aes(x = Infection, y = Length_Week1)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width =0.9)) +
  geom_point(aes(color=Diet), size=1.5, position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.2), alpha = 0.4) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninfected", "Exposed", "Infected"), labels = c("Uninf", "Exp", "Inf")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  ggtitle("Bacterium Experiment") +
  labs(x = "Infection Status", y = bquote("Body Size During Inoculation Period ("* mu* "m)")) +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
p_size
ggsave("figures/PastBodySize_DietxClonexInfection.tiff", plot = p_size, dpi = 600, width = 6.5, height = 4, units = "in", compression="lzw")


## manipulated data to get body size over time by week
p_longdata_size_sum <- p_longdata %>%
  filter(!is.na(WeeklyLength)) %>%
  group_by(Diet, Infection, Clone, Week) %>%
  dplyr::summarize(N  = length(WeeklyLength),
                   size.mean = mean(WeeklyLength, na.rm = T),
                   size.sd   = sd(WeeklyLength, na.rm = T),
                   size.se   = size.sd / sqrt(N)) %>%
  mutate(Treatment = paste0(Diet, "_", Infection, "_", Clone))
p_longdata_size_sum$Week <- as.factor(p_longdata_size_sum$Week)
p_longdata_size_sum$Diet <- factor(p_longdata_size_sum$Diet, levels = c("S", "SM", "M", "M+"))


# Figure 4B: Bacterium expt Body size by week
plot_byWeek_size_past <- ggplot(data = p_longdata_size_sum, aes(x = Week, y = size.mean, group = Treatment, shape = Diet)) +
  geom_errorbar(aes(x=Week, ymin=size.mean-size.se, ymax=size.mean+size.se,color=Diet), width=0.3) + 
  geom_line(aes(color=Diet, group = Treatment), linewidth=1) +
  geom_point(aes(color=Diet), size=2.5, alpha = 0.8) +
  facet_grid(Clone~Infection) +
  ylim(600, 2000) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  annotate("rect", xmin = 0.9, xmax = 1.1, ymin = 600, ymax = 2000,
           alpha = .2,fill = "gray") +
  ggtitle("Bacterium Experiment") +
  labs(x = "Week of Experiment", y = bquote("Body Size ("* mu* "m)")) +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
plot_byWeek_size_past
ggsave(here("figures/PastBodySize_DietxCloneXInfection_byWeek.tiff"), plot = plot_byWeek_size_past, dpi = 300, width = 7, height = 6, units = "in", compression="lzw")





# Figure 4 
Figure4 <- ggarrange(plot_byWeek_size_metsch, plot_byWeek_size_past, labels = "auto", ncol = 2, nrow = 1, widths = c(3.3,4))
Figure4
ggsave("figures/manuscript/Fig4_BodySize_byWeek.tiff", plot = Figure4, dpi = 600, width = 8.5, height = 4, units = "in", compression="lzw")






### Growth rate -----------------------------------------------------------------

### Growth rate between week 1 and 2 for Fungus experiment

m_data_growth <- m_data %>% filter(!is.na(GrowthRate))
hist(m_data_growth$GrowthRate) # wow normal data!! 

mgrowthmod <- lmer(GrowthRate ~ Diet + Infection + Clone + Diet:Infection + Diet:Clone  + Infection:Clone +(1|Block), data = m_data_growth)
# Model had problems when all interactions and random effect included, but none of the interactions are significant
# so I removed the three-way interaction and kept the random effect
summary(mgrowthmod)

# Appendix Table S2
Anova(mgrowthmod)   # Use this model for the manuscript
AIC(mgrowthmod)
plot(mgrowthmod)
qqnorm(resid(mgrowthmod))
qqline(resid(mgrowthmod))
shapiro.test(resid(mgrowthmod)) # not sig, but 0.074

mgrowthmod_contrasts <- emmeans(mgrowthmod, specs = pairwise ~ Diet | Infection + Clone, type = "response")
mgrowthmod_contrasts

# Appendix Table S4: Fungus Expt body size growth rate
mgrowthmod_contrasts2 <- emmeans(mgrowthmod, specs = pairwise ~ Infection | Clone, type = "response")
mgrowthmod_contrasts2

# check if there is a sig interaction for infection x clone when only including uninfected and exposed animals
m_data_growth_noInf <- filter(m_data_growth, Infection != "Infected")
mgrowthmod2 <- lmer(GrowthRate ~ Diet + Infection + Clone + Diet:Infection + Diet:Clone  + Infection:Clone +(1|Block), data = m_data_growth_noInf)
summary(mgrowthmod2)
Anova(mgrowthmod2)
AIC(mgrowthmod2)
plot(mgrowthmod2)
qqnorm(resid(mgrowthmod2))
qqline(resid(mgrowthmod2))
shapiro.test(resid(mgrowthmod2))  # significant, not the best fit

mgrowthmod_contrasts3 <- emmeans(mgrowthmod2, specs = pairwise ~ Infection | Clone, type = "response")
mgrowthmod_contrasts3



# Figure S4C: Growth rate between week 1 and 2 for Fungus experiment
m_growth <- ggplot(m_data_growth, aes(x = Infection, y = GrowthRate)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width =0.9), show.legend = F) +
  geom_point(aes(color=Diet), size=1.5, position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.2), alpha = 0.4, show.legend = F) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninfected", "Exposed", "Infected"), labels = c("Uninf", "Exp", "Inf")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  #ggtitle("Fungus Experiment") +
  labs(x = "Infection Status", y = bquote("Body Size Growth Rate ("* mu* "m" ~day^-1*")")) +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
m_growth
ggsave("figures/MetschGrowthRate_DietxClonexInfection.tiff", plot = m_growth, dpi = 600, width = 6.5, height = 4, units = "in", compression="lzw")




### Growth rate between week 1 and 2 for Bacterium experiment

p_data_growth <- p_data %>% filter(!is.na(GrowthRate))
hist(p_data_growth$GrowthRate)  # wow normal data!! 

pgrowthmod <- lmer(GrowthRate ~ Diet * Infection * Clone + (1|Block), data = p_data_growth)
summary(pgrowthmod)   # Block is singular

# Appendix Table S2
Anova(pgrowthmod)
AIC(pgrowthmod)
plot(pgrowthmod)
qqnorm(resid(pgrowthmod))   # looks good
qqline(resid(pgrowthmod))
shapiro.test(resid(pgrowthmod)) # not sig


pgrowthmod_contrasts <- emmeans(pgrowthmod, specs = pairwise ~ Diet | Infection + Clone, type = "response")
pgrowthmod_contrasts

# Appendix Table S4
pgrowthmod_contrasts2 <- emmeans(pgrowthmod, specs = pairwise ~ Infection | Clone, type = "response")
pgrowthmod_contrasts2

pgrowthmod_contrasts3 <- emmeans(pgrowthmod, specs = pairwise ~ Clone, type = "response")
pgrowthmod_contrasts3



# Figure S4C: Growth rate between week 1 and 2 for Bacterium experiment
p_growth <- ggplot(p_data_growth, aes(x = Infection, y = GrowthRate)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width =0.9)) +
  geom_point(aes(color=Diet), size=1.5, position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.2), alpha = 0.4) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninfected", "Exposed", "Infected"), labels = c("Uninf", "Exp", "Inf")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  #ggtitle("Bacteria Experiment") +
  labs(x = "Infection Status", y = bquote("Body Size Growth Rate ("* mu* "m" ~day^-1*")")) +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
p_growth
ggsave("figures/PastGrowthRate_DietxClonexInfection.tiff", plot = p_growth, dpi = 600, width = 6.5, height = 4, units = "in", compression="lzw")





# Figure S4
FigureS4 <- ggarrange(m_size, p_size, m_growth, p_growth, labels = "auto", ncol = 2, nrow = 2, widths = c(3.3,4))
FigureS4
ggsave("figures/manuscript/FigS4_GrowthRate_BodySize.tiff", plot = FigureS4, dpi = 600, width = 8.5, height = 8, units = "in", compression="lzw")




# Metabolic Rate -------------------------------------------------------------------


## Metsch metabolic rate analyses
m_data_resp <- m_data %>% dplyr::filter(!is.na(metabolic.rate.Week1), metabolic.rate.Week1 < 0.006)  # remove two outliers that cause extra right skew in the data, removing improves qqplot fit and resid plot
hist(as.numeric(m_data_resp$metabolic.rate.Week1)) # More right skewed than I would like
range(m_data_resp$metabolic.rate.Week1)
library(EnvStats)
outlier_test_metabolic <- rosnerTest(m_data_resp$metabolic.rate.Week, k = 4)
outlier_test_metabolic$all.stats


# Metabolic rate during week 1 (exposure period) only
mrespmod <- lmer(metabolic.rate.Week1 ~ Diet + Infection + Clone + Diet:Infection + Diet:Clone + Infection:Clone +(1|Block), data = m_data_resp)
summary(mrespmod)

# Appendix Table S2
Anova(mrespmod)     # use this for the manuscript
AIC(mrespmod)
plot(mrespmod)
qqnorm(resid(mrespmod)) # looks pretty good
qqline(resid(mrespmod))
shapiro.test(resid(mrespmod))   # significant, p = 0.043, But doesn't look too severe from the qqnorm plots (see https://stats.stackexchange.com/questions/32957/testing-normality-assumptions-for-linear-mixed-models-and-mixed-repeated-glm-a)

mrespmod_contrasts <- emmeans(mrespmod, specs = pairwise ~ Diet | Infection + Clone, type = "response")
mrespmod_contrasts

# Appendix Table S4
mrespmod_contrasts2 <- emmeans(mrespmod, specs = pairwise ~ Infection | Clone, type = "response")
mrespmod_contrasts2

mrespmod_contrasts3 <- emmeans(mrespmod, specs = pairwise ~ Clone, type = "response")
mrespmod_contrasts3

mrespmod_contrasts4 <- emmeans(mrespmod, specs = pairwise ~ Diet, type = "response")
mrespmod_contrasts4




# Figure S5A: Metabolic rate in week 1 for Fungus experiment
m_resp_fig <- ggplot(m_data_resp, aes(x = Infection, y = metabolic.rate.Week1)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width =0.9), show.legend = F) +
  geom_point(aes(color=Diet), size=1.5, position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.2), alpha = 0.4, show.legend = F) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninfected", "Exposed", "Infected"), labels = c("Uninf", "Exp", "Inf")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  labs(x = "Infection Status", y = bquote("Metabolic Rate (J/hr)")) +
  ggtitle("Fungus Experiment") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
m_resp_fig
ggsave("figures/MetschMetabolicRate_DietxClonexInfection_Week1.tiff", plot = m_resp_fig, dpi = 600, width = 5.25, height = 4, units = "in", compression="lzw")


# Metabolic rate over whole experiment (time series)
# Figure only, model did not fit well

# Use long data set for time series data of metabolic rate
m_longdata_resp <- m_longdata %>% dplyr::filter(!is.na(metabolic.rate))
hist(as.numeric(m_longdata_resp$metabolic.rate)) # more normal than above, though a bit right skewed
range(m_longdata_resp$metabolic.rate)


m_longdata_resp$Week <- as.factor(m_longdata_resp$Week)
m_longdata_resp$Diet <- factor(m_longdata_resp$Diet, levels = c("S", "SM", "M", "M+"))
m_longdata_resp_sum <- m_longdata_resp %>%
  ungroup() %>% dplyr::group_by(Diet, Infection, Clone, Week) %>%
  dplyr::summarize(N  = length(metabolic.rate),
                   resp.mean = mean(metabolic.rate, na.rm = T),
                   resp.sd   = sd(metabolic.rate, na.rm = T),
                   resp.se   = resp.sd / sqrt(N), .groups = "keep") %>%
  mutate(Treatment = paste0(Diet, "_", Infection, "_", Clone))
m_longdata_resp_sum$Week <- as.factor(m_longdata_resp_sum$Week)
m_longdata_resp_sum$Diet <- factor(m_longdata_resp_sum$Diet, levels = c("S", "SM", "M", "M+"))


# Figure S5C: Plot of metabolism over time by week of the Fungus experiment
plot_byWeek_resp_metsch <- ggplot(data = m_longdata_resp_sum, aes(x = Week, y = resp.mean, group = Treatment, shape = Diet)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_errorbar(aes(x=Week, ymin=resp.mean-resp.se, ymax=resp.mean+resp.se,color=Diet), width=0.2, show.legend = F) + 
  geom_line(aes(color=Diet, group = Treatment), linewidth=1, show.legend = F) +
  geom_point(aes(color=Diet), size=2.5, alpha = 0.8, show.legend = F) +
  facet_grid(Clone~Infection) +
  ylim(-0.002, 0.01) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  annotate("rect", xmin = 0.9, xmax = 1.1, ymin = -0.001, ymax = 0.01,
           alpha = .2,fill = "gray") +
  #ggtitle("Metsch Experiment Metabolic Rate") +
  labs(x = "Week of Experiment", y = "Metabolic Rate (J/hr)") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
plot_byWeek_resp_metsch
ggsave(here("figures/MetschMetabolism_DietxCloneXInfection_byWeek.tiff"), plot = plot_byWeek_resp_metsch, dpi = 300, width = 7, height = 6, units = "in", compression="lzw")




### Pasteuria metabolic rate analyses 

p_data_resp <- p_data %>% dplyr::filter(!is.na(metabolic.rate.Week1))
hist(as.numeric(p_data_resp$metabolic.rate.Week1)) # pretty close to normal, a tad right skewed
range(p_data_resp$metabolic.rate.Week1)

# Metabolic rate during week 1 (exposure period) only
prespmod <- lmer(metabolic.rate.Week1 ~ Diet * Infection * Clone +(1|Block), data = p_data_resp)
summary(prespmod)

# Appendix Table S2
Anova(prespmod)   # use this for the manuscript
AIC(prespmod)
plot(prespmod)
qqnorm(resid(prespmod))  # looks good
qqline(resid(prespmod))
shapiro.test(resid(prespmod)) # not sig


prespmod_contrasts <- emmeans(prespmod, specs = pairwise ~ Diet | Infection + Clone, type = "response")
prespmod_contrasts

# Appendix Table S4
prespmod_contrasts2 <- emmeans(prespmod, specs = pairwise ~ Infection | Clone, type = "response")
prespmod_contrasts2

prespmod_contrasts3 <- emmeans(prespmod, specs = pairwise ~ Diet, type = "response")
prespmod_contrasts3





# Figure S5B: Metabolic rate in week 1 for Bacterium experiment
p_resp_fig <- ggplot(p_data_resp, aes(x = Infection, y = metabolic.rate.Week1)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width =0.9)) +
  geom_point(aes(color=Diet), size=1.5, position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.2), alpha = 0.4) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninfected", "Exposed", "Infected"), labels = c("Uninf", "Exp", "Inf")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  ggtitle("Bacterium Experiment") +
  labs(x = "Infection Status", y = bquote("Metabolic Rate (J/hr)")) +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
p_resp_fig
ggsave("figures/PastMetabolicRate_DietxClonexInfection_Week1.tiff", plot = p_resp_fig, dpi = 600, width = 5.25, height = 4, units = "in", compression="lzw")


# Metabolic rate over whole experiment (time series)
# Figure only, model did not fit well

# Use long data set for time series data of metabolic rate
p_longdata_resp <- p_longdata %>% dplyr::filter(!is.na(metabolic.rate))
hist(as.numeric(p_longdata_resp$metabolic.rate)) # a bit right skewed
range(p_longdata_resp$metabolic.rate)


p_longdata_resp_sum <- p_longdata_resp %>%
  group_by(Diet, Infection, Clone, Week) %>%
  dplyr::summarize(N  = length(metabolic.rate),
                   resp.mean = mean(metabolic.rate, na.rm = T),
                   resp.sd   = sd(metabolic.rate, na.rm = T),
                   resp.se   = resp.sd / sqrt(N)) %>%
  mutate(Treatment = paste0(Diet, "_", Infection, "_", Clone))
p_longdata_resp_sum$Week <- as.factor(p_longdata_resp_sum$Week)
p_longdata_resp_sum$Diet <- factor(p_longdata_resp_sum$Diet, levels = c("S", "SM", "M", "M+"))


# Figure S5D: Plot of metabolism over time by week of the Bacterium experiment
plot_byWeek_resp_past <- ggplot(data = p_longdata_resp_sum, aes(x = Week, y = resp.mean, group = Treatment, shape = Diet)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_errorbar(aes(x=Week, ymin=resp.mean-resp.se, ymax=resp.mean+resp.se,color=Diet), width=0.2) + 
  geom_line(aes(color=Diet, group = Treatment), linewidth=1) +
  geom_point(aes(color=Diet), size=2.5, alpha = 0.8) +
  facet_grid(Clone~Infection) +
  ylim(-0.002, 0.01) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  annotate("rect", xmin = 0.9, xmax = 1.1, ymin = -0.001, ymax = 0.01,
           alpha = .2,fill = "gray") +
  #ggtitle("Past Experiment Metabolic Rate") +
  labs(x = "Week of Experiment", y = "Metabolic Rate (J/hr)") +
  theme_classic() + 
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
plot_byWeek_resp_past
ggsave(here("figures/PastMetabolism_DietxCloneXInfection_byWeek.tiff"), plot = plot_byWeek_resp_past, dpi = 300, width = 7, height = 6, units = "in", compression="lzw")









FigureS5 <- ggarrange(m_resp_fig, p_resp_fig, plot_byWeek_resp_metsch, plot_byWeek_resp_past, labels = "auto",
                      ncol = 2, nrow = 2, widths = c(3.3,4))
FigureS5
ggsave("figures/manuscript/FigS5_MetabolicRate.tiff", plot = FigureS5, dpi = 600, width = 8.5, height = 7, units = "in", compression="lzw")







# Analysis of microcystin-LR  concentrations among treatments and between experiments --------------------

# Appendix Section S1.1

#### Microcystin data (both experiments) - read in data and update order of factors
microcystin <- read.csv("data/Microcystin/Microcystin_ELISA_Samples.csv", stringsAsFactors = F)
head(microcystin)
microcystin$Diet <- as.factor(microcystin$Diet)
microcystin$Diet <- factor(microcystin$Diet, levels = c("S", "SM", "M", "M+"))
microcystin$Parasites <- as.factor(dplyr::recode(microcystin$Parasite, Uninf = "Uninfected"))
microcystin$Infection <- as.factor(dplyr::recode(microcystin$Infection.Status, Uninf = "Uninfected"))
microcystin$Clone <- as.factor(dplyr::recode(microcystin$Clone, MID = "Mid37", STD = "Standard"))
microcystin$Rep <- as.factor(microcystin$Rep)






# look at microcystin concentration across all treatments during parasite exposure (day 8) only
microcystin_exposure <- microcystin %>%
  filter(Day == 8 & !is.na(Microcystin.conc))



hist(microcystin_exposure$Microcystin.conc)
# adjust zero values to be very low but non-zero to be able to run model
microcystin_exposure$Microcystin.conc[ microcystin_exposure$Microcystin.conc == 0] <- 0.01

m_conc_mod2 <- lm(log(Microcystin.conc+1) ~ Experiment * Diet,data = microcystin_exposure)  # use this model
summary(m_conc_mod2)
Anova(m_conc_mod2)   # reported in Appendix Section S1.1
plot(m_conc_mod2)
qqnorm(resid(m_conc_mod2))   # falls off the bottom a bit
qqline(resid(m_conc_mod2))
shapiro.test(resid(m_conc_mod2))


treatment_difs2 <- emmeans(m_conc_mod2, pairwise ~ Diet | Experiment, type = "response")
treatment_difs2  
# no significant differences between S, SM, and M diets within either Past or Metsch experiments
# M+ diet has sig higher Microcystin concentrations compared to all other diets in both experiments

experiment_difs <- emmeans(m_conc_mod2, pairwise ~ Experiment | Diet, type = "response")
experiment_difs  
# no significant differences between concentrations across the experiment at each diet treatment




# fig with boxplots and raw data points for microcystin concentations at day 8 (Not shown in manuscript or appendix)
microcystin_plot <- ggplot(microcystin_exposure, aes(x = Parasite, y = Microcystin.conc+0.001, fill = Diet)) +
  geom_boxplot(position=position_dodge()) +
  geom_jitter(aes(shape = Diet), size = 2, alpha = 0.5, position=position_jitterdodge()) +
  facet_grid(Experiment~Clone, labeller = label_both) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_shape_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  ggtitle("Microcystin concentrations during exposure, Day 8") +
  labs(x = "Parasite Treatments", y = "log microcystin concentation (ug/L)") +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1.0, 10.0), labels = c(0.001, 0.01, 0.1, 1.0, 10.0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
microcystin_plot
#ggsave(here("figures/Microcystin_exposure_at_Day8.tiff"), plot = microcystin_plot, dpi = 300, width = 8, height = 6, units = "in", compression="lzw")


# calculate average microcystin concentration per diet,clone, parasite treatment in each experiment
microcystin_exposure.sum <- microcystin_exposure %>%
  group_by(Experiment, Diet, Clone, Parasites) %>%
  dplyr::summarize(N  = length(Microcystin.conc), 
                   microcystin.avg = mean(Microcystin.conc, na.rm = T),
                   microcystin.sd   = sd(Microcystin.conc, na.rm = T),
                   microcystin.se   = microcystin.sd / sqrt(N))
microcystin_exposure.sum$Exposure <- as.factor(dplyr::recode(microcystin_exposure.sum$Parasites, Metsch = "Parasite", Pasteuria = "Parasite"))
microcystin_exposure.sum$Experiment <- as.factor(dplyr::recode(microcystin_exposure.sum$Experiment, Metsch = "Fungus", Pasteuria = "Bacterium"))
microcystin_exposure.sum
write.csv(microcystin_exposure.sum, "tables/MicrocystinConc_DuringExposure.csv", quote = F, row.names=FALSE)
microcystin_exposure.sum[ 22, "microcystin.avg"] <- 0.01 # alter single value where average is zero (to be able to log transform the y axis in the figure below)

# Figure S1: Avg Microcystin during parasite exposure (day 8)
microcystin_plot2 <- ggplot(microcystin_exposure.sum, aes(x = Exposure, y = microcystin.avg)) +
  geom_point(aes(color = Diet), position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin=microcystin.avg-microcystin.se, ymax=microcystin.avg+microcystin.se,color=Diet), width=0.3, position=position_dodge(width = 0.5)) + 
  #geom_jitter(aes(shape = Diet), size = 2, alpha = 0.5, position=position_jitterdodge()) +
  facet_grid(Experiment~Clone, labeller = label_both) +
  scale_color_discrete(limits = c("S", "SM", "M", "M+")) +
  #scale_shape_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  #scale_shape_manual(values=c(15,16,17,18)) +
  #ggtitle("Microcystin concentrations during exposure, Day 8") +
  labs(x = "Parasite Treatments", y = bquote("log average microcystin concentation (" ~mu *"g/L)")) +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1.0, 10.0), labels = c(0.001, 0.01, 0.1, 1.0, 10.0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
microcystin_plot2
ggsave(here("figures/manuscript/FigS1_Microcystin_exposure_average_at_Day8.tiff"), plot = microcystin_plot2, dpi = 300, width = 4, height = 5, units = "in", compression="lzw")





# Analysis of the effects of the super toxin spike at day 5 of experiment
# in comparison to during exposure (day 8)
super_spike <- microcystin %>%
  filter(Day == 5 | Day == 8) %>%
  filter(Diet == "M" | Diet == "M+") %>%
  filter(Experiment == "Pasteuria") %>%
  filter(!is.na(Microcystin.conc))
# diet color palette
diet_colors_micro <- c("#006837", "#1D91C0")
super_spike$Day <- factor(super_spike$Day)

# Fig S2: comparison of the super spike concentrations to those during exposure in M and M+ treatments
spike_plot <- ggplot(super_spike, aes(x = Diet, y = Microcystin.conc, fill = Diet)) +
  geom_boxplot() +
  geom_jitter(size = 2, alpha = 0.5, width = 0.2) +
  facet_wrap(~Day, labeller = label_both) +
  scale_fill_manual(values = diet_colors_micro) +
  #ggtitle("Microcystin spiked concentrations, Day 5 and 8") +
  labs(x = "Diet", y = bquote("log average microcystin concentation (" ~mu *"g/L)")) +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1.0, 10.0, 100.0), labels = c(0.001, 0.01, 0.1, 1.0, 10.0, 100.0))
spike_plot
ggsave(here("figures/manuscript/FigS2_Microcystin_super_spike_at_Day5.tiff"), plot = spike_plot, dpi = 300, width = 5, height = 4, units = "in", compression="lzw")

# remove outliers for analysis of differences during toxin spike (improve normality of residuals)
super_spike2 <- filter(super_spike, Microcystin.conc < 290)
super_spike2 <- filter(super_spike2, Microcystin.conc < 3 | Microcystin.conc > 5 )

m_conc_mod <- lm(log(Microcystin.conc) ~ Day*Diet, data = super_spike2)
summary(m_conc_mod)
Anova(m_conc_mod)  # reported in Appendix Section S1.1
plot(m_conc_mod)
qqnorm(resid(m_conc_mod))   # falls off the bottom a bit
qqline(resid(m_conc_mod))
shapiro.test(resid(m_conc_mod))  # not significant
treatment_difs <- emmeans(m_conc_mod, pairwise ~ Diet * Day, type = "response")
treatment_difs   # reported in Appendix Section S1.1
# no sig diff between M treatments on days 5 and 8.
# sig differences between M and M+ on day 5, M and M+ on day 8, and M+ on days 5 and 8

