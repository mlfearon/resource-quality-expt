
setwd("C:/Users/mlfea/OneDrive/Documents/PROJECTS/MHMP Daphnia Duffy/Resource Quality/Data and Code")


# load libraries
library(dplyr)
library(tidyr)
library(lubridate)
#library to do survival analyses
library(survival)
#library to make survival curves
library(rms)
library(survminer)
library(emmeans)
library(ggplot2)
library(here)


# set the path to the script relative to the project root directory
here::i_am("scripts/ResourceQuality_Calculations.R")

########################################
### Load Pasteuria data:
########################################

#### lifespan, offspring, infection prevalence
offspring_past <- read.csv("data/ResourceQuality_Pasteuria_Survival_Offspring.csv", stringsAsFactors = F)
head(offspring_past)

# update data classes to factor class
offspring_past$Diet <- as.factor(offspring_past$Diet)
offspring_past$Parasites <- as.factor(offspring_past$Parasites)
offspring_past$Clone <- as.factor(offspring_past$Clone)
offspring_past$Rep <- as.factor(offspring_past$Rep)
offspring_past$Block <- as.factor(offspring_past$Block)
offspring_past$Include_Lifespan <- as.factor(offspring_past$Include_Lifespan)
offspring_past$Infection <- as.factor(offspring_past$Infection)
offspring_past$Treatment <- as.factor(offspring_past$Treatment)
summary(offspring_past)
unique(offspring_past$Treatment)

# set experiment beginning and death date format
emerge <- mdy(offspring_past$DateEmerged)
death <- mdy(offspring_past$DateDied)

time.interval <- emerge %--% death
time.duration <- as.duration(time.interval)
offspring_past$lifespan.days <- as.numeric(time.duration, "days")

offspring_past$SurvObj <- with(offspring_past, Surv(lifespan.days, Status))


# create a long version of the data set based on experimental day
offspring_past_long <- offspring_past %>%
  select(Unique.code:Block,lifespan.days, Infection, InfectionStatus, Babies_7:Babies_38) %>%
  pivot_longer(Babies_7:Babies_38, names_to = "Day", names_prefix = "Babies_", values_to = "Offspring")
View(offspring_past_long)  


# create short version of the data set with one row per replicate
offspring_past_short <- offspring_past %>%
  select(Unique.code:Treatment, Total.Babies, lifespan.days, SurvObj)




### bodysize data
bodysize_past <- read.csv("data/ResourceQuality_Pasteuria_Bodysize.csv", stringsAsFactors = F)
head(bodysize_past)
bodysize_past$Diet <- as.factor(bodysize_past$Diet)
bodysize_past$Parasites <- as.factor(bodysize_past$Parasites)
bodysize_past$Clone <- as.factor(bodysize_past$Clone)
bodysize_past$Rep <- as.factor(bodysize_past$Rep)
bodysize_past$Block <- as.factor(bodysize_past$Block)


# create a long version of the data set based on experimental day
bodysize_past_long <- bodysize_past %>%
  pivot_longer(length_7:length_38, names_to = "Day", names_prefix = "length_", values_to = "Length")
View(bodysize_past_long)  


# create short version of the data set with one row per replicate
        # calculate growth rate between the first and second week of experiment ( 7 to 14 and 8 to 15, respectively)
bodysize_past_short <- bodysize_past %>%
  mutate(GrowthRate1 = (length_14 - length_7)/7,  GrowthRate2 = (length_15 - length_8)/7, GrowthRate = if_else(is.na(GrowthRate1), GrowthRate2, GrowthRate1)) %>%
  select(Unique.code:Block, GrowthRate)





### spore yield data
spores_past <- read.csv("data/ResourceQuality_Pasteuria_SporeCounts.csv", stringsAsFactors = F)
head(spores_past)
spores_past$Diet <- as.factor(spores_past$Diet)
spores_past$Parasites <- as.factor(spores_past$Parasites)
spores_past$Clone <- as.factor(spores_past$Clone)
spores_past$Rep <- as.factor(spores_past$Rep)
spores_past$Block <- as.factor(spores_past$Block)


spores_past_short <- spores_past %>% 
  # create total spores column for each quadrant counted
  mutate(total_1 = rowSums(select(., contains("_1"))), total_2 = rowSums(select(., contains("_2"))), total_3 = rowSums(select(., contains("_3"))), total_4 = rowSums(select(., contains("_4")))) %>%
  # average spore counts across four quadrants for each developmental stage, then multiply by 10,000 to get spores/mL, then multiply by 0.1 mL (10,000 x 0.1 = 1000) = # spores per animal
  mutate(mature_spores = (rowSums(select(., starts_with("mature")))/4) *1000,
         immature_spores = (rowSums(select(., contains("immature")))/4) *1000, 
         caul_spores = (rowSums(select(., contains("caul")))/4) *1000, 
         total_spores = (rowSums(select(., contains("total")))/4) *1000) %>%
  select(Unique.code:Block, mature_spores:total_spores)
View(spores_past_short)




### feeding rate data




#### respiration data 
  
  



#### Microcystin data (both experiments)
microcystin <- read.csv("data/Microcystin/Microcystin_ELISA_Samples.csv", stringsAsFactors = F)
head(microcystin)
microcystin$Diet <- factor(microcystin$Diet, levels = c("S", "SM", "M", "M+"))

# check super toxin spike at day 5 of experiment
super_spike <- microcystin %>%
  filter(Day == 5 | Day == 8) %>%
  filter(Diet == "M" | Diet == "M+") %>%
  filter(Experiment == "Pasteuria") %>%
  filter(!is.na(Microcystin.conc))
# diet color palette
diet_colors_micro <- c("#006837", "#1D91C0")
super_spike$Day <- factor(super_spike$Day)

Day_labels <- c("Day 5, pre-exposure", "Day 8, parasite exposure")
names(Day_labels) <- c(5, 8)

# compare the super spike concentrations to those during exposure 
spike_plot <- ggplot(super_spike, aes(x = Diet, y = Microcystin.conc, fill = Diet)) +
  geom_boxplot() +
  geom_jitter(size = 2, alpha = 0.5, width = 0.2) +
  facet_wrap(~Day, labeller = label_both) +
  scale_fill_manual(values = diet_colors_micro) +
  ggtitle("Microcystin spiked concentrations, Day 5 and 8") +
  labs(x = "Diet", y = "log Microcystin concentation (ug/L)") +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1.0, 10.0, 100.0), labels = c(0.001, 0.01, 0.1, 1.0, 10.0, 100.0))
spike_plot
ggsave(here("figures/Microcystin_super_spike_at_Day5.tiff"), plot = spike_plot, dpi = 300, width = 5, height = 4, units = "in", compression="lzw")

m_conc_mod <- lm(log(Microcystin.conc) ~ Day*Diet, data = super_spike)
summary(m_conc_mod)
plot(m_conc_mod)
treatment_difs <- emmeans(m_conc_mod, pairwise ~ Diet * Day, type = "response")
treatment_difs  
        # no sig diff between M treatments on days 5 and 8.
        # sig differences between M and M+ on day 5, M and M+ on day 8, and M+ on days 5 and 8



# look at microcystin concentration across all treatments during parasite exposure (day 8)
microcystin_exposure <- microcystin %>%
  filter(Day == 8 & !is.na(Microcystin.conc))

# diet color palette
diet_colors <- c("#ADDD8E", "#41AB5D", "#006837", "#1D91C0")

microcystin_plot <- ggplot(microcystin_exposure, aes(x = Parasite, y = Microcystin.conc+0.001, fill = Diet)) +
  geom_boxplot(position=position_dodge()) +
  geom_jitter(aes(shape = Diet), size = 2, alpha = 0.5, position=position_jitterdodge()) +
  facet_grid(Experiment~Clone, labeller = label_both) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_shape_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  ggtitle("Microcystin concentrations during exposure, Day 8") +
  labs(x = "Parasite Treatments", y = "log Microcystin concentation (ug/L)") +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1.0, 10.0), labels = c(0.001, 0.01, 0.1, 1.0, 10.0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
microcystin_plot
ggsave(here("figures/Microcystin_exposure_at_Day8.tiff"), plot = microcystin_plot, dpi = 300, width = 8, height = 6, units = "in", compression="lzw")


# calculate average microcystin concentration per diet,clone, parasite treatment in each experiment
microcystin_exposure.sum <- microcystin_exposure %>%
  group_by(Experiment, Diet, Clone, Parasite) %>%
  dplyr::summarize(N  = length(Microcystin.conc), 
                   microcystin.avg = mean(Microcystin.conc, na.rm = T),
                   microcystin.sd   = sd(Microcystin.conc, na.rm = T),
                   microcystin.se   = microcystin.sd / sqrt(N))
microcystin_exposure.sum
microcystin_exposure.sum[ 22, "microcystin.avg"] <- 0.01 # alter single value where average is zero


microcystin_plot2 <- ggplot(microcystin_exposure.sum, aes(x = Parasite, y = microcystin.avg)) +
  geom_point(aes(color = Diet), position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin=microcystin.avg-microcystin.se, ymax=microcystin.avg+microcystin.se,color=Diet), width=0.3, position=position_dodge(width = 0.5)) + 
  #geom_jitter(aes(shape = Diet), size = 2, alpha = 0.5, position=position_jitterdodge()) +
  facet_grid(Experiment~Clone, labeller = label_both) +
  scale_color_discrete(limits = c("S", "SM", "M", "M+")) +
  #scale_shape_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  #scale_shape_manual(values=c(15,16,17,18)) +
  ggtitle("Microcystin concentrations during exposure, Day 8") +
  labs(x = "Parasite Treatments", y = "log Average Microcystin concentation (ug/L)") +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1.0, 10.0), labels = c(0.001, 0.01, 0.1, 1.0, 10.0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
microcystin_plot2
ggsave(here("figures/Microcystin_exposure_average_at_Day8.tiff"), plot = microcystin_plot2, dpi = 300, width = 8, height = 6, units = "in", compression="lzw")

## table of average Microcystin concentrations per Diet treatment in each experiment
microcystin_exposure.sum2 <- microcystin_exposure %>%
  group_by(Experiment, Diet) %>%
  dplyr::summarize(N  = length(Microcystin.conc), 
                   microcystin.avg = mean(Microcystin.conc, na.rm = T),
                   microcystin.sd   = sd(Microcystin.conc, na.rm = T),
                   microcystin.se   = microcystin.sd / sqrt(N))
microcystin_exposure.sum2


# calculate average microcystin concentration per diet, clone, parasite infection status treatment to join with larger data set
microcystin_exposure.sum3 <- microcystin_exposure %>%
  group_by(Experiment, Diet, Clone, Parasite, Infection.Status) %>%
  dplyr::summarize(N  = length(Microcystin.conc), 
                   microcystin.avg = mean(Microcystin.conc, na.rm = T),
                   microcystin.sd   = sd(Microcystin.conc, na.rm = T),
                   microcystin.se   = microcystin.sd / sqrt(N))
View(microcystin_exposure.sum3)

past_microcystin <- filter(microcystin_exposure.sum3, Experiment == "Pasteuria")
metsch_microcystin <- filter(microcystin_exposure.sum3, Experiment == "Metsch")




##### Join short data sets together
past_data_short <- full_join(offspring_past_short, bodysize_past_short)
past_data_short <- full_join(past_data_short, spores_past_short)
View(past_data_short)
write.csv(past_data_short, "ResourceQuality_Pasteuria_Full.csv", quote = F, row.names=FALSE)


# calculate sample sizes for water samples at each date
past_summary <- past_data_short %>%
  filter(Lifespan >= 17) %>%
  group_by(Diet, Clone, Infection) %>%
  summarise(sample.size.total = n())

#View(past_summary)
write.csv(past_summary, "ResourceQuality_Pasteuria_Samplesize_Day17.csv", quote = F, row.names=FALSE)

past_summary


#### Join long data sets together
past_data_long <- full_join(offspring_past_long, bodysize_past_long)
head(past_data_long)

past_data_long <- past_data_long %>% 
  mutate( Week = case_when(Day %% 7 == 0 ~ "fizz buzz" 1, 
                           Day %% 8 == 1,
                           Day %% 10 == 1,
                           Day %% 12 == 1,
                           Day %% 14 == 2,
                           Day %% 15 == 2,
                           Day %% 16 == 2,
                           Day %% 19 == 2,
                           Day %% 21 == 3,
                           Day %% 22 == 3,
                           Day %% 23 == 3,
                           Day %% 26 == 3,
                           Day %% 28 == 4,
                           Day %% 29 == 4,
                           Day %% 30 == 4,
                           Day %% 33 == 4,
                           Day %% 35 == 5,
                           Day %% 36 == 5,
                           Day %% 37 == 5,
                           Day %% 38 == 5))

write.csv(past_data_long, "ResourceQuality_Pasteuria_Full_ByExptDay.csv", quote = F, row.names=FALSE)





  
########################################
### Load Metsch data:
########################################


#### lifespan, offspring, infection prevalence
offspring_metsch <- read.csv("ResourceQuality_Metsch_Survival_Offspring.csv", stringsAsFactors = F)
head(offspring_metsch)

# update data classes to factor class
offspring_metsch$Diet <- as.factor(offspring_metsch$Diet)
offspring_metsch$Parasites <- as.factor(offspring_metsch$Parasites)
offspring_metsch$Clone <- as.factor(offspring_metsch$Clone)
offspring_metsch$Rep <- as.factor(offspring_metsch$Rep)
offspring_metsch$Block <- as.factor(offspring_metsch$Block)
offspring_metsch$Include_Lifespan <- as.factor(offspring_metsch$Include_Lifespan)
offspring_metsch$Infection <- as.factor(offspring_metsch$Infection)
offspring_metsch$Treatment <- as.factor(offspring_metsch$Treatment)
summary(offspring_metsch)
unique(offspring_metsch$Treatment)

# set experiment beginning and death date format
emerge2 <- mdy(offspring_metsch$DateEmerged)
death2 <- mdy(offspring_metsch$DateDied)

time.interval2 <- emerge2 %--% death2
time.duration2 <- as.duration(time.interval2)
offspring_metsch$lifespan.days <- as.numeric(time.duration2, "days")

offspring_metsch$SurvObj <- with(offspring_metsch, Surv(lifespan.days, Status))


# create a long version of the data set based on experimental day
offspring_metsch_long <- offspring_metsch %>%
  select(Unique.code:Block,lifespan.days, Infection, InfectionStatus, Babies_7:Babies_24) %>%
  pivot_longer(Babies_7:Babies_24, names_to = "Day", names_prefix = "Babies_", values_to = "Offspring")
View(offspring_metsch_long)  


# create short version of the data set with one row per replicate
offspring_metsch_short <- offspring_metsch %>%
  select(Unique.code:Treatment, Total.Babies, lifespan.days, SurvObj)






### bodysize data
bodysize_metsch <- read.csv("ResourceQuality_Metsch_Bodysize.csv", stringsAsFactors = F)
head(bodysize_metsch)
bodysize_metsch$Diet <- as.factor(bodysize_metsch$Diet)
bodysize_metsch$Parasites <- as.factor(bodysize_metsch$Parasites)
bodysize_metsch$Clone <- as.factor(bodysize_metsch$Clone)
bodysize_metsch$Rep <- as.factor(bodysize_metsch$Rep)
bodysize_metsch$Block <- as.factor(bodysize_metsch$Block)


# create a long version of the data set based on experimental day
bodysize_metsch_long <- bodysize_metsch %>%
  pivot_longer(length_7:length_24, names_to = "Day", names_prefix = "length_", values_to = "Length")
View(bodysize_metsch_long)  


# create short version of the data set with one row per replicate
# calculate growth rate between the first and second week of experiment (7 to 14 and 8 to 15, respectively)
bodysize_metsch_short <- bodysize_metsch %>%
  mutate(GrowthRate1 = (length_14 - length_7)/7,  GrowthRate2 = (length_15 - length_8)/7, GrowthRate = if_else(is.na(GrowthRate1), GrowthRate2, GrowthRate1)) %>%
  select(Unique.code:Block, GrowthRate)




# Metsch gut spore data
gutspore_metsch <- read.csv("ResourceQuality_Metsch_GutSpores.csv", stringsAsFactors = F)
head(gutspore_metsch)
gutspore_metsch$Diet <- as.factor(gutspore_metsch$Diet)
gutspore_metsch$Parasites <- as.factor(gutspore_metsch$Parasites)
gutspore_metsch$Clone <- as.factor(gutspore_metsch$Clone)
gutspore_metsch$Rep <- as.factor(gutspore_metsch$Rep)



gutspore_metsch <- gutspore_metsch %>%
  mutate(TotalGutSpores = Embedded.spores + Punctured.spores + Haemoceol.spores,
         TotalPunctSpores = Punctured.spores + Haemoceol.spores,
         HemocytesPerSpore = if_else(TotalPunctSpores > 0, Hemocytes/TotalPunctSpores, 0)) %>%
  filter(Notes != "not exposed") %>%
  select(!Notes)


View(gutspore_metsch)



### spore yield data
spores_metsch <- read.csv("ResourceQuality_Metsch_SporeCounts.csv", stringsAsFactors = F)
head(spores_metsch)
spores_metsch$Diet <- as.factor(spores_metsch$Diet)
spores_metsch$Parasites <- as.factor(spores_metsch$Parasites)
spores_metsch$Clone <- as.factor(spores_metsch$Clone)
spores_metsch$Rep <- as.factor(spores_metsch$Rep)
spores_metsch$Block <- as.factor(spores_metsch$Block)


spores_metsch <- spores_metsch %>% 
  # create total spores column for each quadrant counted
  mutate(total_1 = rowSums(select(., contains("_1"))), total_2 = rowSums(select(., contains("_2"))), total_3 = rowSums(select(., contains("_3"))), total_4 = rowSums(select(., contains("_4")))) %>%
  # average spore counts across four quadrants for each developmental stage, then multiply by 10,000 to get spores/mL, then multiply by 0.1 mL (10,000 x 0.1 = 1000) = # spores per animal
  mutate(mature_spores = (rowSums(select(., starts_with("mature")))/4) * (10000 * Volume),
         immature_spores = (rowSums(select(., contains("immature")))/4) * (10000 * Volume), 
         bud_spores = (rowSums(select(., contains("bud")))/4) * (10000 * Volume), 
         total_spores = (rowSums(select(., contains("total")))/4) * (10000 * Volume))


spores_metsch_short <- spores_metsch %>%
  select(Unique.code:Block, mature_spores:total_spores)
View(spores_metsch_short)


### feeding rate data




#### respiration data 



#### Microcystin data

      # see above in the Pasteuria section for figures and calculations of averages per treatment



##### Join short data sets together
metsch_data_short <- full_join(offspring_metsch_short, bodysize_metsch_short)
metsch_data_short <- full_join(metsch_data_short, spores_metsch_short)
metsch_data_short <- full_join(metsch_data_short, gutspore_metsch)


write.csv(metsch_data_short, "ResourceQuality_Metsch_Full.csv", quote = F, row.names=FALSE)

View(metsch_data_short)


# calculate sample sizes for water samples at each date
metsch_data_short$Rep2 <- as.integer(metsch_data_short$Rep)

metsch_summary <- metsch_data_short %>%
  filter(Rep2 < 11, Lifespan > 6) %>%
  group_by(Diet, Clone, Infection) %>%
  summarise(sample.size.total = n())

#View(metsch_summary)
#write.csv(metsch_summary, "ResourceQuality_Metsch_Samplesize_Day6.csv", quote = F, row.names=FALSE)


#### Join long data sets together
metsch_data_long <- full_join(offspring_metsch_long, bodysize_metsch_long)


write.csv(metsch_data_long, "ResourceQuality_Metsch_Full_ByExptDay.csv", quote = F, row.names=FALSE)

