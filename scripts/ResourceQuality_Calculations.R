
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
library(ggplot2)


########################################
### Load Pasteuria data:
########################################

#### lifespan, offspring, infection prevalence
offspring_past <- read.csv("ResourceQuality_Pasteuria_Survival_Offspring.csv", stringsAsFactors = F)
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
bodysize_past <- read.csv("ResourceQuality_Pasteuria_Bodysize.csv", stringsAsFactors = F)
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
spores_past <- read.csv("ResourceQuality_Pasteuria_SporeCounts.csv", stringsAsFactors = F)
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

