# Resource Quality Calculations and combining data sets

# Code associated with Fearon et al. "Resource Quality Differentially Impacts Daphnia 
# Interactions with Two Parasites" submitted to Ecological Monographs.

# This code combines all sources of experimental data into single data sets for the Resource Quality with Pasteuria and Metschnikowia experiments
# Both short and long (by experimental day) versions of the data are created for use in different analyses.

# Code written by: Michelle Fearon
# Last updated: July 7, 2024

# Terminology Used
    # Bacterium = Pasteuria ramosa
    # Fungus = Metschnikowia bicuspidata (Metsch)



# load libraries ------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(lubridate)
library(survival) #library to do survival analyses
library(rms) #library to make survival curves
library(survminer)  #library to do survival analyses
library(emmeans)
library(ggplot2)
library(car)
library(olsrr)
library(DHARMa)
library(performance)
library(glmmTMB)
library(here)


# set the path to the script relative to the project root directory
here::i_am("scripts/ResourceQuality_Calculations.R")




# Pasteuria (Bacterium) data --------------------------------------------------------


## lifespan, offspring, infection prevalence ---------------------------------------
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
  dplyr::select(Unique.code:Block,lifespan.days, Infection, InfectionStatus, Babies_7:Babies_38) %>%
  pivot_longer(Babies_7:Babies_38, names_to = "Day", names_prefix = "Babies_", values_to = "Offspring")
 

# add week to dataset based on day of experiment
offspring_past_long$ID <- 1:length(offspring_past_long$Block)
offspring_past_long$Week <- NA

for(i in offspring_past_long$ID){
  if(offspring_past_long$Day[i] == 7 | offspring_past_long$Day[i] == 8){
    offspring_past_long$Week[i] <- 1
  }else if(offspring_past_long$Day[i] > 8 | offspring_past_long$Day[i] <= 15){
    offspring_past_long$Week[i] <- 2
  }else if(offspring_past_long$Day[i] > 15 & offspring_past_long$Day[i] <= 22){
    offspring_past_long$Week[i] <- 3
  }else if(offspring_past_long$Day[i] > 22 & offspring_past_long$Day[i] <= 29){
    offspring_past_long$Week[i] <- 4
  }else if(offspring_past_long$Day[i] > 29 & offspring_past_long$Day[i] <= 36){
    offspring_past_long$Week[i] <- 5
  }else if(offspring_past_long$Day[i] > 36){
    offspring_past_long$Week[i] <- 5.5
  }else{
    offspring_past_long$Week[i] <- NA
  }
}



# remove day and offspring by day from this data set
offspring_past_long2 <- dplyr::select(offspring_past_long, Unique.code:InfectionStatus, Week)

# Calculate the number of offspring produced each week and cumulatively by week
offspring_past_longWeek <- offspring_past_long %>%
  group_by(Unique.code, Diet, Parasites, Clone, Rep, Block, Week) %>%
  dplyr::summarize(WeeklyOffspring = sum(Offspring, na.rm = T)) %>%
  mutate(CumulativeWeeklyOffspring = cumsum(WeeklyOffspring))

# join data with other lifespan data
offspring_past_longWeek <- left_join(offspring_past_longWeek, offspring_past_long2)

# re-order columns and remove duplicates
offspring_past_longWeek <- offspring_past_longWeek %>%
  dplyr::select(Unique.code:Week, Infection, InfectionStatus, lifespan.days, WeeklyOffspring, CumulativeWeeklyOffspring) %>%
  distinct()
dim(offspring_past_longWeek)
offspring_past_longWeek$Week <- as.factor(offspring_past_longWeek$Week)
#View(offspring_past_longWeek)


# create short version of the data set with one row per replicate
offspring_past_short <- offspring_past %>%
  dplyr::select(Unique.code:Treatment, Total.Babies, lifespan.days, SurvObj)



## bodysize data ------------------------------------------------------------------
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
  

# add week to dataset based on day of experiment
bodysize_past_long$ID <- 1:length(bodysize_past_long$Block)
bodysize_past_long$Week <- NA

for(i in bodysize_past_long$ID){
  if(bodysize_past_long$Day[i] == 7 | bodysize_past_long$Day[i] == 8){
    bodysize_past_long$Week[i] <- 1
  }else if(bodysize_past_long$Day[i] == 14 | bodysize_past_long$Day[i] == 15){
    bodysize_past_long$Week[i] <- 2
  }else if(bodysize_past_long$Day[i] == 21 | bodysize_past_long$Day[i] == 22){
    bodysize_past_long$Week[i] <- 3
  }else if(bodysize_past_long$Day[i] == 28 | bodysize_past_long$Day[i] == 29){
    bodysize_past_long$Week[i] <- 4
  }else if(bodysize_past_long$Day[i] == 35 | bodysize_past_long$Day[i] == 36){
    bodysize_past_long$Week[i] <- 5
  }else if(bodysize_past_long$Day[i] == 38){
    bodysize_past_long$Week[i] <- 5.5
    }else{
    bodysize_past_long$Week[i] <- NA
  }
}


bodysize_past_longWeek <- bodysize_past_long %>%
  group_by(Unique.code, Diet, Parasites, Clone, Rep, Block, Week) %>%
  dplyr::summarize(WeeklyLength = mean(Length, na.rm = T)) %>%
  mutate(WeeklyLength = ifelse(is.nan(WeeklyLength), NA, WeeklyLength))
bodysize_past_longWeek$Week <- as.factor(bodysize_past_longWeek$Week)



# create short version of the data set with one row per replicate
        # calculate growth rate between the first and second week of experiment ( 7 to 14 and 8 to 15, respectively)
bodysize_past_short <- bodysize_past %>%
  mutate(GrowthRate1 = (length_14 - length_7)/7,  GrowthRate2 = (length_15 - length_8)/7, GrowthRate = if_else(is.na(GrowthRate1), GrowthRate2, GrowthRate1),
         Length_Week1 = if_else(is.na(length_7), length_8, length_7)) %>%
  dplyr::select(Unique.code:Block, Length_Week1, GrowthRate)





## spore yield data ----------------------------------------------------------------
spores_past <- read.csv("data/ResourceQuality_Pasteuria_SporeCounts.csv", stringsAsFactors = F)
head(spores_past)
spores_past$Diet <- as.factor(spores_past$Diet)
spores_past$Parasites <- as.factor(spores_past$Parasites)
spores_past$Clone <- as.factor(spores_past$Clone)
spores_past$Rep <- as.factor(spores_past$Rep)
spores_past$Block <- as.factor(spores_past$Block)


spores_past_short <- spores_past %>% 
  # create total spores column for each quadrant counted
  mutate(total_1 = rowSums(dplyr::select(., contains("_1"))), total_2 = rowSums(dplyr::select(., contains("_2"))), total_3 = rowSums(dplyr::select(., contains("_3"))), total_4 = rowSums(dplyr::select(., contains("_4")))) %>%
  # average spore counts across four quadrants for each developmental stage, then multiply by 10,000 to get spores/mL, then multiply by 0.1 mL (10,000 x 0.1 = 1000) = # spores per animal
  mutate(mature_spores = (rowSums(dplyr::select(., starts_with("mature")))/4) *1000,
         immature_spores = (rowSums(dplyr::select(., contains("immature")))/4) *1000, 
         caul_spores = (rowSums(dplyr::select(., contains("caul")))/4) *1000, 
         total_spores = (rowSums(dplyr::select(., contains("total")))/4) *1000) %>%
  dplyr::select(Unique.code:Block, mature_spores:total_spores)
#View(spores_past_short)


## Spore size ---------------------------------------------------------------------
sporesize_past <- read.csv("data/ResourceQuality_Pasteuria_SporeSize.csv", stringsAsFactors = F)
head(sporesize_past)
sporesize_past$Diet <- as.factor(sporesize_past$Diet)
sporesize_past$Diet <- factor(sporesize_past$Diet, levels = c("S", "SM", "M", "M+"))
sporesize_past$Parasites <- as.factor(sporesize_past$Parasites)
sporesize_past$Clone <- as.factor(sporesize_past$Clone)
sporesize_past$Rep <- as.factor(sporesize_past$Rep)



sporesize_past_sum <- sporesize_past %>%
  mutate(SporeArea = pi  * (HorizontalDiameter/2) * (VerticalDiameter/2)) %>%
  filter(!is.na(SporeArea)) %>%
  dplyr::group_by(Diet, Parasites, Clone, Rep) %>%
  dplyr::summarize(AvgSporeSize = mean(SporeArea))
sporesize_past_sum$Diet <- factor(sporesize_past_sum$Diet, levels = c("S", "SM", "M", "M+"))
hist(sporesize_past_sum$AvgSporeSize)




## feeding rate data -------------------------------------------------------------
feed_data <- read.csv("data/ResourceQuality_FeedingRateCalc.csv", stringsAsFactors = F)  
head(feed_data)
feed_data$Diet <- as.factor(feed_data$Diet)
feed_data$Parasites <- as.factor(dplyr::recode(feed_data$Parasites, Uninf = "Uninfected"))
feed_data$Clone <- as.factor(dplyr::recode(feed_data$Clone, MID = "Mid37", STD = "Standard"))
feed_data$Rep <- as.factor(feed_data$Rep)
feed_data$Block <- as.factor(feed_data$Block)
feed_data$Week <- as.factor(feed_data$Week)

# filter by experiment, and remove dead animals and blanks, generate short and long versions of data
feed_past_data_long <- feed_data %>% 
  filter(Experiment == "Past") %>%
  dplyr::select(Diet, Parasites, Clone, Rep, Block, Week, Clearance, Clearance_rel)

feed_past_data_short <- feed_data %>% 
  filter(Experiment == "Past", Week == 1) %>%
  dplyr::select(Diet, Parasites, Clone, Rep, Block, Clearance, Clearance_rel) %>%
  dplyr::rename(Clearance.Week1 = Clearance, Clearance_rel.Week1 = Clearance_rel)




## respiration data ----------------------------------------------------------------
resp_data <- read.csv("data/ResourceQuality_RespirationRateCalc.csv", stringsAsFactors = F)  
head(resp_data)
resp_data$Diet <- as.factor(resp_data$Diet)
resp_data$Parasites <- as.factor(dplyr::recode(resp_data$Parasites, Uninf = "Uninfected"))
resp_data$Clone <- as.factor(dplyr::recode(resp_data$Clone, MID = "Mid37", STD = "Standard"))
resp_data$Rep <- as.factor(resp_data$Rep)
resp_data$Block <- as.factor(resp_data$Block)
resp_data$Week <- as.factor(resp_data$Week)

# filter by experiment, and remove dead animals, mislabeled samples, and 1 bad read, and remove blanks, generate short and long versions of data
resp_past_data_long <- resp_data %>% 
  filter(Experiment == "Past", Died.After.Resp == "N", Clone != "Blank") %>%
  dplyr::select(Diet, Parasites, Clone, Rep, Block, Week, metabolic.rate)

resp_past_data_short <- resp_data %>% 
  filter(Experiment == "Past", Died.After.Resp == "N", Clone != "Blank", Week == 1) %>%
  dplyr::select(Diet, Parasites, Clone, Rep, Block, metabolic.rate) %>%
  dplyr::rename(metabolic.rate.Week1 = metabolic.rate)




## Microcystin data (both experiments) -------------------------------------------
microcystin <- read.csv("data/Microcystin/Microcystin_ELISA_Samples.csv", stringsAsFactors = F)
head(microcystin)
microcystin$Diet <- as.factor(microcystin$Diet)
microcystin$Diet <- factor(microcystin$Diet, levels = c("S", "SM", "M", "M+"))
microcystin$Parasites <- as.factor(dplyr::recode(microcystin$Parasite, Uninf = "Uninfected"))
microcystin$Infection <- as.factor(dplyr::recode(microcystin$Infection.Status, Uninf = "Uninfected"))
microcystin$Clone <- as.factor(dplyr::recode(microcystin$Clone, MID = "Mid37", STD = "Standard"))
microcystin$Rep <- as.factor(microcystin$Rep)



# look at microcystin concentration across all treatments during parasite exposure (day 8)
microcystin_exposure <- microcystin %>%
  filter(Day == 8 & !is.na(Microcystin.conc))


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



## table of average Microcystin concentrations per Diet treatment in each experiment
microcystin_exposure$Diet2 <- as.factor(dplyr::recode(microcystin_exposure$Diet, S = "non-toxin", SM = "non-toxin", M = "non-toxin"))
microcystin_exposure.sum2 <- microcystin_exposure %>%
  group_by(Experiment, Diet2) %>%
  dplyr::summarize(N  = length(Microcystin.conc), 
                   microcystin.avg = mean(Microcystin.conc, na.rm = T),
                   microcystin.sd   = sd(Microcystin.conc, na.rm = T),
                   microcystin.se   = microcystin.sd / sqrt(N))
View(microcystin_exposure.sum2)
write.csv(microcystin_exposure.sum2, "tables/AverageMicrocystinConc_AllDietsVsToxinDiet" , quote = F, row.names=FALSE)

# calculate average microcystin concentration per diet, clone, parasite infection status treatment to join with larger data set
microcystin_exposure.sum3 <- microcystin_exposure %>%
  group_by(Experiment, Diet, Clone, Parasites, Infection) %>%
  dplyr::summarize(N  = length(Microcystin.conc), 
                   microcystin.avg = mean(Microcystin.conc, na.rm = T),
                   microcystin.sd   = sd(Microcystin.conc, na.rm = T),
                   microcystin.se   = microcystin.sd / sqrt(N))
View(microcystin_exposure.sum3)

past_microcystin <- microcystin_exposure.sum3 %>% 
  filter(Experiment == "Pasteuria") %>%
  ungroup() %>%
  dplyr::select(Diet:Infection, microcystin.avg)






## Join bacterium short data sets together ------------------------------------
past_data_short <- full_join(offspring_past_short, bodysize_past_short)
past_data_short <- full_join(past_data_short, spores_past_short)
past_data_short <- full_join(past_data_short, sporesize_past_sum)
past_data_short <- full_join(past_data_short, feed_past_data_short)
past_data_short <- full_join(past_data_short, resp_past_data_short)
past_data_short <- full_join(past_data_short, past_microcystin)
View(past_data_short)
dim(past_data_short)  # should have 160 rows, 1 for each animal in the study
write.csv(past_data_short, "data/ResourceQuality_Pasteuria_Full.csv", quote = F, row.names=FALSE)



## Join bacterium long data sets together -------------------------------------
past_data_longWeek <- full_join(offspring_past_longWeek, bodysize_past_longWeek)
past_data_longWeek <- full_join(past_data_longWeek, feed_past_data_long)
past_data_longWeek <- full_join(past_data_longWeek, resp_past_data_long)
head(past_data_longWeek)
dim(past_data_longWeek)
View(past_data_longWeek)

write.csv(past_data_longWeek, "data/ResourceQuality_Pasteuria_Full_ByExptWeek.csv", quote = F, row.names=FALSE)





  
########################################
# Metsch (Fungus) data --------------------------------------------------------------



#### lifespan, offspring, infection prevalence ---------------------------------------
offspring_metsch <- read.csv("data/ResourceQuality_Metsch_Survival_Offspring.csv", stringsAsFactors = F)
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
  dplyr::select(Unique.code:Block,lifespan.days, Infection, InfectionStatus, Babies_7:Babies_24) %>%
  pivot_longer(Babies_7:Babies_24, names_to = "Day", names_prefix = "Babies_", values_to = "Offspring")
View(offspring_metsch_long)  



# add week to dataset based on day of experiment
offspring_metsch_long$ID <- 1:length(offspring_metsch_long$Block)
offspring_metsch_long$Week <- NA

for(i in offspring_metsch_long$ID){
  if(offspring_metsch_long$Day[i] == 7 | offspring_metsch_long$Day[i] == 8){
    offspring_metsch_long$Week[i] <- 1
  }else if(offspring_metsch_long$Day[i] > 8 | offspring_metsch_long$Day[i] <= 15){
    offspring_metsch_long$Week[i] <- 2
  }else if(offspring_metsch_long$Day[i] > 15 & offspring_metsch_long$Day[i] <= 22){
    offspring_metsch_long$Week[i] <- 3
  }else if(offspring_metsch_long$Day[i] > 22){
    offspring_metsch_long$Week[i] <- 3.5
  }else{
    offspring_metsch_long$Week[i] <- NA
  }
}



# remove day and offspring by day from this data set
offspring_metsch_long2 <- dplyr::select(offspring_metsch_long, Unique.code:InfectionStatus, Week)

# Calculate the number of offspring produced each week and cumulatively by week
offspring_metsch_longWeek <- offspring_metsch_long %>%
  group_by(Unique.code, Diet, Parasites, Clone, Rep, Block, Week) %>%
  dplyr::summarize(WeeklyOffspring = sum(Offspring, na.rm = T)) %>%
  mutate(CumulativeWeeklyOffspring = cumsum(WeeklyOffspring))

# join data with other lifespan data
offspring_metsch_longWeek <- left_join(offspring_metsch_longWeek, offspring_metsch_long2)

# re-order columns and remove duplicates
offspring_metsch_longWeek <- offspring_metsch_longWeek %>%
  dplyr::select(Unique.code:Week, Infection, InfectionStatus, lifespan.days, WeeklyOffspring, CumulativeWeeklyOffspring) %>%
  distinct()
dim(offspring_metsch_longWeek)
offspring_metsch_longWeek$Week <- as.factor(offspring_metsch_longWeek$Week)
#View(offspring_metsch_longWeek)



# create short version of the data set with one row per replicate
offspring_metsch_short <- offspring_metsch %>%
  dplyr::select(Unique.code:Treatment, Total.Babies, lifespan.days, SurvObj)






## bodysize data ----------------------------------------------------------------
bodysize_metsch <- read.csv("data/ResourceQuality_Metsch_Bodysize.csv", stringsAsFactors = F)
head(bodysize_metsch)
bodysize_metsch$Diet <- as.factor(bodysize_metsch$Diet)
bodysize_metsch$Parasites <- as.factor(bodysize_metsch$Parasites)
bodysize_metsch$Clone <- as.factor(bodysize_metsch$Clone)
bodysize_metsch$Rep <- as.factor(bodysize_metsch$Rep)
bodysize_metsch$Block <- as.factor(bodysize_metsch$Block)


# create a long version of the data set based on experimental day
bodysize_metsch_long <- bodysize_metsch %>%
  pivot_longer(length_7:length_24, names_to = "Day", names_prefix = "length_", values_to = "Length")
  


# add week to dataset based on day of experiment
bodysize_metsch_long$ID <- 1:length(bodysize_metsch_long$Block)
bodysize_metsch_long$Week <- NA

for(i in bodysize_metsch_long$ID){
  if(bodysize_metsch_long$Day[i] == 7 | bodysize_metsch_long$Day[i] == 8){
    bodysize_metsch_long$Week[i] <- 1
  }else if(bodysize_metsch_long$Day[i] == 14 | bodysize_metsch_long$Day[i] == 15){
    bodysize_metsch_long$Week[i] <- 2
  }else if(bodysize_metsch_long$Day[i] == 21 | bodysize_metsch_long$Day[i] == 22){
    bodysize_metsch_long$Week[i] <- 3
  }else if(bodysize_metsch_long$Day[i] == 24){
    bodysize_metsch_long$Week[i] <- 3.5
  }else{
    bodysize_metsch_long$Week[i] <- NA
  }
}


bodysize_metsch_longWeek <- bodysize_metsch_long %>%
  group_by(Unique.code, Diet, Parasites, Clone, Rep, Block, Week) %>%
  dplyr::summarize(WeeklyLength = mean(Length, na.rm = T)) %>%
  mutate(WeeklyLength = ifelse(is.nan(WeeklyLength), NA, WeeklyLength))

bodysize_metsch_longWeek$Week <- as.factor(bodysize_metsch_longWeek$Week)


# create short version of the data set with one row per replicate
# calculate growth rate between the first and second week of experiment (7 to 14 and 8 to 15, respectively)
bodysize_metsch_short <- bodysize_metsch %>%
  mutate(GrowthRate1 = (length_14 - length_7)/7,  GrowthRate2 = (length_15 - length_8)/7, GrowthRate = if_else(is.na(GrowthRate1), GrowthRate2, GrowthRate1),
         Length_Week1 = if_else(is.na(length_7), length_8, length_7)) %>%
  dplyr::select(Unique.code:Block, Length_Week1, GrowthRate)




## Metsch gut spore data -----------------------------------------------------------
gutspore_metsch <- read.csv("data/ResourceQuality_Metsch_GutSpores.csv", stringsAsFactors = F)
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
  dplyr::select(!Notes)





## spore yield data -----------------------------------------------------------------
spores_metsch <- read.csv("data/ResourceQuality_Metsch_SporeCounts.csv", stringsAsFactors = F)
head(spores_metsch)
spores_metsch$Diet <- as.factor(spores_metsch$Diet)
spores_metsch$Parasites <- as.factor(spores_metsch$Parasites)
spores_metsch$Clone <- as.factor(spores_metsch$Clone)
spores_metsch$Rep <- as.factor(spores_metsch$Rep)
spores_metsch$Block <- as.factor(spores_metsch$Block)


spores_metsch <- spores_metsch %>% 
  # create total spores column for each quadrant counted
  mutate(total_1 = rowSums(dplyr::select(., contains("_1"))), total_2 = rowSums(dplyr::select(., contains("_2"))), total_3 = rowSums(dplyr::select(., contains("_3"))), total_4 = rowSums(dplyr::select(., contains("_4")))) %>%
  # average spore counts across four quadrants for each developmental stage, then multiply by 10,000 to get spores/mL, then multiply by 0.1 mL (10,000 x 0.1 = 1000) = # spores per animal
  mutate(mature_spores = (rowSums(dplyr::select(., starts_with("mature")))/4) * (10000 * Volume),
         immature_spores = (rowSums(dplyr::select(., contains("immature")))/4) * (10000 * Volume), 
         bud_spores = (rowSums(dplyr::select(., contains("bud")))/4) * (10000 * Volume), 
         total_spores = (rowSums(dplyr::select(., contains("total")))/4) * (10000 * Volume))


spores_metsch_short <- spores_metsch %>%
  dplyr::select(Unique.code:Block, mature_spores:total_spores)
View(spores_metsch_short)


## feeding rate data -------------------------------------------------------------
        # data was loaded in above, just need to filter by Metsch experiment 

# filter by experiment, and remove dead animals and blanks, generate short and long versions of data
feed_metsch_data_long <- feed_data %>% 
  filter(Experiment == "Metsch") %>%
  dplyr::select(Diet, Parasites, Clone, Rep, Block, Week, Clearance, Clearance_rel)

feed_metsch_data_short <- feed_data %>% 
  filter(Experiment == "Metsch", Week == 1) %>%
  dplyr::select(Diet, Parasites, Clone, Rep, Block, Clearance, Clearance_rel) %>%
  dplyr::rename(Clearance.Week1 = Clearance, Clearance_rel.Week1 = Clearance_rel)


## respiration data ----------------------------------------------------------------
      # data was loaded in above, just need to filter by Metsch experiment 
unique(resp_data$Experiment)
# filter by experiment, and remove dead animals and blanks, generate short and long versions of data
resp_metsch_data_long <- resp_data %>% 
  filter(Experiment == "Metsch", Died.After.Resp != "Y",Clone != "Blank") %>%
  dplyr::select(Diet, Parasites, Clone, Rep, Block, Week, metabolic.rate)

resp_metsch_data_short <- resp_data %>% 
  filter(Experiment == "Metsch", Died.After.Resp == "N", Clone != "Blank", Week == 1) %>%
  dplyr::select(Diet, Parasites, Clone, Rep, Block, metabolic.rate) %>%
  dplyr::rename(metabolic.rate.Week1 = metabolic.rate)




## Microcystin data -------------------------------------------------------------
metsch_microcystin <- microcystin_exposure.sum3 %>% 
  filter(Experiment == "Metsch") %>%
  ungroup() %>%
  dplyr::select(Diet:Clone, Infection, microcystin.avg)
unique(metsch_microcystin$Clone)
View(metsch_microcystin)





## Join fungus short data sets together --------------------------------------------------
metsch_data_short <- full_join(offspring_metsch_short, bodysize_metsch_short)
metsch_data_short <- full_join(metsch_data_short, spores_metsch_short)
metsch_data_short <- full_join(metsch_data_short, gutspore_metsch)
metsch_data_short <- full_join(metsch_data_short, feed_metsch_data_short)
metsch_data_short <- full_join(metsch_data_short, resp_metsch_data_short)
metsch_data_short <- left_join(metsch_data_short, metsch_microcystin, by = c(c("Diet", "Clone", "Infection")))
unique(metsch_data_short$Clone)
dim(metsch_data_short)   # should be 240 rows, one for each animal in the experiment
write.csv(metsch_data_short, "data/ResourceQuality_Metsch_Full.csv", quote = F, row.names=FALSE)

View(metsch_data_short)


## Join fungus long data sets together -----------------------------------------------------
metsch_data_longWeek <- full_join(offspring_metsch_longWeek, bodysize_metsch_longWeek)
metsch_data_longWeek <- full_join(metsch_data_longWeek, feed_metsch_data_long)
metsch_data_longWeek <- full_join(metsch_data_longWeek, resp_metsch_data_long)
head(metsch_data_longWeek)
View(metsch_data_longWeek)
write.csv(metsch_data_longWeek, "data/ResourceQuality_Metsch_Full_ByExptWeek.csv", quote = F, row.names=FALSE)

