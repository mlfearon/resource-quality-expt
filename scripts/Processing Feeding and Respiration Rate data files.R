# Processing Feeding rate and Respiration Rate data files

#Michelle's working directory


library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyverse)
library(lubridate)
library(EnvStats)
library(emmeans)
library(fs)
#install_github('colin-olito/LoLinR')
# https://colin-olito.github.io/LoLinR/vignettes/LoLinR.html  # how to use this package link
library(LoLinR)
library(here)



# set the path to the script relative to the project root directory
here::i_am("scripts/Processing Feeding and Respiration Rate data files.R")


################################################
# Feeding Rate
################################################


# returns the file names as a character vector
file_paths <- fs::dir_ls(path = "./data/FeedingRate/")


# read in all files in the directory, remove first row, add block, day and plate metadata for each dataset, combine all rows
FeedingRate <- file_paths %>%
  map( function(path) {
    read_csv(path, skip=1)} %>% # read in all files and remove the first row
      mutate(Name = str_remove_all(path, "./data/FeedingRate/"))) %>%   # add file name to each row
  bind_rows %>% # bind all rows together
  separate(Name, into=c("Block", "Day", "Plate"), sep="_") %>%  # separate the file name into block, day and plate numbers
  separate(Plate, into=c("Plate", "csv"), sep="[.]") %>%   # separate the plate number and csv file extension
  dplyr::rename(Date.Time = "Reading Date/Time", chlorophyll = "chlorphyll:435,665") %>%
  extract(Well, into=c("Row", "Column"), regex = "(^[[:alpha:]])([[:digit:]]+)", convert = TRUE, remove = FALSE) %>% # separate well row and column info
  select(Well:Plate)

#View(FeedingRate)

# recode block, day and plate to numbers
FeedingRate$Block <- dplyr::recode(FeedingRate$Block, Block1 = 1, Block2 = 2, Block3 = 3, Block4 = 4)
FeedingRate$Day <- dplyr::recode(FeedingRate$Day, Day7 = 7, Day8 = 8, Day14 = 14, Day15 = 15, Day21 = 21, Day22 = 22, Day28 = 28, Day29 = 29, Day35 = 35, Day36 = 36)
FeedingRate$Plate <- dplyr::recode(FeedingRate$Plate, Plate1 = 1, Plate2 = 2, Plate3 = 3, Plate4 = 4, Plate5 = 5, Plate6 = 6, Plate7 = 7, Plate8 = 8)


# need to join meta data to the feeding rate data
past_feed_meta <- read.csv(here("data/ResourceQuality_Pasteuria_FeedingRateMetadata.csv"), header = T, stringsAsFactors = F)
metsch_feed_meta <- read.csv(here("data/ResourceQuality_Metsch_FeedingRateMetadata.csv"), header = T, stringsAsFactors = F)

feeding_meta <- bind_rows(past_feed_meta, metsch_feed_meta)

FeedingRate <- full_join(FeedingRate, feeding_meta, by = c("Block", "Day", "Plate", "Column"))



# Calculate end time for Pasteuria experiment (subtract 20 minutes from the plate read time stamp to estimate end time of feeding trial)
plateread <- mdy_hm(FeedingRate$Date.Time)
FeedingRate$newplateread <- format(as.POSIXct(plateread - minutes(20)), format = "%H:%M")
FeedingRate$EndTime <- format(FeedingRate$EndTime, format = "%H:%M")
FeedingRate$ID <- 1:length(FeedingRate$Well.x)

for(i in FeedingRate$ID){
  if(FeedingRate$EndTime[i] == "     "){
    FeedingRate$EndTime[i] <- format(FeedingRate$newplateread[i], format = "%H:%M")
  }
}


# calculate the duration of the feeding trial from start and end times (in hours)
FeedingRate$StartTime <- as.POSIXct(FeedingRate$StartTime, tz = "EST", format = "%H:%M") # adds in today's date, but only need the start/end times relative to each other
FeedingRate$EndTime <- as.POSIXct(FeedingRate$EndTime, tz = "EST", format = "%H:%M")     # adds in today's date, but only need the start/end times relative to each other

time.interval <- interval(start = FeedingRate$StartTime, end = FeedingRate$EndTime)
time.duration <- as.duration(time.interval)
FeedingRate$Time <- as.numeric(time.duration, "hours")

# for instances where the calculated time is less than 6 hours, adjust the total time to 6 hours
# because none of the plates were stopped before 6 hours
for(j in FeedingRate$ID){
  if(FeedingRate$Time[j] < 6){
    FeedingRate$Time[j] <- 6.0
  }
}

sum(FeedingRate$Time < 6)

# filter out fluorescence readings that don't have a sample associated with it
FeedingRate <- filter(FeedingRate, !is.na(Clone))




## Checking for outliers
hist(FeedingRate$chlorophyll[FeedingRate$Diet == "S"])
hist(FeedingRate$chlorophyll[FeedingRate$Diet == "SM"])
hist(FeedingRate$chlorophyll[FeedingRate$Diet == "M"])
hist(FeedingRate$chlorophyll[FeedingRate$Diet == "M+"])
hist(FeedingRate$chlorophyll[FeedingRate$Clone == "NoDaphnia"])
hist(FeedingRate$chlorophyll[FeedingRate$Clone == "NoAlgae"])


# get dataset with only treatments
FeedingRate_treatments <- filter(FeedingRate, Clone != "NoAlgae")
FeedingRate_treatments <- filter(FeedingRate_treatments, Clone != "NoDaphnia")

# look for outliers in data
boxplot(FeedingRate_treatments$chlorophyll ~ FeedingRate_treatments$Diet)
          # looks like there are at least a couple of outliers in some of the diet treatments

# get data set for each diet to check for outliers within each diet treatment
FeedingRate_S <- filter(FeedingRate_treatments, Diet == "S")
FeedingRate_SM <- filter(FeedingRate_treatments, Diet == "SM")
FeedingRate_M <- filter(FeedingRate_treatments, Diet == "M")
FeedingRate_M_toxin <- filter(FeedingRate_treatments, Diet == "M+")

# test for outlier in S diet
outlier_test_S <- rosnerTest(FeedingRate_S$chlorophyll, k = 10)
outlier_test_S$all.stats
      # six outliers in treatment data, with values between 73 and 58

# test for outlier in SM diet
outlier_test_SM <- rosnerTest(FeedingRate_SM$chlorophyll, k = 5)
outlier_test_SM$all.stats
    # one outlier in treatment data, value at 40

# test for outlier in M diet
outlier_test_M <- rosnerTest(FeedingRate_M$chlorophyll, k = 5)
outlier_test_M$all.stats
    # one outlier in treatment data, value at 21

# test for outlier in M+ diet
outlier_test_Mtoxin <- rosnerTest(FeedingRate_M_toxin$chlorophyll, k = 5)
outlier_test_Mtoxin$all.stats
        # test indicates that there are two outliers, values of 45 and 20.

# Look at updated boxplot by diet treatment
boxplot(FeedingRate$chlorophyll ~ FeedingRate$Diet)
      

###################### REMOVED ONLY THE MOST EXTRME OUTLIERS -- NEED TO DECIDE ON HOW MANY OTHERS TO REMOVE FROM THE TESTS ABOVE!

# remove  10 outliers from treatment data
#FeedingRate <- as.data.frame(FeedingRate)
dim(FeedingRate)
FeedingRate[ FeedingRate$Diet == "S" & FeedingRate$chlorophyll == 73, "chlorophyll"] <- NA
FeedingRate[ FeedingRate$Diet == "S" & FeedingRate$chlorophyll == 61 & !is.na(FeedingRate$chlorophyll), "chlorophyll"] <- NA
FeedingRate[ FeedingRate$Diet == "S" & FeedingRate$chlorophyll == 59 & !is.na(FeedingRate$chlorophyll), "chlorophyll"] <- NA
FeedingRate[ FeedingRate$Diet == "S" & FeedingRate$chlorophyll == 58 & !is.na(FeedingRate$chlorophyll), "chlorophyll"] <- c(NA, NA, NA)
FeedingRate[ FeedingRate$Diet == "SM" & FeedingRate$chlorophyll == 40, "chlorophyll"] <- NA
FeedingRate[ FeedingRate$Diet == "M" & FeedingRate$chlorophyll == 21, "chlorophyll"] <- NA
FeedingRate[ FeedingRate$Diet == "M+" & FeedingRate$chlorophyll == 45, "chlorophyll"] <- NA
FeedingRate[ FeedingRate$Diet == "M+" & FeedingRate$chlorophyll == 20 & !is.na(FeedingRate$chlorophyll), "chlorophyll"] <- NA
FeedingRate <- filter(FeedingRate, !is.na(chlorophyll))
dim(FeedingRate) # checked that it worked, should have 9214 rows


# get data sets with NoDaphnia and No algae controls
FeedingRate_NoAlgae <- filter(FeedingRate, Clone == "NoAlgae")
FeedingRate_NoDaphnia <- filter(FeedingRate, Clone == "NoDaphnia")

# look for outliers
boxplot(FeedingRate_NoDaphnia$chlorophyll ~ FeedingRate_NoDaphnia$Diet)  # The no daphnia measurements vary based on diet treatment, so check for outliers within each diet type
boxplot(FeedingRate_NoAlgae$chlorophyll ~ FeedingRate_NoAlgae$Diet)

# test for outlier in NoAlgae data
outlier_test_noalgae <- rosnerTest(FeedingRate_NoAlgae$chlorophyll, k = 10)
outlier_test_noalgae$all.stats
# several outliers in the No Algae controls - we don't use this for any calculations so it is already left out


# filter No Daphnia controls by diet treatment to test for outliers within each diet type b/c 
FeedingRate_NoDaphnia_S <- filter(FeedingRate_NoDaphnia, Diet == "S")
FeedingRate_NoDaphnia_SM <- filter(FeedingRate_NoDaphnia, Diet == "SM")
FeedingRate_NoDaphnia_M <- filter(FeedingRate_NoDaphnia, Diet == "M")
FeedingRate_NoDaphnia_M_toxin <- filter(FeedingRate_NoDaphnia, Diet == "M+")

# test for outlier in NoDaphnia data (ALL CONTROLS)
outlier_test_nodaph <- rosnerTest(FeedingRate_NoDaphnia$chlorophyll, k = 5)
outlier_test_nodaph$all.stats
      # no outliers in the No Daphnia controls

# test for outliers in S diet NoDaphnia data
outlier_test_nodaph_S <- rosnerTest(FeedingRate_NoDaphnia_S$chlorophyll, k = 5)
outlier_test_nodaph_S$all.stats
      # no outliers in the S diet No Daphnia controls

# test for outlier in SM diet NoDaphnia data
outlier_test_nodaph_SM <- rosnerTest(FeedingRate_NoDaphnia_SM$chlorophyll, k = 5)
outlier_test_nodaph_SM$all.stats
      # no outliers in the SM diet No Daphnia controls

# test for outlier in M diet NoDaphnia data
outlier_test_nodaph_M <- rosnerTest(FeedingRate_NoDaphnia_M$chlorophyll, k = 5)
outlier_test_nodaph_M$all.stats
      # no outliers in the M diet No Daphnia controls

# test for outlier in M+ diet NoDaphnia data
outlier_test_nodaph_M_toxin <- rosnerTest(FeedingRate_NoDaphnia_M_toxin$chlorophyll, k = 5)
outlier_test_nodaph_M_toxin$all.stats
      # no outliers in the M+ diet No Daphnia controls




# remove daphnia that died during feeding rate trial
FeedingRate_died <- filter(FeedingRate, Died.After.Feed == "Y")
#View(FeedingRate_died)
FeedingRate <- filter(FeedingRate, Died.After.Feed != "Y")
dim(FeedingRate)

# quick figures to check the distribution of raw data points for each treatment across weeks
FeedingRate_diets <- filter(FeedingRate, Clone != "NoDaphnia", Clone != "NoAlgae")
check_data <- ggplot(data = FeedingRate_diets, aes(x = Day, y = chlorophyll, color = Diet, shape = Diet)) +
  geom_jitter(alpha = 0.8, position = position_dodge(width = 1)) +
  facet_grid(Experiment ~ Parasites)
check_data
ggsave(here("figures/FeedingRateRAW_DietxParasitexExperiment.tiff"), plot = check_data, dpi = 300, width = 8, height = 6, units = "in", compression="lzw")


# quick figures to check the distribution of raw data points for controls for each treatment across weeks
FeedingRate_controls <- filter(FeedingRate, Clone == "NoDaphnia")
check_controls <- ggplot(data = FeedingRate_controls, aes(x = Day, y = chlorophyll, color = Diet, shape = Diet)) +
  geom_jitter(alpha = 0.8, position = position_dodge(width = 1)) +
  facet_grid(Experiment ~ Parasites)
check_controls
ggsave(here("figures/FeedingRateCONTROLS_DietxParasitexExperiment.tiff"), plot = check_controls, dpi = 300, width = 8, height = 6, units = "in", compression="lzw")




# Summarize data by block, day, plate, column and date/time
    # calculate a mean for the chlorophyll reading based on the groups so that 
    # all 8 fluorescence reading replicates are averaged together. Also calculate the median and St. Dev.

FeedingRateCalc <- FeedingRate %>%
  group_by(Experiment, Treatment, Block, Day, Plate, Column, Time, Diet, Parasites, Clone, Rep) %>%
  dplyr::summarize(chlorophyll.avg = mean(chlorophyll),
            chlorophyll.med = median(chlorophyll),
            chlorophyll.sd = sd(chlorophyll))

FeedingRateCalc_NoDaphnia <- filter(FeedingRateCalc, Clone == "NoDaphnia")
FeedingRateCalc_NoAlgae <- filter(FeedingRateCalc, Clone == "NoAlgae")
FeedingRateCalc_Treatment <- filter(FeedingRateCalc, Clone != "NoDaphnia", Clone != "NoAlgae")


#View(FeedingRateCalc_Treatment)
class(FeedingRateCalc_Treatment$Day)
FeedingRateCalc_Treatment$Rep <- as.integer(FeedingRateCalc_Treatment$Rep)


# add in bodysize for each daphnia
### bodysize data
bodysize_past <- read.csv(here("data/ResourceQuality_Pasteuria_Bodysize.csv"), stringsAsFactors = F)
head(bodysize_past)
bodysize_metsch <- read.csv(here("data/ResourceQuality_Metsch_Bodysize.csv"), stringsAsFactors = F)
head(bodysize_metsch)



# create a long version of the data set based on experimental day
bodysize_past_long <- bodysize_past %>%
  pivot_longer(length_7:length_38, names_to = "Day", names_prefix = "length_", values_to = "Length") %>%
  filter(!is.na(Length), Day != 38) %>%  # remove days when there are no bodysize measurements, remove last date where there are no corresponding feeding rate
  select(!Unique.code)
#View(bodysize_past_long)


bodysize_metsch_long <- bodysize_metsch %>%
  pivot_longer(length_7:length_24, names_to = "Day", names_prefix = "length_", values_to = "Length") %>%
  filter(!is.na(Length), Day != 24) %>%  # remove days when there are no bodysize measurements, remove last date where there are no corresponding feeding rate
  select(!Unique.code)

bodysize_long <- bind_rows(bodysize_past_long, bodysize_metsch_long)
 

# recode Mid37 and Std labels to match feeding rate
bodysize_long$Clone <- dplyr::recode(bodysize_long$Clone, Mid37 = "MID", Standard = "STD")
bodysize_long$Parasites <- dplyr::recode(bodysize_long$Parasites, Uninfected = "Uninf")
bodysize_long$Day <- as.numeric(bodysize_long$Day)
dim(bodysize_long)

FeedingRateCalc_Treatment <- full_join(FeedingRateCalc_Treatment, bodysize_long)


#### remove daphnia with missing bodysize data
# first check why these samples are missing data
filter(FeedingRateCalc_Treatment, is.na(Length)) # missing animal, error in photo, or dead animal removed from data
filter(FeedingRateCalc_Treatment, is.na(chlorophyll.avg))  # missing feeding data due to dead animals removed or extra bodysize measurements that don't have feeding data associated

# remove samples with missing data
FeedingRateCalc_Treatment <- filter(FeedingRateCalc_Treatment, !is.na(Length), !is.na(chlorophyll.avg))


# add Control NoDaphnia treatments back into FeedingRateCalc dataset as a new column
FeedingRateCalc_NoDaphnia <- FeedingRateCalc_NoDaphnia %>%
  dplyr::rename(chlorophyll.nodaphnia = "chlorophyll.avg", chlorophyll.nodaphnia.med = "chlorophyll.med", chlorophyll.nodaphnia.sd = "chlorophyll.sd") %>%
  ungroup() %>%
  select(-Column, -Clone, -Rep, -Treatment)

FeedingRateCalc_Treatment <- full_join(FeedingRateCalc_Treatment, FeedingRateCalc_NoDaphnia, by = c("Experiment", "Block", "Day", "Plate", "Diet"))


# add week to dataset based on day of experiment
FeedingRateCalc_Treatment$ID <- 1:length(FeedingRateCalc_Treatment$Block)
FeedingRateCalc_Treatment$Week <- NA

for(i in FeedingRateCalc_Treatment$ID){
  if(FeedingRateCalc_Treatment$Day[i] == 7 | FeedingRateCalc_Treatment$Day[i] == 8){
    FeedingRateCalc_Treatment$Week[i] <- 1
  }else if(FeedingRateCalc_Treatment$Day[i] == 14 | FeedingRateCalc_Treatment$Day[i] == 15){
    FeedingRateCalc_Treatment$Week[i] <- 2
  }else if(FeedingRateCalc_Treatment$Day[i] == 21 | FeedingRateCalc_Treatment$Day[i] == 22){
    FeedingRateCalc_Treatment$Week[i] <- 3
  }else if(FeedingRateCalc_Treatment$Day[i] == 28 | FeedingRateCalc_Treatment$Day[i] == 29){
    FeedingRateCalc_Treatment$Week[i] <- 4
  }else if(FeedingRateCalc_Treatment$Day[i] == 35 | FeedingRateCalc_Treatment$Day[i] == 36){
    FeedingRateCalc_Treatment$Week[i] <- 5
  }else{
    FeedingRateCalc_Treatment$Week[i] <- NA
  }
}


# clean up columns
FeedingRateCalc_Treatment <- FeedingRateCalc_Treatment %>%
  select(Experiment:Column, Week, Diet:Rep, Time.x, Length, chlorophyll.avg, chlorophyll.med, chlorophyll.sd, chlorophyll.nodaphnia, chlorophyll.nodaphnia.med, chlorophyll.nodaphnia.sd) %>%
  dplyr::rename(Time = "Time.x", Parasites = "Parasites.x")



##################################
# Calculate feeding rate (per ml, per hr)
    # ln (mean no daphnia control/ mean remaining food in sample) * (Volume 10 mL/ length of time of the assay in hr)


# add in total volume of feeding vials (10 mL)
FeedingRateCalc_Treatment$volume <- 10

# (ln (mean no daphnia control) - ln (mean remaining food in sample)) * (Volume 10 mL/ length of time of the assay in hr)
FeedingRateCalc_Treatment$Clearance <- (log(FeedingRateCalc_Treatment$chlorophyll.nodaphnia) - log(FeedingRateCalc_Treatment$chlorophyll.avg)) * (FeedingRateCalc_Treatment$volume / FeedingRateCalc_Treatment$Time)
FeedingRateCalc_Treatment$Clearance2 <- (log(FeedingRateCalc_Treatment$chlorophyll.nodaphnia.med) - log(FeedingRateCalc_Treatment$chlorophyll.med)) * (FeedingRateCalc_Treatment$volume / FeedingRateCalc_Treatment$Time)


# calculate the feeding rate for relative bodysize (mL -hr -mm)
FeedingRateCalc_Treatment$Clearance_rel <- FeedingRateCalc_Treatment$Clearance / ((FeedingRateCalc_Treatment$Length/1000) ^ 2)
FeedingRateCalc_Treatment$Clearance_rel2 <- FeedingRateCalc_Treatment$Clearance2 / ((FeedingRateCalc_Treatment$Length/1000) ^ 2)

ggplot(data = FeedingRateCalc_Treatment, aes( x= Diet, y = Clearance, color = Parasites)) +
  geom_boxplot() +
  geom_jitter() +
  facet_grid(Clone~Week) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_hline(yintercept = -0.4, color = "red", linetype = "dashed") 
boxplot(FeedingRateCalc_Treatment$Clearance ~ FeedingRateCalc_Treatment$Diet)


write.csv(FeedingRateCalc_Treatment, "data/ResourceQuality_FeedingRateCalc.csv", quote = F, row.names=FALSE)

hist(FeedingRateCalc_Treatment$chlorophyll.sd)
hist(FeedingRateCalc_Treatment$chlorophyll.nodaphnia.sd)

filter(FeedingRateCalc_Treatment, chlorophyll.nodaphnia.sd > 7)


# diet color palette
diet_colors <- c("#ADDD8E", "#41AB5D", "#006837", "#1D91C0")

## Plot of Algae clearance and relative clearnace over time by week
FeedingRateCalc_Sum.clean <- FeedingRateCalc_Treatment %>%
  group_by(Diet, Parasites, Clone, Week, Treatment, Experiment) %>%
  summarize(N  = length(Clearance),
            feed.mean = mean(Clearance, na.rm = T),
            feed.sd   = sd(Clearance, na.rm = T),
            feed.se   = feed.sd / sqrt(N), 
            feed.rel.mean = mean(Clearance_rel, na.rm = T),
            feed.rel.sd   = sd(Clearance_rel, na.rm = T),
            feed.rel.se   = feed.rel.sd / sqrt(N))
FeedingRateCalc_Sum.clean$Week <- as.factor(FeedingRateCalc_Sum.clean$Week)
FeedingRateCalc_Sum.clean$Diet <- factor(FeedingRateCalc_Sum.clean$Diet, levels = c("S", "SM", "M", "M+"))

View(FeedingRateCalc_Sum.clean_past)
str(FeedingRateCalc_Sum.clean_past)
# split figures by past and metsch experiments
FeedingRateCalc_Sum.clean_past <- as.data.frame(filter(FeedingRateCalc_Sum.clean, Experiment == "Past"))
FeedingRateCalc_Sum.clean_metsch <- filter(FeedingRateCalc_Sum.clean, Experiment == "Metsch")

plot_byWeek_feed_past <- ggplot(data = FeedingRateCalc_Sum.clean_past, aes(x = Week, y = feed.mean, group = Treatment, shape = Diet)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_errorbar(aes(x=Week, ymin=feed.mean-feed.se, ymax=feed.mean+feed.se,color=Diet), width=0.2) + 
  geom_line(aes(color=Diet, group = Treatment), linewidth=1.5) +
  geom_point(aes(color=Diet), size=3, alpha = 0.8) +
  facet_grid(Clone~Parasites) +
  #ylim(-0.001, 0.009) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  annotate("rect", xmin = 0.9, xmax = 1.1, ymin = -0.6, ymax = 0.5,
           alpha = .2,fill = "gray") +
  ggtitle("Past Experiment Feeding Rate") +
  labs(x = "Week of Experiment", y = "Average Clearance Rate (mL per hr)") +
  theme_classic()
plot_byWeek_feed_past
ggsave(here("figures/PastFeedingRate_DietxCloneXParasite.tiff"), plot = plot_byWeek_feed_past, dpi = 300, width = 6, height = 6, units = "in", compression="lzw")


plot_byWeek_feed_metsch <- ggplot(data = FeedingRateCalc_Sum.clean_metsch, aes(x = Week, y = feed.mean, group = Treatment, shape = Diet)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_errorbar(aes(x=Week, ymin=feed.mean-feed.se, ymax=feed.mean+feed.se,color=Diet), width=0.2) + 
  geom_line(aes(color=Diet, group = Treatment), linewidth=1.5) +
  geom_point(aes(color=Diet), size=3, alpha = 0.8) +
  facet_grid(Clone~Parasites) +
  #ylim(-0.001, 0.009) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  annotate("rect", xmin = 0.9, xmax = 1.1, ymin = -0.4, ymax = 0.6,
           alpha = .2,fill = "gray") +
  ggtitle("Metsch Experiment Feeding rate") +
  labs(x = "Week of Experiment", y = "Average Clearance Rate (mL per hr)") +
  theme_classic()
plot_byWeek_feed_metsch
ggsave(here("figures/MetschFeedingRate_DietxCloneXParasite.tiff"), plot = plot_byWeek_feed_metsch, dpi = 300, width = 6, height = 6, units = "in", compression="lzw")


Feeding_fig <- grid_arrange_shared_legend(plot_byWeek_feed_past, plot_byWeek_feed_metsch, nrow=1, ncol = 2, position = "right")
ggsave(here("figures/Feeding_byExpt.tiff"), plot = Feeding_fig, dpi = 300, width = 9, height = 6, units = "in", compression="lzw")


### Plots of Relative clearance rate
plot_byWeek_feedrel_past <- ggplot(data = FeedingRateCalc_Sum.clean_past, aes(x = Week, y = feed.rel.mean, group = Treatment, shape = Diet)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_errorbar(aes(x=Week, ymin=feed.rel.mean-feed.rel.se, ymax=feed.rel.mean+feed.rel.se,color=Diet), width=0.2) + 
  geom_line(aes(color=Diet, group = Treatment), linewidth=1.5) +
  geom_point(aes(color=Diet), size=3, alpha = 0.8) +
  facet_grid(Clone~Parasites) +
  #ylim(-0.001, 0.009) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  annotate("rect", xmin = 0.9, xmax = 1.1, ymin = -0.6, ymax = 0.5,
           alpha = .2,fill = "gray") +
  ggtitle("Past Experiment Relative Feeding Rate") +
  labs(x = "Week of Experiment", y = "Average Relative Clearance Rate (mL per hr per mm^2)") +
  theme_classic()
plot_byWeek_feedrel_past
ggsave(here("figures/PastRelativeFeedingRate_DietxCloneXParasite.tiff"), plot = plot_byWeek_feedrel_past, dpi = 300, width = 6, height = 6, units = "in", compression="lzw")


plot_byWeek_feedrel_metsch <- ggplot(data = FeedingRateCalc_Sum.clean_metsch, aes(x = Week, y = feed.rel.mean, group = Treatment, shape = Diet)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_errorbar(aes(x=Week, ymin=feed.rel.mean-feed.rel.se, ymax=feed.rel.mean+feed.rel.se,color=Diet), width=0.2) + 
  geom_line(aes(color=Diet, group = Treatment), linewidth=1.5) +
  geom_point(aes(color=Diet), size=3, alpha = 0.8) +
  facet_grid(Clone~Parasites) +
  #ylim(-0.001, 0.009) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  annotate("rect", xmin = 0.9, xmax = 1.1, ymin = -0.4, ymax = 0.6,
           alpha = .2,fill = "gray") +
  ggtitle("Metsch Experiment Relative Feeding rate") +
  labs(x = "Week of Experiment", y = "Average Relative Clearance Rate (mL per hr per mm^2)") +
  theme_classic()
plot_byWeek_feedrel_metsch
ggsave(here("figures/MetschRelativeFeedingRate_DietxCloneXParasite.tiff"), plot = plot_byWeek_feedrel_metsch, dpi = 300, width = 6, height = 6, units = "in", compression="lzw")


FeedingRelative_fig <- grid_arrange_shared_legend(plot_byWeek_feedrel_past, plot_byWeek_feedrel_metsch, nrow=1, ncol = 2, position = "right")
ggsave(here("figures/RelativeFeeding_byExpt.tiff"), plot = FeedingRelative_fig, dpi = 300, width = 9, height = 6, units = "in", compression="lzw")


################################################
# Respiration Rate
################################################

# read in metadata for respiration rate trials
# need to join meta data to the feeding rate data
past_resp_meta <- read.csv(here("data/ResourceQuality_Pasteuria_RespirationRateMetadata.csv"), header = T, stringsAsFactors = F)
metsch_resp_meta <- read.csv(here("data/ResourceQuality_Metsch_RespirationRateMetadata.csv"), header = T, stringsAsFactors = F)

RespRateCalc <- bind_rows(past_resp_meta, metsch_resp_meta)
unique(RespRateCalc$Parasites)
unique(RespRateCalc$Clone)
unique(RespRateCalc$Diet)
### Need to convert all Oxygen files from excel to csv, remove first 12 rows and last 7 columns (27-32).
## I tried to do this in a fancy way, and I was able to select just the oxygen files, but ran into problems with the 
## excel file type and some of the weird symbols that were in the header of the excel file prevented it from downloading
## Anyways I got frustrated and just did it by hand in the end... try again another day
      # code to select just the Oxygen docs, problem is with reading in the docs as they are due to special symbols in the header

# returns the file names as a character vector
file_paths2 <- fs::dir_ls(path = "./data/RespirationRate/Oxygen_Sat/", regexp = ('(Oxygen)(.csv)$')) # need to select just the oxygen docs
length(file_paths2)


# read in all files in the directory, remove first row, add block, day and plate metadata for each data set, combine all rows
RespRate <- file_paths2 %>%
  map( function(path) {
    read_csv(path, skip=0, col_types = "cncccccccccccccccccccccccc")} %>% # read in all files and remove the first row
      mutate(Name = str_remove_all(path, "./data/RespirationRate/Oxygen_Sat/"))) %>% # add file name to each row
  bind_rows %>% # bind all rows together
  separate(Name, into=c("Block", "Day", "Plate", "results", "Oxygen"), sep="_") %>% # separate the file name into block, day and plate numbers
  dplyr::rename(Date.Time = "Date/Time", Time.Min = "Time/Min.") %>%
  #extract(Well, into=c("Row", "Column"), regex = "(^[[:alpha:]])([[:digit:]]+)", convert = TRUE, remove = FALSE) %>% # separate well row and column info
  select(Block:Plate, Date.Time:D6)

#View(RespRate)

# recode block, day and plate to numbers
RespRate$Block <- dplyr::recode(RespRate$Block, Block1 = 1, Block2 = 2, Block3 = 3, Block4 = 4)
RespRate$Day <- dplyr::recode(RespRate$Day, Day7 = 7, Day8 = 8, Day14 = 14, Day15 = 15, Day21 = 21, Day22 = 22, Day28 = 28, Day29 = 29, Day35 = 35, Day36 = 36)
RespRate$Plate <- dplyr::recode(RespRate$Plate, plate1 = 1, plate2 = 2, plate3 = 3, plate4 = 4)

# convert Date.Time column to a date and time format
RespRate$Date.Time <- as.POSIXct(RespRate$Date.Time, tz = "EST", format = "%d.%m.%Y %H:%M%OS") # adds in today's date, but only need the start/end times relative to each other


# convert all "No Sensor" occurrences to NAs, convert all Oxygen saturation values to numeric
RespRate <- RespRate %>% 
  mutate_if(is.character,stringr::str_replace_all, pattern = "No Sensor",replacement = "NA") %>%
  mutate(across(where(is.character), as.numeric))
str(RespRate)

# add week to dataset based on day of experiment
RespRate$ID <- 1:length(RespRate$Block)
RespRate$Week <- NA

for(i in RespRate$ID){
  if(RespRate$Day[i] == 7 | RespRate$Day[i] == 8){
    RespRate$Week[i] <- 1
  }else if(RespRate$Day[i] == 14 | RespRate$Day[i] == 15){
    RespRate$Week[i] <- 2
  }else if(RespRate$Day[i] == 21 | RespRate$Day[i] == 22){
    RespRate$Week[i] <- 3
  }else if(RespRate$Day[i] == 28 | RespRate$Day[i] == 29){
    RespRate$Week[i] <- 4
  }else if(RespRate$Day[i] == 35 | RespRate$Day[i] == 36){
    RespRate$Week[i] <- 5
  }else{
    RespRate$Week[i] <- NA
  }
}


# remove the extra time over 140 minutes of sampling to standardize the amount of time sampled for each animal
        # some respiration trials ran over the 140 min, so we are removing the "extra" time so that regressions are computed over the same time window for all animals
RespRate <- filter(RespRate, Time.Min < 141)


# Samples from block 3, Day 22 plates 1-4 were samples on 1 min intervals instead of 2 min intervals
# Need to thin this data to match with all other plates in the data set
library(ds4psy)
RespRate_Block3Day22 <- RespRate %>% filter(Block == 3, Day == 22) %>% as.data.frame()
head(RespRate_Block3Day22)
dim(RespRate_Block3Day22)
range(RespRate_Block3Day22$ID)
RespRate_Block3Day22_thinned <- RespRate_Block3Day22 %>%
  filter(!is_wholenumber(ID/2))  # remove even rows so that timepoints match other plates (0, 2, 4 min, etc)
head(RespRate_Block3Day22_thinned)
tail(RespRate_Block3Day22_thinned)
dim(RespRate_Block3Day22_thinned)


# replace block 3, Day 22 plates 1-4 with thinned data
RespRate <- RespRate %>% filter(ID < 3433 | ID > 4017) # remove based on line IDs
dim(RespRate)  # should be 4242 rows (4806-564)
RespRate <- rbind(RespRate, RespRate_Block3Day22_thinned)
dim(RespRate)  # should be 4525 rows(4242+283)


# remove 3 outlier blanks from the data set 
RespRate[ RespRate$Block == 1 & RespRate$Day == 29 & RespRate$Plate == 2, "D1"] <- NA
RespRate[ RespRate$Block == 4 & RespRate$Day == 14 & RespRate$Plate == 4, "B6"] <- NA
RespRate[ RespRate$Block == 3 & RespRate$Day == 15 & RespRate$Plate == 3, "C5"] <- NA


# loop through respiration calculation for each block, day, plate, and well
    # compile into new data frame where each row represents a single block, day, plate, and well combination with the respiration rate calculated
    # join with the metadata to compare patterns across different treatments

# setting up vectors to cycle through
block <- c(1:4)
plate <- c(1:4)
wells <- colnames(select(RespRate, A1:D6))

# Set up empty columns to add data to in the for loop below
RespRateCalc$O2.sat.per.hr <- NA
RespRateCalc$mean.control.O2.sat <- NA
RespRateCalc$sd.control.O2.sat <- NA
RespRateCalc$VO2 <- NA
RespRateCalc$metabolic.rate <- NA


### Loop takes a VERY long time to run ###

for(i in block){
  thisblock <- filter(RespRate, Block == i) # filter data by each block
  if(i > 2){        # block 3 and 4 (Metsch expt) only have 3 weeks, while blocks 1 and 2 (Past) have 5 weeks
    week <- c(1:3)
  } else {
    week <- c(1:5)
  }
  for(j in week){
    thisweek <- filter(thisblock, Week == j) # filter block data by each week of the data set
    
    for(k in plate){
      thisplate <- filter(thisweek, Plate == k) # filter week data by each plate
      
      for(l in wells){
        if(sum(is.na(thisplate[ , l])) > 0){
          RespRateCalc[ RespRateCalc$Block == i & RespRateCalc$Week == j & RespRateCalc$Plate == k & RespRateCalc$Well == l, "O2.sat.per.hr"] <- NA
        } else{
          ox.sat.values <- unlist(c((thisplate[ , l])), use.names = F)
          # run all the local regressions to estimate the monotonic slope for each animal
          # use alpha threshold of at least 30% of the data included in the regression
          # Generally the L% method (pc) seems to fit the data best for most animals
          daphniaRegs  <-  rankLocReg(xall=thisplate$Time.Min, yall=ox.sat.values, alpha=0.3, # (ma = rate of change in oxygen saturation for an experimental treatment)
                                      method="pc", verbose=TRUE) 
          # get the specific row for the animal to store associated data
          RespRateCalc[ RespRateCalc$Block == i & RespRateCalc$Week == j & RespRateCalc$Plate == k & RespRateCalc$Well == l, "O2.sat.per.hr"] <- daphniaRegs$allRegs[1, "b1"] * 60  # extract the best monotonic slope (per min), multiply by 60 to get the per hour slope, store the Oxygen saturation per hour
        }

      }
      # filter data to get only this plate 
      temp_plate <-  filter(RespRateCalc, Block == i & Week == j & Plate == k)
      
      # filter data to get only blank controls 
      temp_controls <-  filter(temp_plate, Parasites == "Blank") 
      
      # calculated the mean of the wells with controls only (mc = per run average rate of change for the blank controls)
      control <- mean(temp_controls$O2.sat.per.hr)
      control.sd <- sd(temp_controls$O2.sat.per.hr)
      
      # Calculate the rate of oxygen consumption (VO2) (mL O2/hr)
      # VO2=−1×[(ma−mc)∕100]×V×𝛽O2
      for(m in wells){
        RespRateCalc[ RespRateCalc$Block == i & RespRateCalc$Week == j & RespRateCalc$Plate == k & RespRateCalc$Well == m, "VO2"] <- 
          -1 * ((RespRateCalc[ RespRateCalc$Block == i & RespRateCalc$Week == j & RespRateCalc$Plate == k & RespRateCalc$Well == m, "O2.sat.per.hr"] - control) / 100) * 0.002 * 6.15 
        RespRateCalc[ RespRateCalc$Block == i & RespRateCalc$Week == j & RespRateCalc$Plate == k & RespRateCalc$Well == m, "mean.control.O2.sat"] <- control
        RespRateCalc[ RespRateCalc$Block == i & RespRateCalc$Week == j & RespRateCalc$Plate == k & RespRateCalc$Well == m, "sd.control.O2.sat"] <- control.sd
      }
      # V = volume of water in vials in mL (we had 2 mL vials)  = 0.002 L
      # 𝛽O2 = oxygen capacitance of air-saturated water at 20°C = 6.40 (Cameron, 1986) in mg/L (I think)   
      
      ################## NOTE: my samples appear to be run with an internal temp of 21.75 to 22.5°C, 
      ################# so maybe need to change this number for my calculations to be more accurate?
              # 6.15 at 22°C (volume solubility of O2) mL O2 per L of water (Cameron, 1986) 
      
    }
  }
}
str(thisplate)

# calculate the metabolic rate (converted to metabolic rate (J/hr) using the calorific conversion factor of 20.08 J/ml O2 (Lighton, 2008))
RespRateCalc <- RespRateCalc %>%
  select(Block:Died.After.Resp, O2.sat.per.hr:metabolic.rate, Notes) %>%
  mutate(metabolic.rate = VO2 * 20.13)  # metabolic rate is J/hr (calculated as 16+5.164(RQ))
# where RQ, the respiratory quotient, is the ratio of VCO2 to VO2, which can be assumed to be 0.8 when not measured (Lighton 2008)
# Problem is that calculates to 20.1312. 
# to get 20.08, the RQ would need to 0.79
View(RespRateCalc)

# Set the Diet factor order
RespRateCalc$Diet <- factor(RespRateCalc$Diet, levels = c("S", "SM", "M", "M+"))

# save CSV file
write.csv(RespRateCalc, here("data/ResourceQuality_RespirationRateCalc.csv"), quote = F, row.names=FALSE)



# check the metabolic rate values for the animals that died during or within 24 hours of respiration trial
Dead_animals <- filter(RespRateCalc, Died.After.Resp == "Y")
      # 6 of the 19 animals marked in Past experiment have negative values for metabolic rate, most are low values
      # 3 of the 12 animals marked in the Metsch experiment have negative values for metabolic rate


##############################
# check the Oxygen saturation per hr slope calculation for a few extreme values, 
##############################

# Block 1 day 15 plate 1 S Blank B2 - positive rate
test1 <- RespRate %>% filter(Block == 1, Day == 15, Plate == 1) %>% as.data.frame()
daphniaRegs_test1  <-  rankLocReg(xall=test1$Time.Min, yall=test1$B2, alpha=0.3, 
                            method="pc", verbose=TRUE) 
summary(daphniaRegs_test1) 
plot(daphniaRegs_test1, rank = 1)       # L% is the best fit
outputRankLocRegPlot(daphniaRegs_test1) # yes, slight positive slope but not much change in O2 sat over time

# Block 1 day 8 plate 2 SM Blank A1 - positive rate
test2 <- RespRate %>% filter(Block == 1, Day == 8, Plate == 2) %>% as.data.frame()
daphniaRegs_test2  <-  rankLocReg(xall=test2$Time.Min, yall=test2$A1, alpha=0.3, 
                                  method="pc", verbose=TRUE) 
summary(daphniaRegs_test2)
plot(daphniaRegs_test2, rank = 1)       # L% is the best fit
outputRankLocRegPlot(daphniaRegs_test2) # yes, slight positive slope but not much change in O2 sat over time

# Block 3 day 15 plate 2 SM Uninf STD B2 - very negative rate
test3 <- RespRate %>% filter(Block == 3, Day == 15, Plate == 2) %>% as.data.frame()
daphniaRegs_test3  <-  rankLocReg(xall=test3$Time.Min, yall=test3$B2, alpha=0.3, 
                                  method="pc", verbose=TRUE) 
summary(daphniaRegs_test3)
plot(daphniaRegs_test3, rank = 1)       # L% is the best fit
outputRankLocRegPlot(daphniaRegs_test3) # yes, very negative slope but fits the data well, not an error

# Block 3 day 15 plate 3 M Uninf MID D4 - very negative rate
test4 <- RespRate %>% filter(Block == 3, Day == 15, Plate == 3) %>% as.data.frame()
daphniaRegs_test4  <-  rankLocReg(xall=test4$Time.Min, yall=test4$D4, alpha=0.3, 
                                  method="pc", verbose=TRUE) 
summary(daphniaRegs_test4)
plot(daphniaRegs_test4, rank = 1)       # L% is the best fit
outputRankLocRegPlot(daphniaRegs_test4) # yes, very negative slope but fits the data well, not an error

# Block 2 day 14 plate 1 S Uninf MID a4 - sample had an error in recording initially (in C2), but moved to new well (A4)
test5 <- RespRate %>% filter(Block == 2, Day == 14, Plate == 1) %>% filter(!is.na(A4)) %>% as.data.frame()
daphniaRegs_test5  <-  rankLocReg(xall=test5$Time.Min, yall=test5$A4, alpha=0.3, 
                                  method="pc", verbose=TRUE) 
summary(daphniaRegs_test5)
plot(daphniaRegs_test5, rank = 1)       # L% is the best fit
outputRankLocRegPlot(daphniaRegs_test5) # Very BAD fit compared to other animals in the same treatment, do not include this point
        # since there are NAs in part of the time series, the rankLocReg function in the loop does not calculate
        # a slope, so it is already an NA in the calculated metabolic rates

# Block 3 day 15 plate 3 M Blank C5 - check extreme blank control (-8.58 slope)
test4 <- RespRate %>% filter(Block == 3, Day == 15, Plate == 3) %>% as.data.frame()
daphniaRegs_test4a  <-  rankLocReg(xall=test4$Time.Min, yall=test4$C5, alpha=0.3, 
                                  method="pc", verbose=TRUE) 
summary(daphniaRegs_test4a)
plot(daphniaRegs_test4a, rank = 1)       # L% is the best fit
outputRankLocRegPlot(daphniaRegs_test4a) # yes, very negative slope but fits the data well, not an error
        # range of measurements is much lower than all other blanks, but similar to some other wells nearby  

# Block 4 day 14 plate 4 M+ Blank B6 - check extreme blank control (-5.82 slope)
test6 <- RespRate %>% filter(Block == 4, Day == 14, Plate == 4) %>% as.data.frame()
daphniaRegs_test6  <-  rankLocReg(xall=test6$Time.Min, yall=test6$B6, alpha=0.3, 
                                  method="pc", verbose=TRUE) 
summary(daphniaRegs_test6)
plot(daphniaRegs_test6, rank = 1)       # L% is the best fit
outputRankLocRegPlot(daphniaRegs_test6) # yes, very negative slope but fits the data well, not an error
# range of measurements is much lower than all other wells 

# Block 1 day 29 plate 2 SM Blank D1 - check extreme blank control (-4.68 slope)
test7 <- RespRate %>% filter(Block == 1, Day == 29, Plate == 2) %>% as.data.frame()
daphniaRegs_test7  <-  rankLocReg(xall=test7$Time.Min, yall=test7$D1, alpha=0.3, 
                                  method="pc", verbose=TRUE) 
summary(daphniaRegs_test7)
plot(daphniaRegs_test7, rank = 1)       # L% is the best fit
outputRankLocRegPlot(daphniaRegs_test7) # yes, very negative slope but fits the data well, not an error
# range of measurements is much lower than all other wells 



###############################
# Initial figures to explore the respiration data
##############################

# function to make a grid plot with a shared legend
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}


# diet color palette
diet_colors <- c("#ADDD8E", "#41AB5D", "#006837", "#1D91C0")

# filter day to remove NAs and by week 1
RespRateCalc_plot_Week1 <- RespRateCalc %>%
  filter(!is.na(Parasites)) %>%
  filter(Week == 1) 

plot_week1 <- ggplot(data = RespRateCalc_plot_Week1, aes(x = Parasites, y = metabolic.rate)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width=0.8)) +
  geom_point(aes(color=Diet), size=2, position=position_jitterdodge(dodge.width=0.8, jitter.width = 0.05), alpha = 0.4) +
  facet_grid(~Clone) +
  scale_x_discrete(limits = c("Uninf", "Metsch", "Pasteuria", "Blank")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  ggtitle("Test plot of Respiration Rate parasite exposure x diet x host clone") +
  labs(x = "Parasite Exposure", y = "Metabolic Rate (J/hr)") +
  theme_classic()
plot_week1



plot_week1_oxysat <- ggplot(data = RespRateCalc_plot_Week1, aes(x = Parasites, y = O2.sat.per.hr)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width=0.8)) +
  geom_point(aes(color=Diet), size=2, position=position_jitterdodge(dodge.width=0.8, jitter.width = 0.05), alpha = 0.4) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninf", "Metsch", "Pasteuria", "Blank")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  ggtitle("Test plot of O2 saturation per hour parasite exposure x diet x host clone") +
  labs(x = "Parasite Exposure", y = "Oxygen Saturation Per hour (%/hr)") +
  theme_classic()
plot_week1_oxysat

View(RespRateCalc_plot_Week1)

# filter day to remove NAs and by week 2
RespRateCalc_plot_Week2 <- RespRateCalc %>%
  filter(!is.na(Parasites)) %>%
  filter(Week == 2) 
str(RespRateCalc_plot_Week2)

plot_week2 <- ggplot(data = RespRateCalc_plot_Week2, aes(x = Parasites, y = metabolic.rate)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width=0.8)) +
  geom_point(aes(color=Diet), size=2, position=position_jitterdodge(dodge.width=0.8, jitter.width = 0.05), alpha = 0.4) +
  facet_wrap(~Clone) +
  #scale_x_discrete(limits = c("Uninf", "Metsch", "Pasteuria", "Blank")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  ggtitle("Test plot of Respiration Rate parasite exposure x diet x host clone") +
  labs(x = "Parasite Exposure", y = "Metabolic Rate (J/hr)") +
  theme_classic()
plot_week2



plot_week2_oxysat <- ggplot(data = RespRateCalc_plot_Week2, aes(x = Parasites, y = O2.sat.per.hr)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width=0.8)) +
  geom_point(aes(color=Diet), size=2, position=position_jitterdodge(dodge.width=0.8, jitter.width = 0.05), alpha = 0.4) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninf", "Metsch", "Pasteuria", "Blank")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  ggtitle("Test plot of O2 saturation per hour parasite exposure x diet x host clone") +
  labs(x = "Parasite Exposure", y = "Oxygen Saturation Per hour (%/hr)") +
  theme_classic()
plot_week2_oxysat


# filter day to remove NAs and by week 3
RespRateCalc_plot_Week3 <- RespRateCalc %>%
  filter(!is.na(Parasites)) %>%
  filter(Week == 3)  %>%
  filter(Died.After.Resp != "Y")
str(RespRateCalc_plot_Week3)
sum(is.na(RespRateCalc_plot_Week3$O2.sat.per.hr))

plot_week3 <- ggplot(data = RespRateCalc_plot_Week3, aes(x = Parasites, y = metabolic.rate)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width=0.8)) +
  geom_point(aes(color=Diet), size=2, position=position_jitterdodge(dodge.width=0.8, jitter.width = 0.05), alpha = 0.4) +
  facet_wrap(~Clone) +
  #scale_x_discrete(limits = c("Uninf", "Metsch", "Pasteuria", "Blank")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  ggtitle("Test plot of Respiration Rate parasite exposure x diet x host clone") +
  labs(x = "Parasite Exposure", y = "Metabolic Rate (J/hr)") +
  theme_classic()
plot_week3



plot_week3_oxysat <- ggplot(data = RespRateCalc_plot_Week3, aes(x = Parasites, y = O2.sat.per.hr)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width=0.8)) +
  geom_point(aes(color=Diet), size=2, position=position_jitterdodge(dodge.width=0.8, jitter.width = 0.05), alpha = 0.4) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninf", "Metsch", "Pasteuria", "Blank")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  ggtitle("Test plot of O2 saturation per hour parasite exposure x diet x host clone") +
  labs(x = "Parasite Exposure", y = "Oxygen Saturation Per hour (%/hr)") +
  theme_classic()
plot_week3_oxysat


## Plot of O2 saturation and metabolism over time by week
RespRateCalc_Sum.clean <- RespRateCalc %>%
  filter(!is.na(Parasites), Parasites != "Blank", Died.After.Resp == "N") %>%
  group_by(Diet, Parasites, Clone, Week, Treatment, Experiment) %>%
  summarize(N  = length(O2.sat.per.hr),
            O2.mean = mean(O2.sat.per.hr, na.rm = T),
            O2.sd   = sd(O2.sat.per.hr, na.rm = T),
            O2.se   = O2.sd / sqrt(N),
            resp.mean = mean(metabolic.rate, na.rm = T),
            resp.sd   = sd(metabolic.rate, na.rm = T),
            resp.se   = resp.sd / sqrt(N))
RespRateCalc_Sum.clean$Week <- as.factor(RespRateCalc_Sum.clean$Week)
RespRateCalc_Sum.clean$Diet <- factor(RespRateCalc_Sum.clean$Diet, levels = c("S", "SM", "M", "M+"))


#  Grouped figures (both experiments)
plot_byWeek <- ggplot(data = RespRateCalc_Sum.clean, aes(x = Week, y = O2.mean, group = Treatment, shape = Diet)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_errorbar(aes(x=Week, ymin=O2.mean-O2.se, ymax=O2.mean+O2.se,color=Diet), width=0.2) + 
  geom_line(aes(color=Diet, group = Treatment), linewidth=2) +
  geom_point(aes(color=Diet), size=4, alpha = 0.8) +
  facet_grid(Clone~Parasites) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  ggtitle("Test plot of O2 saturation per hour by week and parasite x diet x host clone") +
  labs(x = "Week of Experiment", y = "Oxygen Saturation Per hour (%/hr)") +
  theme_classic()
plot_byWeek

plot_byWeek_resp <- ggplot(data = RespRateCalc_Sum.clean, aes(x = Week, y = resp.mean, group = Treatment, shape = Diet)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_errorbar(aes(x=Week, ymin=resp.mean-resp.se, ymax=resp.mean+resp.se,color=Diet), width=0.2) + 
  geom_line(aes(color=Diet, group = Treatment), linewidth=2) +
  geom_point(aes(color=Diet), size=4, alpha = 0.8) +
  facet_grid(Clone~Parasites) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  ggtitle("Test plot of metabolic rate per hour by week and parasite x diet x host clone") +
  labs(x = "Week of Experiment", y = "Metabolic Rate (J/hr)") +
  theme_classic()
plot_byWeek_resp

# split figures by past and metsch experiments
RespRateCalc_Sum.Past <- filter(RespRateCalc_Sum.clean, Experiment == "Past")
RespRateCalc_Sum.Metsch <- filter(RespRateCalc_Sum.clean, Experiment == "Metsch")

plot_byWeek_resp_past <- ggplot(data = RespRateCalc_Sum.Past, aes(x = Week, y = resp.mean, group = Treatment, shape = Diet)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_errorbar(aes(x=Week, ymin=resp.mean-resp.se, ymax=resp.mean+resp.se,color=Diet), width=0.2) + 
  geom_line(aes(color=Diet, group = Treatment), linewidth=1.5) +
  geom_point(aes(color=Diet), size=3, alpha = 0.8) +
  facet_grid(Clone~Parasites) +
  ylim(-0.001, 0.009) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  annotate("rect", xmin = 0.9, xmax = 1.1, ymin = -0.001, ymax = 0.009,
           alpha = .2,fill = "gray") +
  ggtitle("Past Experiment Metabolic Rate") +
  labs(x = "Week of Experiment", y = "Metabolic Rate (J/hr)") +
  theme_classic()
plot_byWeek_resp_past
ggsave(here("figures/PastMetabolism_DietxCloneXParasite_NoBlankOutliers_NoDead.tiff"), plot = plot_byWeek_resp_past, dpi = 300, width = 6, height = 6, units = "in", compression="lzw")


plot_byWeek_resp_metsch <- ggplot(data = RespRateCalc_Sum.Metsch, aes(x = Week, y = resp.mean, group = Treatment, shape = Diet)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_errorbar(aes(x=Week, ymin=resp.mean-resp.se, ymax=resp.mean+resp.se,color=Diet), width=0.2) + 
  geom_line(aes(color=Diet, group = Treatment), linewidth=1.5) +
  geom_point(aes(color=Diet), size=3, alpha = 0.8) +
  facet_grid(Clone~Parasites) +
  ylim(-0.001, 0.009) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  annotate("rect", xmin = 0.9, xmax = 1.1, ymin = -0.001, ymax = 0.009,
           alpha = .2,fill = "gray") +
  ggtitle("Metsch Experiment Metabolic rate") +
  labs(x = "Week of Experiment", y = "Metabolic Rate (J/hr)") +
  theme_classic()
plot_byWeek_resp_metsch
ggsave(here("figures/MetschMetabolism_DietxCloneXParasite_NoBlankOutliers_NoDead.tiff"), plot = plot_byWeek_resp_metsch, dpi = 300, width = 6, height = 6, units = "in", compression="lzw")


Metabolism_fig <- grid_arrange_shared_legend(plot_byWeek_resp_past, plot_byWeek_resp_metsch, nrow=1, ncol = 2, position = "right")
ggsave(here("figures/Metabolism_byExpt_NoBlankOutliers_NoDead.tiff"), plot = Metabolism_fig, dpi = 300, width = 9, height = 6, units = "in", compression="lzw")


plot_byWeek_O2_past <- ggplot(data = RespRateCalc_Sum.Past, aes(x = Week, y = O2.mean, group = Treatment, shape = Diet)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_errorbar(aes(x=Week, ymin=O2.mean-O2.se, ymax=O2.mean+O2.se,color=Diet), width=0.2) + 
  geom_line(aes(color=Diet, group = Treatment), linewidth=1.5) +
  geom_point(aes(color=Diet), size=3, alpha = 0.8) +
  facet_grid(Clone~Parasites) +
  ylim(-6.5, 0) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  ggtitle("Past Experiment O2 Saturation Rate") +
  labs(x = "Week of Experiment", y = "Oxygen Saturation Per hour (%/hr)") +
  theme_classic()
plot_byWeek_O2_past
ggsave(here("figures/PastO2Sat_DietxCloneXParasite_NoBlankOutliers.tiff"), plot = plot_byWeek_O2_past, dpi = 300, width = 6, height = 6, units = "in", compression="lzw")


plot_byWeek_O2_metsch <- ggplot(data = RespRateCalc_Sum.Metsch, aes(x = Week, y = O2.mean, group = Treatment, shape = Diet)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  geom_errorbar(aes(x=Week, ymin=O2.mean-O2.se, ymax=O2.mean+O2.se,color=Diet), width=0.2) + 
  geom_line(aes(color=Diet, group = Treatment), linewidth=1.5) +
  geom_point(aes(color=Diet), size=3, alpha = 0.8) +
  facet_grid(Clone~Parasites) +
  ylim(-6.5, 0) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  ggtitle("Metsch Experiment O2 Saturation rate") +
  labs(x = "Week of Experiment", y = "Oxygen Saturation Per hour (%/hr)") +
  theme_classic()
plot_byWeek_O2_metsch
ggsave(here("figures/MetschO2Sat_DietxCloneXParasite_NoBlankOutliers.tiff"), plot = plot_byWeek_O2_metsch, dpi = 300, width = 6, height = 6, units = "in", compression="lzw")


O2_fig <- grid_arrange_shared_legend(plot_byWeek_O2_past, plot_byWeek_O2_metsch, nrow=1, ncol = 2, position = "right")
ggsave(here("figures/O2_Sat_byExpt_NoBlankOutliers.tiff"), plot = O2_fig, dpi = 300, width = 9, height = 6, units = "in", compression="lzw")



## check of variation in blank controls within and across plates
Blank_data <- RespRateCalc %>%
  filter(!is.na(Parasites), Parasites == "Blank") %>%
  mutate(block_plate = paste0(Block, "_", Plate))
Blank_data$Plate <- as.factor(Blank_data$Plate)
Blank_data$Week <- as.factor(Blank_data$Week)
Blank_data$Diet <- factor(Blank_data$Diet, levels = c("S", "SM", "M", "M+"))

blankplot_byWeek <- ggplot(data = Blank_data, aes(x = Week, y = O2.sat.per.hr, group = Diet, shape = Diet)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  #geom_errorbar(aes(x=Week, ymin=O2.mean-O2.se, ymax=O2.mean+O2.se,color=Diet), width=0.2) + 
  #geom_line(aes(color=Diet, group = Diet), linewidth=2) +
  geom_jitter(aes(color=Diet), size=3, alpha = 0.8, width = 0.2) +
  facet_wrap(~Block) +
  #scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  ggtitle("Test plot of CONTROLS O2 saturation per hour by week and diet") +
  labs(x = "Week of Experiment", y = "Oxygen Saturation Per hour (%/hr)") +
  theme_classic()
blankplot_byWeek
ggsave(here("figures/RespBlanks_byWeek_Block_plate.tiff"), plot = blankplot_byWeek, dpi = 300, width = 9, height = 6, units = "in", compression="lzw")

# read in the full data set without anything removed
data <- read.csv(here("data/ResourceQuality_RespirationRateCalc.csv"))
data_blank <- filter(data, Parasites == "Blank")

# Split the blank data by block
Blank_data_block1 <- filter(data_blank, Block == "1")
Blank_data_block2 <- filter(data_blank, Block == "2")
Blank_data_block3 <- filter(data_blank, Block == "3")
Blank_data_block4 <- filter(data_blank, Block == "4")

# test for outliers in each block
# block 1
outlier_test_block1 <- rosnerTest(Blank_data_block1$O2.sat.per.hr, k = 5)
outlier_test_block1$all.stats
      # one outlier = -4.7015417
Blank_data_block1[ Blank_data_block1$O2.sat.per.hr < -4, ]

# block 2
outlier_test_block2 <- rosnerTest(Blank_data_block2$O2.sat.per.hr, k = 5)
outlier_test_block2$all.stats
      # no outliers

# block 3
outlier_test_block3 <- rosnerTest(Blank_data_block3$O2.sat.per.hr, k = 5)
outlier_test_block3$all.stats
      # one outlier = -8.621880
Blank_data_block3[ Blank_data_block3$O2.sat.per.hr < -8, ]

# block 4
outlier_test_block4 <- rosnerTest(Blank_data_block4$O2.sat.per.hr, k = 5)
outlier_test_block4$all.stats
# one outlier = -5.851114 
Blank_data_block4[ Blank_data_block4$O2.sat.per.hr < -5, ]



# Split the blank data by experiment
Blank_data_past <- filter(Blank_data, Experiment == "Past")
Blank_data_metsch <- filter(Blank_data, Experiment == "Metsch")

# test for outliers in the Past blank data
outlier_test_past <- rosnerTest(Blank_data_past$O2.sat.per.hr, k = 5)
outlier_test_past$all.stats
          # one outlier

# test for outliers in the Metsch blank data
outlier_test_metsch <- rosnerTest(Blank_data_metsch$O2.sat.per.hr, k = 5)
outlier_test_metsch$all.stats
        # one outlier

# test for outliers in the whole blank data set
outlier_test <- rosnerTest(Blank_data$O2.sat.per.hr, k = 5)
outlier_test$all.stats
      # three values are considered outliers in the blank O2 saturation data

# filter out outliers from blank data
Blank_data_NoOutliers <- filter(Blank_data, O2.sat.per.hr > -5.2)

# identify the outlier points to remove ahead of the loop above
Blank_data_Outliers <- filter(Blank_data, O2.sat.per.hr < -5.2)

# model of all blank data - sig diff in Diet*Week and Block*Week
mod <- aov(O2.sat.per.hr ~ Diet*Block*Week, data = Blank_data)
summary(mod)

# model of blank data without outliers - all interactions are significant now...
mod2 <- aov(O2.sat.per.hr ~ Diet*Block*Week, data = Blank_data_NoOutliers)
summary(mod2)

# plot of blank data without outliers
blankplot_byWeek2 <- ggplot(data = Blank_data_NoOutliers, aes(x = Week, y = O2.sat.per.hr, group = Diet, shape = Diet)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray") +
  #geom_errorbar(aes(x=Week, ymin=O2.mean-O2.se, ymax=O2.mean+O2.se,color=Diet), width=0.2) + 
  #geom_line(aes(color=Diet, group = Diet), linewidth=2) +
  geom_jitter(aes(color=Diet), size=3, alpha = 0.8, width = 0.2) +
  facet_wrap(~Block) +
  #scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_color_manual(values = diet_colors) +
  scale_shape_manual(values=c(15,16,17,18)) +
  ggtitle("Plot of CONTROLS O2 saturation/hr by week and diet, No Outliers") +
  labs(x = "Week of Experiment", y = "Oxygen Saturation Per hour (%/hr)") +
  theme_classic()
blankplot_byWeek2
ggsave(here("figures/RespBlanks_byWeek_Block_plate_NoOutliers.tiff"), plot = blankplot_byWeek2, dpi = 300, width = 9, height = 6, units = "in", compression="lzw")


a <- emmeans(mod, specs = pairwise ~ Block | Week, type = "response")

mod2 <- aov(O2.sat.per.hr ~ Block*Week, data = Blank_data)
summary(mod2)
