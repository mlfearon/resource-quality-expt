# Processing Feeding rate and Respiration Rate data files

#Michelle's working directory


library(ggplot2)
library(dplyr)
library(tidyverse)
library(lubridate)
library(fs)
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

View(FeedingRate)

# recode block, day and plate to numbers
FeedingRate$Block <- recode(FeedingRate$Block, Block1 = 1, Block2 = 2, Block3 = 3, Block4 = 4)
FeedingRate$Day <- recode(FeedingRate$Day, Day7 = 7, Day8 = 8, Day14 = 14, Day15 = 15, Day21 = 21, Day22 = 22, Day28 = 28, Day29 = 29, Day35 = 35, Day36 = 36)
FeedingRate$Plate <- recode(FeedingRate$Plate, Plate1 = 1, Plate2 = 2, Plate3 = 3, Plate4 = 4, Plate5 = 5, Plate6 = 6, Plate7 = 7, Plate8 = 8)


# need to join meta data to the feeding rate data
setwd("C:/Users/mlfearon/OneDrive - Umich/PROJECTS/MHMP Daphnia Duffy/Resource Quality/Data and Code")


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


# filter out fluorescence readings that don't have a sample associated with it
FeedingRate <- na.omit(FeedingRate)




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


# remove outliers
dim(FeedingRate)
#FeedingRate2 <- FeedingRate[FeedingRate$Block != 1 & FeedingRate$Day != 21 & FeedingRate$Plate != 1 & FeedingRate$Diet != "S" & FeedingRate$Well.x != "H5"]

FeedingRate <- FeedingRate[-5913, ] # not elegant but it works to remove the outlier...

FeedingRate <- FeedingRate %>%
  filter(chlorophyll < 70)

dim(FeedingRate) # success!
boxplot(FeedingRate$chlorophyll ~ FeedingRate$Diet)


# get data sets with NoDaphnia and No algae controls
FeedingRate_NoAlgae <- filter(FeedingRate, Clone == "NoAlgae")
FeedingRate_NoDaphnia <- filter(FeedingRate, Clone == "NoDaphnia")

# look for outliers
boxplot(FeedingRate_NoDaphnia$chlorophyll ~ FeedingRate_NoDaphnia$Diet)
boxplot(FeedingRate_NoAlgae$chlorophyll ~ FeedingRate_NoAlgae$Diet)


# remove daphnia that died during feeding rate trial
FeedingRate_died <- filter(FeedingRate, Notes == "died")
View(FeedingRate_died)
FeedingRate <- filter(FeedingRate, Notes != "died")
dim(FeedingRate)


# Summarize data by block, day, plate, column and date/time
    # calculate a mean for the chlorophyll reading based on the groups so that 
    # all 8 fluorescence reading replicates are averaged together

FeedingRateCalc <- FeedingRate %>%
  group_by(Block, Day, Plate, Column, Time, Diet, Parasites, Clone, Rep) %>%
  summarize(chlorophyll.avg = mean(chlorophyll))

FeedingRateCalc_NoDaphnia <- filter(FeedingRateCalc, Clone == "NoDaphnia")
FeedingRateCalc_NoAlgae <- filter(FeedingRateCalc, Clone == "NoAlgae")
FeedingRateCalc_Treatment <- filter(FeedingRateCalc, Clone != "NoDaphnia", Clone != "NoAlgae")


View(FeedingRateCalc_Treatment)
class(FeedingRateCalc_Treatment$Day)
FeedingRateCalc_Treatment$Rep <- as.integer(FeedingRateCalc_Treatment$Rep)


# add in bodysize for each daphnia
### bodysize data
bodysize_past <- read.csv("ResourceQuality_Pasteuria_Bodysize.csv", stringsAsFactors = F)
head(bodysize_past)
bodysize_metsch <- read.csv("ResourceQuality_Metsch_Bodysize.csv", stringsAsFactors = F)
head(bodysize_metsch)



# create a long version of the data set based on experimental day
bodysize_past_long <- bodysize_past %>%
  pivot_longer(length_7:length_38, names_to = "Day", names_prefix = "length_", values_to = "Length") %>%
  filter(!is.na(Length), Day != 38) %>%  # remove days when there are no bodysize measurements, remove last date where there are no corresponding feeding rate
  select(!Unique.code)
View(bodysize_past_long)


bodysize_metsch_long <- bodysize_metsch %>%
  pivot_longer(length_7:length_24, names_to = "Day", names_prefix = "length_", values_to = "Length") %>%
  filter(!is.na(Length), Day != 24) %>%  # remove days when there are no bodysize measurements, remove last date where there are no corresponding feeding rate
  select(!Unique.code)

bodysize_long <- bind_rows(bodysize_past_long, bodysize_metsch_long)
 

# recode Mid37 and Std labels to match feeding rate
bodysize_long$Clone <- recode(bodysize_long$Clone, Mid37 = "MID", Standard = "STD")
bodysize_long$Parasites <- recode(bodysize_long$Parasites, Uninfected = "Uninf")
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
  dplyr::rename(chlorophyll.nodaphnia = "chlorophyll.avg") %>%
  ungroup() %>%
  select(-Column, -Clone, -Rep)

FeedingRateCalc_Treatment <- full_join(FeedingRateCalc_Treatment, FeedingRateCalc_NoDaphnia, by = c("Block", "Day", "Plate", "Time", "Diet"))


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


# Calculate feeding rate (per ml, per hr)
    # ln (mean no daphnia control/ mean remaining food in sample) * (Volume 10 mL/ length of time of the assay in hr)


# add in total volume of feeding vials (10 mL)
FeedingRateCalc_Treatment$volume <- 10

# (ln (mean no daphnia control) - ln (mean remaining food in sample)) * (Volume 10 mL/ length of time of the assay in hr)
FeedingRateCalc_Treatment$Clearance <- (log(FeedingRateCalc_Treatment$chlorophyll.nodaphnia) - log(FeedingRateCalc_Treatment$chlorophyll.avg))* (FeedingRateCalc_Treatment$volume / FeedingRateCalc_Treatment$Time)

# calculate the feeding rate for relative bodysize (mL -hr -mm)
FeedingRateCalc_Treatment$Clearance_rel <- FeedingRateCalc_Treatment$Clearance / ((FeedingRateCalc_Treatment$Length/1000) ^ 2)

ggplot(data = FeedingRateCalc_Treatment, aes( x= Diet, y = Clearance, color = Parasites.x)) +
  geom_boxplot() +
  geom_jitter() +
  facet_grid(Clone~Week) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_hline(yintercept = -0.4, color = "red", linetype = "dashed") 
boxplot(FeedingRateCalc_Treatment$Clearance ~ FeedingRateCalc_Treatment$Diet)





write.csv(FeedingRateCalc_Treatment, "data/ResourceQuality_FeedingRateCalc.csv", quote = F, row.names=FALSE)



################################################
# Respiration Rate
################################################


### Need to convert all Oxygen files from excel to csv, remove first 12 rows and last 7 columns (27-32).
## I tried to do this in a fancy way, and I was able to select just the oxygen files, but ran into problems with the 
## excel file type and some of the weird symbols that were in the header of the excel file prevented it from downloading
## Anyways I got frustrated and just did it by hand in the end... try again another day
      # code to select just the Oxygen docs, problem is with reading in the docs as they are due to special symbols in the header

# returns the file names as a character vector
file_paths2 <- fs::dir_ls(path = "./data/RespirationRate/Oxygen_Sat/", regexp = ('(Oxygen)(.csv)$')) # need to select just the oxygen docs
length(file_paths2)

df <- read.csv(here("data/RespirationRate/Oxygen_sat/Block1_Day22_plate3_results_Oxygen.csv"), stringsAsFactors = F,)
dim(df)
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

View(RespRate)

# recode block, day and plate to numbers
RespRate$Block <- recode(RespRate$Block, Block1 = 1, Block2 = 2, Block3 = 3, Block4 = 4)
RespRate$Day <- recode(RespRate$Day, Day7 = 7, Day8 = 8, Day14 = 14, Day15 = 15, Day21 = 21, Day22 = 22, Day28 = 28, Day29 = 29, Day35 = 35, Day36 = 36)
RespRate$Plate <- recode(RespRate$Plate, plate1 = 1, plate2 = 2, plate3 = 3, plate4 = 4)

# convert Date.Time column to a date and time format
RespRate$Date.Time <- as.POSIXct(RespRate$Date.Time, tz = "EST", format = "%d.%m.%Y %H:%M%OS") # adds in today's date, but only need the start/end times relative to each other


# convert all "No Sensor" occurrences to NAs, convert all Oxygen saturation values to numeric
RespRate <- RespRate %>% 
  mutate_if(is.character,stringr::str_replace_all, pattern = "No Sensor",replacement = "NA") %>%
  mutate(across(where(is.character), as.numeric))
str(RespRate)


# loop through respiration calculation for each block, day, plate, and well
    # compile into new data frame where each row represents a single block, day, plate, and well combination with the respiration rate calculated
    # join with the metadata to compare patterns across different treatments










# directory where csv files are located
path <- file.path(getwd())


#-------------------
# make a list of all file names with csv or txt ext.
#-------------------
# $         : end of file name
# (csv|txt) : multiple file extentions
# \\.       : avoid unwanted cases such as .ccsv
#-------------------
v.filename <- list.files(path, pattern="\\.(csv|txt)$", 
                         ignore.case = TRUE, 
                         full.names = FALSE)



# read and make each data.frame
for(fn in v.filename) {
  
  df.each = read.csv(fn)
  
  print(fn); print(df.each)
}




df <- list.files(full.names = TRUE) %>% 
  lapply( function(x) {read_csv(x, skip=1 )}) %>%
  mutate(Name = str_remove_all(x, ".csv", "")) %>%
  bind_rows

??str_remove_all
?map_df

map_df(dat_files, ~read_csv(.x) %>%
         mutate(month_year = str_remove_all(.x, ".csv", "")) %>%
         separate(month_year, into=c("Month", "Year"), sep=" ")
