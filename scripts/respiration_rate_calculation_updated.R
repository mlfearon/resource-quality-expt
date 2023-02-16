# Calculations of metabolic rate 

# methods and code were based on  Nørgaard et al 2020. Energetic scaling across different host densities and its consequences for pathogen proliferation. Functional Ecology
# DOI: 10.1111/1365-2435.13721 

# https://colin-olito.github.io/LoLinR/vignettes/LoLinR.html

library(tidyverse)
library(here)
library(devtools)
install_github('colin-olito/LoLinR')
library(LoLinR)

# set the path to the script relative to the project root directory
here::i_am("scripts/respiration_rate_calculation.R")

# read in respiration file for a single plate
resp <- read.csv(here("data/Day10_plate3_results_Oxygen_example.csv"), stringsAsFactors = F, header = T)
# try with one of my plates
resp <- read.csv(here("data/RespirationRate/Oxygen_Sat/Block2_Day7_plate1_results_Oxygen.csv"))

# remove temp, pressure, salinity metadata (constant during sampling)
resp <- resp[,1:26] 

# remove the first 20 minutes of sampling while animals are acclimating
resp <- resp[which(resp$Time.Min.>20),] 

resp.rate <- matrix(0,24,3)

# calculate the slope of oxygen over time for each well (ma = rate of change in oxygen saturation for an experimental treatment)
for(i in 3:length(resp)){
  resp.rate[i-2,2] <- lm(data = resp, resp[,i] ~ Time.Min.)[1]$coefficients[2] * 60 # multiply by 60 to convert rate from per minute to per hour
}

# calculated the mean of the wells with controls only (mc = per run average rate of change for the blank controls)
control <- mean(c(resp.rate[1,2],resp.rate[6,2],resp.rate[11,2],resp.rate[16,2]))
control <- mean(c(resp.rate[17,2],resp.rate[6,2],resp.rate[23,2],resp.rate[16,2]))

# Calculate the rate of oxygen consumption (VO2) (mL O2/hr)

# VO2=−1×[(ma−mc)∕100]×V×𝛽O2

for(i in 1:length(resp.rate[,2])){
  resp.rate[i,3] <- -1 * ((resp.rate[i,2] - control) / 100) * 2 * 6.4 
}
# V = volume of water in vials in mL (we had 2 mL vials)
# 𝛽O2 = oxygen capacitance of air-saturated water at 20°C = 6.40 (Cameron, 1986)   

################## NOTE: my samples appear to be run with an internal temp of 21.75 to 22.5°C

resp.rate <- as.data.frame(resp.rate)
resp.rate[,1] <- colnames(resp)[3:26]
colnames(resp.rate) <- c("well", "O2.sat.hr", "VO2")

# calculate the metabolic rate (converted to metabolic rate (J/hr) using the calorific conversion factor of 20.08 J/ml O2 (Lighton, 2008))
resp.rate <- resp.rate %>%
  mutate(metabolic.rate = VO2 * 20.08) # metabolic rate is J/hr


# try again with one I know the samples to see what that looks like.
# also need to fix that the rate of oxygen consumption should be per hour, but currently is per minute?
meta_example <- read.csv(here("data/Block2_Day7_plate1_metadata_example.csv"))

resp.rate <- full_join(resp.rate, meta_example, by = c("well" = "Well"))

resp.rate_plot <- resp.rate %>%
  filter(Parasites != "Blank" & !is.na(Parasites))
diet_colors <- c("#ADDD8E", "#41AB5D", "#006837", "#1D91C0")
plot_test <- ggplot(data = resp.rate_plot, aes(x = Parasites, y = metabolic.rate)) +
  geom_boxplot(aes(fill=Diet),position=position_dodge(width=0.8)) +
  geom_point(aes(color=Diet), size=2, position=position_jitterdodge(dodge.width=0.8, jitter.width = 0.25), alpha = 0.4) +
  facet_wrap(~Clone) +
  scale_x_discrete(limits = c("Uninf", "Pasteuria")) +
  scale_fill_discrete(limits = c("S", "SM", "M", "M+")) +
  scale_fill_manual(values = diet_colors) +
  scale_color_manual(values = c(rep("black", 4))) +
  ggtitle("Test plot of Respiration Rate parasite exposure x diet x host clone for Past Expt") +
  labs(x = "Parasite Exposure", y = "Metabolic Rate (J/hr)") +
  theme_classic()
plot_test




# test of calculations with LoLinR package to estimate the monotonic slope for each animal


daphniaRegs  <-  rankLocReg(xall=resp$Time.Min., yall=resp$B1, alpha=0.2, 
                           method="pc", verbose=TRUE) 
summary(daphniaRegs)
plot(daphniaRegs, rank = 1)
outputRankLocRegPlot(daphniaRegs)
  ## Need to figure out how to determine which regression to use, are the ones in the table ranked and how?
  ## plan to loop through this code, extract the best regression slope to get the change in oxygen saturation per min
