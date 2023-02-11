# Calculations of metabolic rate 

# methods and code were based on  Nørgaard et al 2020. Energetic scaling across different host densities and its consequences for pathogen proliferation. Functional Ecology
# DOI: 10.1111/1365-2435.13721 

https://colin-olito.github.io/LoLinR/vignettes/LoLinR.html

library(tidyverse)
library(here)
library(devtools)
install_github('colin-olito/LoLinR')
library(LoLinR)

# set the path to the script relative to the project root directory
here::i_am("scripts/respiration_rate_calculation.R")

# read in respiration file for a single plate
resp <- read.csv(here("data/Day10_plate3_results_Oxygen_example.csv"), stringsAsFactors = F, header = T)


# remove temp, pressure, salinity metadata (constant during sampling)
resp <- resp[,1:26] 

# remove the first 20 minutes of sampling while animals are accimating
resp <- resp[which(resp$Time.Min.>20),] 

resp.rate <- matrix(0,24,3)

# calculate the slope of oxygen over time for each well (ma = rate of change in oxygen saturation for an experimental treatment)
for(i in 3:length(resp)){
  resp.rate[i-2,2] <- lm(data = resp, resp[,i] ~ Time.Min.)[1]$coefficients[2] * 60 # multiply by 60 to convert rate from per minute to per hour
}

# calculated the mean of the wells with controls only (mc = per run average rate of change for the blank controls)
control <- mean(c(resp.rate[1,2],resp.rate[6,2],resp.rate[11,2],resp.rate[16,2]))

# Calculate the rate of oxygen consumption (VO2) (mL O2/hr)

# VO2=−1×[(ma−mc)∕100]×V×𝛽O2

for(i in 1:length(resp.rate[,2])){
  resp.rate[i,3] <- -1 * ((resp.rate[i,2] - control) / 100) * 2 * 6.4 
}
# V = volume of water in vials in mL (we had 2 mL vials)
# 𝛽O2 = oxygen capacitance of air-saturated water at 20°C = 6.40 (Cameron, 1986)

resp.rate <- as.data.frame(resp.rate)
resp.rate[,1] <- colnames(resp)[3:26]
colnames(resp.rate) <- c("well", "O2.sat.hr", "VO2")

# calculate the metabolic rate (converted to metabolic rate (J/hr) using the calorific conversion factor of 20.08 J/ml O2 (Lighton, 2008))
resp.rate <- resp.rate %>%
  mutate(metabolic.rate = VO2 * 20.08) # metabolic rate is J/hr


# try again with one I know the samples to see what that looks like.
# also need to fix that the rate of oxygen consumption should be per hour, but currently is per minute?





# test of calculations with LoLinR package to estimate the monotonic slope for each animal


daphniaRegs  <-  rankLocReg(xall=resp$Time.Min., yall=resp$A2, alpha=0.2, 
                           method="pc", verbose=TRUE) 
summary(daphniaRegs)
plot(daphniaRegs, rank = 1)
outputRankLocRegPlot(daphniaRegs)
