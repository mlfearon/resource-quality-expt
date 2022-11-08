library(here)

resp <- read.csv("./respiration_metsch/day10_plate3_results_Oxygen.csv",header = T)

resp <- resp[,1:26]

resp <- resp[which(resp$Time.Min.>20),]

resp.rate <- matrix(0,24,2)

for(i in 3:length(resp)){
  resp.rate[i-2,1] <- lm(data = resp, resp[,i] ~ Time.Min.)[1]$coefficients[2]
}


control <- mean(c(resp.rate[1,1],resp.rate[6,1],resp.rate[11,1],resp.rate[16,1]))

for(i in 1:length(resp.rate[,2])){
  resp.rate[i,2] <- -1 * ((resp.rate[i,1] - control) / 100) * 0.002 * 6.4
}


