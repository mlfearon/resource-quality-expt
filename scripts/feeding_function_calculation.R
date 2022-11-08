library(here)
library(R2jags)
library(fitdistrplus)
library(mcmcplots)

# Set working directory
setwd("C:/Users/mlfea/OneDrive/Documents/Projects/MHMP Daphnia Duffy/Resource Quality/Data and Code")

RFU <- read.csv("feeding_function_trial.csv",header = T)

colnames(RFU) <- c("treatment","replicate","group","length","fluor1","fluor2","fluor3",
                   "fluor4","fluor5","fluor6","fluor7","fluor8","fluor.av","time")

#convert ml food/Liter to cells algae per ml
RFU$treatment <- ((RFU$treatment / 1000) * 1000000)

#OK, clearance rates in ml calculated as 
#(ln(initial concentration) - ln(final concentration))* Volume / time

##First we should get linear relationship between fluorescence and concentration

control.full <- RFU[which(RFU$group == "full"),]

##SWITCH TERMS
full.reg <- lm(data = control.full, treatment ~ fluor.av)

summary(full.reg)

plot(y = control.full$treatment, x = control.full$fluor.av)

##Now calculate clearance rate for each individual

feeders <- RFU[which(RFU$group == "feed"),]

for(i in 1:length(feeders$treatment)){
  feeders$initial.rfu[i] <- control.full$fluor.av[floor((i-1)/10)+1]
}

feeders$initial.con <- full.reg$coefficients[1] + full.reg$coefficients[2] * feeders$initial.rfu
feeders$final.con <- full.reg$coefficients[1] + full.reg$coefficients[2] * feeders$fluor.av

feeders$volume <- 10


feeders$clearance <- (log(feeders$initial.con) - log(feeders$final.con))* feeders$volume / feeders$time

feeders$clearance_rel <- feeders$clearance / (feeders$length ^ 2)
##OK, now that we have clearance, we can calculate clearance based on length
##and the parameters we want to estimate

feeders <- na.omit(feeders)

#feeding.constant * (L.u[i]^2) * ((cells/volume)/(saturation.constant + (cells/volume)))

#OK, let's first test to see whether there's a humped shape relationship

#test distribution. Because has negatives, will be normal

type_decision <- lm(data = feeders, clearance_rel ~ initial.con + I(initial.con^2))

summary(type_decision)

plot(x=feeders$initial.con,y=feeders$clearance_rel)

#filter out negative clearance values
feeders <- feeders[which(feeders$clearance > 0),]

#make a column for algae cells filtered
feeders$intake <- feeders$clearance * feeders$initial.con

plot(x = feeders$initial.con, y = feeders$intake)
#ok that does look saturating
##OK the function is u shaped, so the exact opposite of type three

#Ok, need to make sure it's food eaten and not clearance
model.loc = ("simple_model.txt")

jagsscript = cat("
                 model {  
                  
                 #Parameters to estimate
                 feeding.constant ~ dunif(0,100000)
                 saturation.constant ~ dunif(0,100000)

                 intake.var ~ dunif(0,100000000) 
                 intake.tau <- 1/(intake.var ^ 2) 
                 
                  
                  for(i in 1:67){
                  
                  intake.estimate[i] <- (feeding.constant * (length[i]^2) * 
                  ((concentration[i])/(saturation.constant + (concentration[i]))))
                  
                  intake.real[i] ~ dnorm(intake.estimate[i],intake.tau)

                  }

                 }
                 
                 ", 
                 file = model.loc)


jags.data = list(length = feeders$length,
                 intake.real = feeders$intake,
                 concentration = feeders$initial.con)

jags.params = c("feeding.constant","saturation.constant")


type_two_model = jags(jags.data, parameters.to.save = jags.params, 
                      model.file = model.loc, n.chains = 3, n.burnin = 1000, n.thin = 5, 
                      n.iter = 10000, DIC = TRUE)#, inits = inits)

type_two_model$BUGSoutput$summary

mcmcplot(type_two_model)

save(type_two_model, file = "type2.Rdata")


##OK, so the half saturation constnt is 10695 cells/ml
##AKA aproximately 10 ml/Liter
##Option 1: 40, 12, 4 ml/L 
##Option 2: 30, 10, 3.5 ml/L. Let's go with that one

# 30/1000 = x/30, 30*(30/1000), x = 0.9 ml
#x = 0.3
# x = 0.1


