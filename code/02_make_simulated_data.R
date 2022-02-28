#########################
#  
#  This script generates some simulated seroprevalence
#  and death data of different locations with time
#  changing IFR, in order to test the statistical models
#  
#########################

library(dplyr)
library(tidyr)
source("./functions_auxiliary.R")

slope1 <- 0.13
intercept1 <- -12.9
interceptChange <- -0.7
interceptSigma <- 1.11
interceptChangeSigma <- 0.1

ageBins <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69",
             "70-79", "80+")
ageMeans <- mid_bin_age(ageBins)
nAges <- length(ageBins)

nLocations <- 12

populations <- round(10^rnorm(nLocations, mean=5, sd=0.6))

propInfected1 <- runif(nLocations, 0.05, 0.2) #infected first wave
propInfectedInter <- runif(nLocations, 0.05, 0.3) #infected in interval
propInfected2 <- propInfected1 + propInfectedInter #infected total

seroprevalenceSamples <- sample(50:250, nLocations, replace=T)

allDf <- NULL
for (l in c(1:nLocations)) {
  # Get the IFR for each age and period for this location 
  intercept1Loc <- rnorm(1, intercept1, interceptSigma)
  interceptChangeLoc <- rnorm(1, interceptChange, interceptChangeSigma)
  intercept2Loc <- intercept1Loc + interceptChangeLoc
  lin1 <- intercept1Loc + slope1*ageMeans
  ifr1 <- exp(lin1)/(1 + exp(lin1))
  lin2 <- intercept2Loc + slope1*ageMeans
  ifr2 <- exp(lin2)/(1 + exp(lin2))

  # Sample the number of deaths
  agePops <- rep(round(populations[l]/nAges), nAges)
  realInfections1 <- round(agePops*propInfected1[l])
  realInfectionsInter <- round(agePops*propInfectedInter[l])
  deaths1 <- rbinom(nAges, realInfections1, ifr1)
  deathsInter <- rbinom(nAges, realInfectionsInter, ifr2)
  deaths2 <- deaths1 + deathsInter

  # Simulate the seroprevalence study
  positives1 <- rbinom(nAges, seroprevalenceSamples[l], propInfected1[l])
  seroprev1 <- signif(positives1/seroprevalenceSamples[l]*100, digits=3)
  seroprevCI1 <- binomial_confint(rep(seroprevalenceSamples[l], nAges), positives1)
  seroprevL1 <- signif(seroprevCI1$lower*100, digits=3)
  seroprevH1 <- signif(seroprevCI1$upper*100, digits=3)
  positives2 <- rbinom(nAges, seroprevalenceSamples[l], propInfected2[l])
  seroprev2 <- signif(positives2/seroprevalenceSamples[l]*100, digits=3)
  seroprevCI2 <- binomial_confint(rep(seroprevalenceSamples[l], nAges), positives2)
  seroprevL2 <- signif(seroprevCI2$lower*100, digits=3)
  seroprevH2 <- signif(seroprevCI2$upper*100, digits=3)

  # Put data into dataframe
  locationDf <- data.frame(Age=ageBins, Population=agePops,
                           Prevalence1=seroprev1, Prevalence1L=seroprevL1,
                           Prevalence1H=seroprevH1, Deaths1=deaths1,
                           Prevalence2=seroprev2, Prevalence2L=seroprevL2,
                           Prevalence2H=seroprevH2, Deaths2=deaths2,
                           Location=l)
  allDf <- rbind(allDf, locationDf)
}
  
write.csv(allDf, file="../data/raw_data/simulated_data.csv", row.names=FALSE)

