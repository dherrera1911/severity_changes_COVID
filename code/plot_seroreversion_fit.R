library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(rstan)
library(bayesplot)
library(tidybayes)
library(viridis)
library(RColorBrewer)
library(ggthemes)
library(lemon)
library(cowplot)
library(ggpubr)
source("./functions_auxiliary.R")


regLineSize=1
ribbonAlpha=0.2
locLineSize=0.4
locLineAlpha=0.4
locPointSize=1.6
sampleLineSize=0.1
sampleLineAlpha=0.1
dataAlpha <- 0.5


############################
### Load data
############################

# model fit
posteriorTraces <- read.csv("../data/analysis_results/sensitivity_decay.csv",
                            stringsAsFactors=FALSE)

# raw data
seroKnown <- read.csv("../data/raw_data/seroreversion_prior_positives_known_retest.csv",
                     stringsAsFactors=FALSE) %>%
  dplyr::mutate(., test.name=str_replace(test.name, " $", ""),
                number.of.seropositives.among.prior.positives=
                  as.integer(number.of.seropositives.among.prior.positives),
                number.of.prior.positives=as.integer(number.of.prior.positives)) %>%
  as_tibble(.)
# Change variable names to be more friendly
newVarNames <- c("phase_id", "country", "location", "sampleType", "startDate",
                 "endDate", "midpointDate", "testTime", "testName",
                 "citationID", "nSeropositives", "nSamples",
                 "sensitivityMean", "sensitivityL", "sensitivityH",
                 "includedTable", "notes")
names(seroKnown) <- newVarNames

# Relabel test times that give intervals into single times
seroKnown$testTime[seroKnown$testTime=="0.5 - 2"] <- "1"
seroKnown$testTime[seroKnown$testTime=="1 - 3"] <- "2"
seroKnown$testTime[seroKnown$testTime=="2 - 4"] <- "3"
seroKnown$testTime[seroKnown$testTime=="2 - 3"] <- "2.5"
seroKnown$testTime[seroKnown$testTime=="4 - 6"] <- "5"
seroKnown$testTime[seroKnown$testTime=="4 - 6"] <- "5"
seroKnown$testTime[seroKnown$testTime==">6"] <- "7"
# Convert test times to integer
seroKnown$testTime <- as.numeric(seroKnown$testTime)

# Remove mixed assays
seroKnown <- dplyr::filter(seroKnown, !stringr::str_detect(testName, "OR"))

# Give confidence intervals to datapoints lacking them (only N samples and N positive)
for (r in c(1:nrow(seroKnown))) {
  if (is.na(seroKnown[[r,"sensitivityL"]])) {
    seroprevConfint <- binomial_confint(seroKnown[[r,"nSamples"]],
                                        seroKnown[[r, "nSeropositives"]])
    seroKnown[[r,"sensitivityL"]] <- signif(seroprevConfint$lower*100, digits=4)
    seroKnown[[r,"sensitivityH"]] <- signif(seroprevConfint$upper*100, digits=4)
    seroKnown[[r,"sensitivityMean"]] <- signif(seroKnown[[r,"nSeropositives"]]/
                                              seroKnown[[r,"nSamples"]]*100, digits=4)
  }
}

# Fit beta distribution to studies without raw data, to estimate raw data
nonRaw <- is.na(seroKnown$nSamples)
fittedBetas <- fit_beta_ci(meanEstimate=seroKnown$sensitivityMean[nonRaw]/100,
                            lower=seroKnown$sensitivityL[nonRaw]/100,
                            upper=seroKnown$sensitivityH[nonRaw]/100)

totalCount <- with(fittedBetas, shape1+shape2)
meanSensitivity <- with(fittedBetas, shape1/(shape1+shape2)*100)
confintSensitivity <- binomial_confint(round(totalCount), round(fittedBetas$shape1))

# Add estimated raw data to the dataframe
seroKnown$nSamples[nonRaw] <- with(fittedBetas, round(shape1+shape2))
seroKnown$nSeropositives[nonRaw] <- round(fittedBetas$shape1)

# filter out increasing sensitivity test
seroKnown <- dplyr::filter(seroKnown, !(testName=="Vitros Ortho total Ig anti-spike"
                                        & citationID=="n-132"))

testNames <- unique(seroKnown$testName)
timeVec <- seq(0.5, 8, 0.2)
assayFitDf <- NULL
for (tN in testNames) {
  # get posterior 
  testTrace <- dplyr::filter(posteriorTraces, assay==tN)
  sensitivityPosterior <- seroreversion_samples(testTrace, timeVec,
                                                slopeName="assaySlope",
                                                interceptName="assayIntercept",
                                                timeNormalization=4)
  testDf <- data.frame(time=timeVec,
                       sensitivity=sensitivityPosterior$prop_mean,
                       sensitivityL=sensitivityPosterior$prop_L,
                       sensitivityH=sensitivityPosterior$prop_H,
                       testName=tN)
  assayFitDf <- rbind(assayFitDf, testDf)
}
  

seroprevKnownPlot <- dplyr::filter(seroKnown) %>%
  ggplot(aes(x=testTime, y=sensitivityMean, color=citationID,
             shape=sampleType)) +
  geom_pointrange(aes(ymin=sensitivityL, ymax=sensitivityH),
                  position=position_jitter(width=0.15, height=0)) +
  geom_line(data=assayFitDf, aes(x=time, y=sensitivity*100),
            color="black", linetype="solid", size=regLineSize,
            inherit.aes=FALSE) +
  geom_ribbon(data=assayFitDf,
              aes(x=time, ymin=sensitivityL*100, ymax=sensitivityH*100),
              alpha=ribbonAlpha, colour=NA, show.legend=FALSE,
              inherit.aes=FALSE) +
  facet_wrap(.~testName, ncol=4) +
  theme_bw() +
  theme(legend.position="top") +
  guides(color=FALSE, shape=guide_legend("Sample type")) +
  xlim(0, 12) +
  xlab("Diagnosis to test (months)") +
  ylab("Sensitivity (%)")




