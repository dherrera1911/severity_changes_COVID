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
posteriorTraces <- read.csv("../data/analysis_results/sensitivity_decay_all.csv",
                            stringsAsFactors=FALSE)


# Plot model parameters vs their priors
sigmas <- c("interceptSigma", "slopeSigma", "studySigma")
modelSigmas <- dplyr::filter(posteriorTraces, .variable==sigmas)

gammaPars <- list(c(4,4), c(4,8), c(4,4))
priorDf <- NULL
for (ss in c(1:length(sigmas))) {
  x <- seq(0, 2.5, 0.01)
  y <- dgamma(x, gammaPars[[ss]][1], gammaPars[[ss]][2])
  priorDf <- rbind(priorDf, data.frame(x=x, y=y, .variable=sigmas[ss]))
}

sigmasPlot <- ggplot(data=modelSigmas, aes(x=.value)) +
  geom_density() +
  facet_wrap(.~.variable) +
  geom_line(data=priorDf, aes(x=x, y=y), color="red")


# raw data
#seroKnown <- read.csv("../data/raw_data/PCR_to_serotest_known_times.csv",
#                     stringsAsFactors=FALSE) %>%
#  dplyr::mutate(., test.name=str_replace(test.name, " $", ""),
#                number.of.seropositives.among.prior.positives=
#                  as.integer(number.of.seropositives.among.prior.positives),
#                number.of.prior.positives=as.integer(number.of.prior.positives)) %>%
#  as_tibble(.)
## Change variable names to be more friendly
#newVarNames <- c("phase_id", "country", "location", "sampleType", "startDate",
#                 "endDate", "midpointDate", "testTime", "testName",
#                 "citationID", "nSeropositives", "nSamples",
#                 "sensitivityMean", "sensitivityL", "sensitivityH",
#                 "includedTable", "notes")
#names(seroKnown) <- newVarNames
#
## Relabel test times that give intervals into single times
#seroKnown$testTime[seroKnown$testTime=="0.5 - 2"] <- "1"
#seroKnown$testTime[seroKnown$testTime=="1 - 3"] <- "2"
#seroKnown$testTime[seroKnown$testTime=="2 - 3"] <- "2.5"
#seroKnown$testTime[seroKnown$testTime=="2 - 4"] <- "3"
#seroKnown$testTime[seroKnown$testTime=="4 - 5"] <- "4.5"
#seroKnown$testTime[seroKnown$testTime=="4 - 6"] <- "5"
#seroKnown$testTime[seroKnown$testTime=="5 - 7"] <- "6"
#seroKnown$testTime[seroKnown$testTime==">6"] <- "7"
## Convert test times to integer
#seroKnown$testTime <- as.numeric(seroKnown$testTime)
#
## Remove mixed assays
#seroKnown <- dplyr::filter(seroKnown, !stringr::str_detect(testName, "OR"))
#
## Give confidence intervals to datapoints lacking them (only N samples and N positive)
#for (r in c(1:nrow(seroKnown))) {
#  if (is.na(seroKnown[[r,"sensitivityL"]])) {
#    seroprevConfint <- binomial_confint(seroKnown[[r,"nSamples"]],
#                                        seroKnown[[r, "nSeropositives"]])
#    seroKnown[[r,"sensitivityL"]] <- signif(seroprevConfint$lower*100, digits=4)
#    seroKnown[[r,"sensitivityH"]] <- signif(seroprevConfint$upper*100, digits=4)
#    seroKnown[[r,"sensitivityMean"]] <- signif(seroKnown[[r,"nSeropositives"]]/
#                                              seroKnown[[r,"nSamples"]]*100, digits=4)
#  }
#}
#
## Fit beta distribution to studies without raw data, to estimate raw data
#nonRaw <- is.na(seroKnown$nSamples)
#fittedBetas <- fit_beta_ci(meanEstimate=seroKnown$sensitivityMean[nonRaw]/100,
#                            lower=seroKnown$sensitivityL[nonRaw]/100,
#                            upper=seroKnown$sensitivityH[nonRaw]/100)
#
#totalCount <- with(fittedBetas, shape1+shape2)
#meanSensitivity <- with(fittedBetas, shape1/(shape1+shape2)*100)
#confintSensitivity <- binomial_confint(round(totalCount), round(fittedBetas$shape1))
#
## Add estimated raw data to the dataframe
#seroKnown$nSamples[nonRaw] <- with(fittedBetas, round(shape1+shape2))
#seroKnown$nSeropositives[nonRaw] <- round(fittedBetas$shape1)
#
## filter out increasing sensitivity test
#seroKnown <- dplyr::filter(seroKnown, !(testName=="Vitros Ortho total Ig anti-spike"
#                                        & citationID=="n-132"))

seroAll <- read.csv("../data/analysis_results/PCR_to_serotest_estimated_times.csv",
                    stringsAsFactors=FALSE)

seroFitted <- seroAll

testNames <- unique(seroFitted$testName)
timeVec <- seq(0.5, 10, 0.2)
assayFitDf <- NULL
assaySummaryDf <- NULL
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
  summaryStats <- group_by(testTrace, .variable) %>%
    summarize(., meanVal=mean(.value), sdVal=sd(.value))
  assayIntercept <- summaryStats$meanVal[summaryStats$.variable=="assayIntercept"]
  assaySlope <- summaryStats$meanVal[summaryStats$.variable=="assaySlope"]
  assayInterceptSD <- summaryStats$sdVal[summaryStats$.variable=="assayIntercept"]
  assaySlopeSD <- summaryStats$sdVal[summaryStats$.variable=="assaySlope"]
  summaryDf <- data.frame(intercept=assayIntercept,
                          slope=assaySlope,
                          interceptSD=assayInterceptSD,
                          slopeSD=assaySlopeSD, testName=tN)
  assaySummaryDf <- rbind(assaySummaryDf, summaryDf)
}
  
tests1 <- unique(seroFitted$testName)[1:15]
tests2 <- unique(seroFitted$testName)[16:33]

seroprevKnownPlot1 <- dplyr::filter(seroKnown, testName %in% tests1) %>%
  ggplot(aes(x=testTime, y=sensitivityMean, color=citationID,
             shape=sampleType)) +
  geom_pointrange(aes(ymin=sensitivityL, ymax=sensitivityH),
                  position=position_jitter(width=0.15, height=0)) +
  geom_line(data=dplyr::filter(assayFitDf, testName %in% tests1),
            aes(x=time, y=sensitivity*100), color="black", linetype="solid",
            size=regLineSize, inherit.aes=FALSE) +
  geom_ribbon(data=dplyr::filter(assayFitDf, testName %in% tests1),
              aes(x=time, ymin=sensitivityL*100, ymax=sensitivityH*100),
              alpha=ribbonAlpha, colour=NA, show.legend=FALSE,
              inherit.aes=FALSE) +
  facet_wrap(.~testName, ncol=3) +
  theme_bw() +
  theme(legend.position="top") +
  guides(color=FALSE, shape=guide_legend("Sample type")) +
  xlim(0, 12) +
  xlab("Diagnosis to test (months)") +
  ylab("Sensitivity (%)")

seroprevKnownPlot2 <- dplyr::filter(seroKnown, testName %in% tests2) %>%
  ggplot(aes(x=testTime, y=sensitivityMean, color=citationID,
             shape=sampleType)) +
  geom_pointrange(aes(ymin=sensitivityL, ymax=sensitivityH),
                  position=position_jitter(width=0.15, height=0)) +
  geom_line(data=dplyr::filter(assayFitDf, testName %in% tests2),
            aes(x=time, y=sensitivity*100), color="black", linetype="solid",
            size=regLineSize, inherit.aes=FALSE) +
  geom_ribbon(data=dplyr::filter(assayFitDf, testName %in% tests2),
              aes(x=time, ymin=sensitivityL*100, ymax=sensitivityH*100),
              alpha=ribbonAlpha, colour=NA, show.legend=FALSE,
              inherit.aes=FALSE) +
  facet_wrap(.~testName, ncol=3) +
  theme_bw() +
  theme(legend.position="top") +
  guides(color=FALSE, shape=guide_legend("Sample type")) +
  xlim(0, 12) +
  xlab("Diagnosis to test (months)") +
  ylab("Sensitivity (%)")

ggsave("../data/figures/seroreversion_fit1.png", seroprevKnownPlot1,
  units="cm", width=24, height=30)
ggsave("../data/figures/seroreversion_fit2.png", seroprevKnownPlot2,
  units="cm", width=24, height=30)


# Plot parameters
slope_intercept <- ggplot(data=assaySummaryDf, aes(x=intercept, y=slope)) +
  geom_point() +
  theme_bw()
ggsave("../data/figures/seroreversion_assays_correlation.png", slope_intercept,
  units="cm", width=20, height=20)




