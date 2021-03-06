##################################
##################################
# 
# This script fits a model to data on time-varying
# sensitivity of different serology assays. Only data
# of serology testing on previously diagnosed individuals
# is used.
# 
# For the studies without known time from diagnosis to
# sero-testing, the diagnosis time was estimated (in
# script 05).
# 
# The script also does some fitting leaving out some of
# the data.
# 
##################################
##################################

library(dplyr)
library(tidyr)
library(lubridate)
library(rstan)
library(bayesplot)
library(tidybayes)

source("./functions_auxiliary.R")

# Fitting parameters
nChains <- 4
nCores <- 2
nIter <- 4000

############
# 1) Load data generated in script 06
############
seroFitted <- read.csv("../data/analysis_results/PCR_to_serotest_all.csv",
                       stringsAsFactors=FALSE)

############
# 2) Load file with assay characteristics
############
assayChars <- read.csv("../data/raw_data/assay_characteristics.csv",
                       stringsAsFactors=FALSE, na.strings=c("", "--")) %>%
  dplyr::rename(., testName=test_name_long, antigen_target=antigen.target,
  assay_type=assay.type, isotype_target=isotype.target)

# add assay characteristics column to the data on sensitivity
antigenVec <- NULL
antibodyVec <- NULL
techniqueVec <- NULL
for (r in c(1:nrow(seroFitted))) {
  charRow <- assayChars$testName == seroFitted$testName[r]
  if (any(charRow )) {
    antigenVec <- c(antigenVec, assayChars$antigen_target[charRow])
    antibodyVec <- c(antibodyVec, assayChars$isotype_target[charRow])
    techniqueVec <- c(techniqueVec, assayChars$assay_type[charRow])
  } else {
    antigenVec <- c(antigenVec, NA)
    antibodyVec <- c(antibodyVec, NA)
    techniqueVec <- c(techniqueVec, NA)
  }
}
seroFitted$antigenTarget <- antigenVec
seroFitted$antibodyTarget <- antibodyVec
seroFitted$technique <- techniqueVec

# Add binary columns indicating presence or absence of characteristics
seroFitted$spikeAntigen <- str_detect(seroFitted$antigenTarget, "S") |
  str_detect(seroFitted$antigenTarget, "RBD")
seroFitted$rbdAntigen <- str_detect(seroFitted$antigenTarget, "RBD")
seroFitted$nuceloAntigen <- str_detect(seroFitted$antigenTarget, "N")
seroFitted$LFA <- str_detect(seroFitted$technique, "LFIA")
seroFitted$ELISA <- str_detect(seroFitted$technique, "ELISA")
seroFitted$chemlum <- (str_detect(seroFitted$technique, "CMIA") |
                       str_detect(seroFitted$technique, "CLIA"))
seroFitted$IgM <- str_detect(seroFitted$antibodyTarget, "IgM")
seroFitted$IgA <- str_detect(seroFitted$antibodyTarget, "IgA")

### Filter rows without value for the characteristic of interest
seroFitted <- dplyr::filter(seroFitted, !is.na(seroFitted$antigenTarget))

#nChars <- 4
#charMatrix <- matrix(as.integer(c(seroFitted$spikeAntigen, seroFitted$nuceloAntigen,
#                                  seroFitted$rbdAntigen,
#                                  seroFitted$spikeAntigen*seroFitted$nuceloAntigen)),
#                     ncol=nChars)
#charsName <- c("S", "N", "RBD", "SN")
#fileIdentifier <- "antigenTarget"

#nChars <- 2
##charMatrix <- matrix(as.integer(c(seroFitted$LFA, seroFitted$ELISA,
##                                  seroFitted$chemlum)), ncol=nChars)
#charMatrix <- matrix(as.integer(c(seroFitted$LFA, !seroFitted$LFA)), ncol=nChars)
#charsName <- c("LFA", "Rest")
#fileIdentifier <- "technique"

nChars <- 4
charMatrix <- matrix(as.integer(c(rep(1, nrow(seroFitted)),
                                  seroFitted$IgM, seroFitted$IgA,
                                  seroFitted$IgM & seroFitted$IgA)), ncol=nChars)
charsName <- c("IgG", "IgG+IgM", "IgG+IgA", "Total")
fileIdentifier <- "antibody"


###################
# 3) Get initial values for the fitting
###################

### Estimate some initial values for the fit
# compute all logits
logitVals <- with(seroFitted, log((sensitivityMean/100)/(1-sensitivityMean/100)))
# mean logits at time 1
meanIntercept <- mean(logitVals[seroFitted$testTime<=1 & logitVals!=Inf])
sdIntercept <- sd(logitVals[seroFitted$testTime<=1] & logitVals!=Inf)
# slopes
seroFitted$normalizedTime <- with(seroFitted, testTime/mean(testTime))
# fit glm and get slope
logReg <- glm(cbind(nSeropositives, nSamples-nSeropositives) ~ normalizedTime,
    data=seroFitted, family="binomial")
meanSlope <- logReg$coefficients[2]

## Define function to make list of initial values, for STAN
initial_values <- function(nChains, paramListName, lowerUni, upperUni, paramSize) {
  initList <- list()
  for (ch in c(1:nChains)) {
    initList[[ch]] <- list()
    for (p in c(1:length(paramListName))) {
      initList[[ch]][[paramListName[p]]] <- runif(paramSize[p], min=lowerUni[p],
                                                  max=upperUni[p])
    }
  }
  return(initList)
}
# Set the ranges for the initial values of the parameters
paramListName <- c("timeSlope", "intercept", "slopeSigma", "interceptSigma",
  "assaySlope", "assayIntercept", "studyIntercept", "charSlope")
lowerUni <- c(-1, meanIntercept*0.9, 0.2, sdIntercept*0.9,
              -0.6, meanIntercept*0.9, -0.1, -0.1)
upperUni <- c(0, meanIntercept*1, 0.25, sdIntercept*1.1,
              -0.5, meanIntercept*1, 0.1, 0.1)
nTests <- length(unique(seroFitted$testName))
nStudies <- length(unique(paste(seroFitted$testName, seroFitted$citationID)))
paramSize <- c(1, 1, 1, 1, nTests, nTests, nStudies, nChars)

# Sample the initial values to use
initList <- initial_values(nChains=nChains, paramListName=paramListName,
                           lowerUni=lowerUni, upperUni=upperUni,
                           paramSize=paramSize)

###################
# 5) Fit the model to the data
###################

# compile model
outcome_reg <- rstan::stan_model("./sensitivity_change_assay_chars.stan",
                                 model_name="time_change_sensitivity",
                                 warn_pedantic=TRUE)

# Make a list with the input we pass to STAN
assayVec <- as.factor(seroFitted$testName)
studies <- as.factor(paste(seroFitted$testName, seroFitted$citationID))
assayDataList <- list(N=nrow(seroFitted),
                  K=length(unique(seroFitted$testName)),
                  M=length(unique(studies)),
                  C=nChars,
                  characteristics=charMatrix,
                  assay=as.integer(assayVec),
                  study=as.integer(studies),
                  timeVec=seroFitted$testTime/4,
                  nPositive=seroFitted$nSeropositives,
                  nTested=seroFitted$nSamples)

# Fit model
model <- rstan::sampling(outcome_reg, data=assayDataList,
                           chains=nChains, iter=nIter, refresh=0,
                           verbose=TRUE, cores=nCores, init=initList)

###################
# 6) Extract and save the posterior samples of the fit
###################
posteriorTraces <- tidybayes::gather_draws(model,
                                           intercept, slopeSigma,
                                           interceptSigma,
                                           studySigma,
                                           assayIntercept[loc],
                                           assaySlope[loc],
                                           studyIntercept[stu],
                                           charSlope[par])

posteriorTraces$parName <- charsName[posteriorTraces$par]

#charParams <- dplyr::filter(posteriorTraces, .variable=="charSlope") %>%
#  pivot_wider(., id_cols=.draw, names_from=parName, values_from=.value)
#mean(charParams$RBD > charParams$S)

posteriorTraces$assay <- levels(assayVec)[posteriorTraces$loc]
studyList <- strsplit(as.character(levels(studies)), " n-")
assayVec <- NULL
studyVec <- NULL
for (i in c(1:length(studyList))) {
  assayVec <- c(assayVec, studyList[[i]][1])
  studyVec <- c(studyVec, studyList[[i]][2])
}

posteriorTraces$studyAssay <- assayVec[posteriorTraces$stu]
posteriorTraces$studyCitation <- studyVec[posteriorTraces$stu]

write.csv(posteriorTraces, paste("../data/analysis_results/7_assay_characteristics_",
                              fileIdentifier, ".csv", sep=""), row.names=FALSE)


####################
## 7) Fit the model to subsets of the data, to test performance on left out data
####################
#
## Tests with multiple time points, which we will use for testing
#multiSamples <- group_by(seroFitted, testName) %>%
#  summarize(., nTimes=length(unique(testTime)) > 1) %>%
#  dplyr::filter(., nTimes)
#
## Set the testing groups
#nLeaveOut <- 4
#nGroups <- ceiling(nrow(multiSamples)/nLeaveOut)
#valGroup <- rep(c(1:nGroups), nLeaveOut)[1:nrow(multiSamples)]
#valGroup <- sample(valGroup)
#
#validationDf <- NULL
#
#for (ng in c(1:nGroups)) {
#  leaveOutTests <- multiSamples$testName[valGroup==ng]
#  meanTimes <- dplyr::filter(seroFitted, testName %in% leaveOutTests) %>%
#    group_by(., testName) %>%
#    summarize(., medianTime= mean(testTime))
#
#  # Test on the first half of the assay data, and the second half separately
#  for (testHalf in c(0, 1)) {
#    # Create fitting and testing data frames
#    valSampleDf <- NULL
#    fittingDf <- seroFitted
#    for (lo in c(1:nrow(meanTimes))) {
#      tName <- meanTimes$testName[lo]
#      timesVec <- fittingDf[fittingDf$testName==tName, "testTime"]
#      if (testHalf == 0) {
#        validationRows <- which((fittingDf$testName == tName &
#                                  fittingDf$testTime <= meanTimes$medianTime[lo]))
#      } else {
#        validationRows <- which((fittingDf$testName == tName &
#                                  fittingDf$testTime > meanTimes$medianTime[lo]))
#      } 
#      valSampleDf <- rbind(valSampleDf, fittingDf[validationRows,])
#      fittingDf <- fittingDf[-validationRows,]
#    }
#    valSampleDf$validationTimeGroup <- testHalf
#
#    # Prepare data and initializations for fitting data
#    nStudies <- length(unique(paste(fittingDf$testName, fittingDf$citationID)))
#    paramSize <- c(1, 1, 1, 1, nTests, nTests, nStudies)
#    # Sample the initial values to use
#    initList <- initial_values(nChains=nChains, paramListName=paramListName,
#                               lowerUni=lowerUni, upperUni=upperUni,
#                               paramSize=paramSize)
#
#    # Make a list with the input we pass to STAN
#    assayVec <- as.factor(fittingDf$testName)
#    studies <- as.factor(paste(fittingDf$testName, fittingDf$citationID))
#    assayDataList <- list(N=nrow(fittingDf),
#                      K=length(unique(fittingDf$testName)),
#                      M=length(unique(studies)),
#                      assay=as.integer(assayVec),
#                      study=as.integer(studies),
#                      timeVec=fittingDf$testTime/4,
#                      nPositive=fittingDf$nSeropositives,
#                      nTested=fittingDf$nSamples)
#
#    # Fit model
#    modelVal <- rstan::sampling(outcome_reg, data=assayDataList,
#                               chains=nChains, iter=nIter, refresh=0,
#                               verbose=TRUE, cores=nCores, init=initList)
#
#    posteriorTraces <- tidybayes::gather_draws(modelVal, timeSlope,
#                                               intercept, slopeSigma,
#                                               interceptSigma,
#                                               studySigma,
#                                               assayIntercept[loc],
#                                               assaySlope[loc],
#                                               studyIntercept[stu])
#
#    # Put assay labels to the samples of the assay (not any specific study)
#    posteriorTraces$assay <- levels(assayVec)[posteriorTraces$loc]
#    studyList <- strsplit(as.character(levels(studies)), " n-")
#    # Put assay and study labels to the samples of specific studies
#    studyAssayVec <- NULL
#    studyVec <- NULL
#    for (i in c(1:length(studyList))) {
#      studyAssayVec <- c(studyAssayVec, studyList[[i]][1])
#      studyVec <- c(studyVec, studyList[[i]][2])
#    }
#    posteriorTraces$studyAssay <- studyAssayVec[posteriorTraces$stu]
#    posteriorTraces$studyCitation <- studyVec[posteriorTraces$stu]
#
#    # Extract each row of the validation Df, and see the prediction
#    # of its sensitivity by the model
#    sensitivityMean <- NULL
#    sensitivityL <- NULL
#    sensitivityH <- NULL
#    for (valRow in c(1:nrow(valSampleDf))) {
#      tn <- valSampleDf$testName[valRow]
#      id <- strsplit(valSampleDf$citationID[valRow], "-")[[1]][2]
#      time <- valSampleDf$testTime[valRow]
#      # Extract posterior samples relevant to this row (assay, study)
#      valRowSamples <- dplyr::filter(posteriorTraces, (studyAssay==tn) |
#                                     (assay==tn) | (.variable=="studySigma"))
#      valRowSensitivity <- NULL # vector to put sensitivity samples
#      if (any(valRowSamples$studyCitation == id, na.rm=TRUE)) {
#        for (dr in unique(valRowSamples$.draw)) {
#          studyDrawDf <- dplyr::filter(valRowSamples, .draw==dr)
#          assaySlope <- with(studyDrawDf, .value[.variable=="assaySlope"])
#          assayIntercept <- with(studyDrawDf, .value[.variable=="assayIntercept"])
#          studyIntercept <- with(studyDrawDf, .value[.variable=="studyIntercept"])
#          lin <- assayIntercept + studyIntercept + time/4 * assaySlope
#          valRowSensitivity[length(valRowSensitivity)+1] <- 1/(1+exp(-lin))
#        }
#      } else {
#        for (dr in unique(valRowSamples$.draw)) {
#          studyDrawDf <- dplyr::filter(valRowSamples, .draw==dr)
#          assaySlope <- with(studyDrawDf, .value[.variable=="assaySlope"])
#          assayIntercept <- with(studyDrawDf, .value[.variable=="assayIntercept"])
#          studySigma <- with(studyDrawDf, .value[.variable=="studySigma"])
#          studyIntercept <- rnorm(1, mean=0, sd=studySigma)
#          lin <- assayIntercept + studyIntercept + time/4 * assaySlope
#          valRowSensitivity[length(valRowSensitivity)+1] <- 1/(1+exp(-lin))
#        }
#      }
#      sensitivityMean[length(sensitivityMean)+1] <- mean(valRowSensitivity)
#      confintSensitivity <- quantile(valRowSensitivity, probs=c(0.025, 0.975))
#      sensitivityL[length(sensitivityL)+1] <- confintSensitivity[1]
#      sensitivityH[length(sensitivityH)+1] <- confintSensitivity[2]
#    }
#    valSampleDf$sensitivityMeanPred<- sensitivityMean
#    valSampleDf$sensitivityLPred<- sensitivityL
#    valSampleDf$sensitivityHPred<- sensitivityH
#    validationDf <- rbind(validationDf, valSampleDf)
#  }
#}
#
#write.csv(validationDf, paste("../data/analysis_results/7_assay_characteristics_",
#                              fileIdentifier, ".csv", sep=""), row.names=FALSE)
#
#
#
##sensitivityTrace <- dplyr::mutate(posteriorTraces, .chain=factor(.chain))  %>%
##  ggplot(., aes(x=.iteration, y=.value, color=.chain)) +
##  geom_line() +
##  facet_wrap(.~.variable, scales="free") +
##  theme_bw()
##
#pairsPlot <- pairs(model, pars=c("timeSlope", "intercept", "slopeSigma",
#                                 "interceptSigma", "studySigma"))
##png("../data/figures/simulated_data_fit_pairs_plot.png")
#
