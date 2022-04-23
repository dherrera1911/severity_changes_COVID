library(dplyr)
library(tidyr)
library(lubridate)
library(rstan)
library(bayesplot)
library(tidybayes)

source("./functions_auxiliary.R")

# Load the seroassay used in each location with multiple timepoints
ifrSeroassays <- read.csv("../data/raw_data/location_seroassays.csv",
                          stringsAsFactors=FALSE)

############
############
# Data of known time between diagnosis and serosurvey
############
############
seroKnown <- read.csv("../data/raw_data/PCR_to_serotest_known_times.csv",
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
                 "includedTable", "notes", "usedIFR")
names(seroKnown) <- newVarNames

# Relabel test times that give intervals into single times
seroKnown$testTime[seroKnown$testTime=="0.5 - 2"] <- "1"
seroKnown$testTime[seroKnown$testTime=="1 - 3"] <- "2"
seroKnown$testTime[seroKnown$testTime=="2 - 3"] <- "2.5"
seroKnown$testTime[seroKnown$testTime=="2 - 4"] <- "3"
seroKnown$testTime[seroKnown$testTime=="4 - 5"] <- "4.5"
seroKnown$testTime[seroKnown$testTime=="4 - 6"] <- "5"
seroKnown$testTime[seroKnown$testTime=="5 - 7"] <- "6"
seroKnown$testTime[seroKnown$testTime==">6"] <- "7"
# Convert test times to integer
seroKnown$testTime <- as.numeric(seroKnown$testTime)

############
############
# Data of unknown time between diagnosis and serosurvey
############
############
seroEstimated <- read.csv("../data/analysis_results/PCR_to_serotest_estimated_times.csv",
                          stringsAsFactors=FALSE)
seroEstimated$includedTable <- NA
seroEstimated <- dplyr::filter(seroEstimated, !(is.na(nSamples) & is.na(sensitivityL)))


### Put together known and unknown serotests
seroAll <- rbind(seroKnown, seroEstimated)

##### Select which dataset to fit
seroFitted <- seroAll

# Remove mixed assays
seroFitted <- dplyr::filter(seroFitted, !stringr::str_detect(testName, "OR"))

# Give confidence intervals to datapoints lacking them (only N samples and N positive)
for (r in c(1:nrow(seroFitted))) {
  if (is.na(seroFitted[[r,"sensitivityL"]])) {
    seroprevConfint <- binomial_confint(seroFitted[[r,"nSamples"]],
                                        seroFitted[[r, "nSeropositives"]])
    seroFitted[[r,"sensitivityL"]] <- signif(seroprevConfint$lower*100, digits=4)
    seroFitted[[r,"sensitivityH"]] <- signif(seroprevConfint$upper*100, digits=4)
    seroFitted[[r,"sensitivityMean"]] <- signif(seroFitted[[r,"nSeropositives"]]/
                                              seroFitted[[r,"nSamples"]]*100, digits=4)
  }
}

# Fit beta distribution to studies without raw data, to estimate raw data
nonRaw <- is.na(seroFitted$nSamples)
fittedBetas <- fit_beta_ci(meanEstimate=seroFitted$sensitivityMean[nonRaw]/100,
                            lower=seroFitted$sensitivityL[nonRaw]/100,
                            upper=seroFitted$sensitivityH[nonRaw]/100)

totalCount <- with(fittedBetas, shape1+shape2)
meanSensitivity <- with(fittedBetas, shape1/(shape1+shape2)*100)
confintSensitivity <- binomial_confint(round(totalCount), round(fittedBetas$shape1))

# Add estimated raw data to the dataframe
seroFitted$nSamples[nonRaw] <- with(fittedBetas, round(shape1+shape2))
seroFitted$nSeropositives[nonRaw] <- round(fittedBetas$shape1)

# filter out increasing sensitivity test
seroFitted <- dplyr::filter(seroFitted, !(testName=="Vitros Ortho total Ig anti-spike"
                                        & citationID=="n-132"))

write.csv(seroFitted, "../data/analysis_results/PCR_to_serotest_estimated_times.csv",
          row.names=FALSE)

###################
###################
# Fit model
###################
###################

# compile model
outcome_reg <- rstan::stan_model("./sensitivity_change.stan",
                                 model_name="time_change_sensitivity",
                                 warn_pedantic=TRUE)

# Fitting parameters
nChains <- 4
nCores <- 2
nIter <- 4000

### Estimate some initial values for the fit
# compute all logits
logitVals <- with(seroFitted, log((sensitivityMean/100)/(1-sensitivityMean/100)))
# mean logits at time 1
meanIntercept <- mean(logitVals[seroFitted$testTime<=1 & logitVals!=Inf])
sdIntercept <- sd(logitVals[seroFitted$testTime<=1] & logitVals!=Inf)
# slopes
seroFitted$normalizedTime <- with(seroFitted, testTime/mean(testTime))

logReg <- glm(cbind(nSeropositives, nSamples-nSeropositives) ~ normalizedTime,
    data=seroFitted, family="binomial")
meanSlope <- logReg$coefficients[2]

####################
# define function to make list of initial values, for STAN
####################
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
  "assaySlope", "assayIntercept", "studyIntercept")
lowerUni <- c(-1, meanIntercept*0.9, 0.2, sdIntercept*0.9,
              -0.6, meanIntercept*0.9, -0.1)
upperUni <- c(0, meanIntercept*1, 0.25, sdIntercept*1.1,
              -0.5, meanIntercept*1, 0.1)
nTests <- length(unique(seroFitted$testName))
nStudies <- length(unique(paste(seroFitted$testName, seroFitted$citationID)))
paramSize <- c(1, 1, 1, 1, nTests, nTests, nStudies)

# Sample the initial values to use
initList <- initial_values(nChains=nChains, paramListName=paramListName,
                           lowerUni=lowerUni, upperUni=upperUni,
                           paramSize=paramSize)


# Make a list with the input we pass to STAN
assayVec <- as.factor(seroFitted$testName)
studies <- as.factor(paste(seroFitted$testName, seroFitted$citationID))
assayDataList <- list(N=nrow(seroFitted),
                  K=length(unique(seroFitted$testName)),
                  M=length(unique(studies)),
                  assay=as.integer(assayVec),
                  study=as.integer(studies),
                  timeVec=seroFitted$testTime/4,
                  nPositive=seroFitted$nSeropositives,
                  nTested=seroFitted$nSamples)

# Fit model
model <- rstan::sampling(outcome_reg, data=assayDataList,
                           chains=nChains, iter=nIter, refresh=0,
                           verbose=TRUE, cores=nCores, init=initList)

######
# Extract model data, and check diagnostics
######
posteriorTraces <- tidybayes::gather_draws(model, timeSlope,
                                           intercept, slopeSigma,
                                           interceptSigma,
                                           studySigma,
                                           assayIntercept[loc],
                                           assaySlope[loc],
                                           studyIntercept[stu])

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

write.csv(posteriorTraces, "../data/analysis_results/sensitivity_decay_all.csv",
          row.names=FALSE)

sensitivityTrace <- dplyr::mutate(posteriorTraces, .chain=factor(.chain))  %>%
  ggplot(., aes(x=.iteration, y=.value, color=.chain)) +
  geom_line() +
  facet_wrap(.~.variable, scales="free") +
  theme_bw()

pairsPlot <- pairs(model, pars=c("timeSlope", "intercept", "slopeSigma",
                                 "interceptSigma", "studySigma"))
#png("../data/figures/simulated_data_fit_pairs_plot.png")

