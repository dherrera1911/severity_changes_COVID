library(dplyr)
library(tidyr)
library(lubridate)
library(rstan)
library(bayesplot)
library(tidybayes)

source("./functions_auxiliary.R")

############
# Data of unknown time between diagnosis and serosurvey
############
seroUnknown <- read.csv("../data/raw_data/seroreversion_prior_positives_unknown_retest.csv",
                     stringsAsFactors=FALSE) %>%
  dplyr::mutate(., test.name.for.later.sampling=
                str_replace(test.name.for.later.sampling, " $", ""),
                number.of.later.seropositives.among.initial.seropositives=
                  as.integer(number.of.later.seropositives.among.initial.seropositives),
                number.of.initial.seropositives=as.integer(number.of.initial.seropositives)) %>%
  as_tibble(.)

# Change variable names to be more friendly
newVarNames <- c("phase_id", "country", "location", "sampleType",
                 "midpointInitial", "midpointFollowup", "testTime",
                 "testName", "citationID", "nSeropositives",
                 "nSamples", "sensitivityMean", "sensitivityL",
                 "sensitivityH", "includedTable", "notes")

names(seroUnknown) <- newVarNames

# Remove mixed assays
seroUnknown <- dplyr::filter(seroUnknown, !stringr::str_detect(testName, "OR"))

# Give confidence intervals to datapoints lacking them (only N samples and N positive)
seroUnknown$sensitivityL <- as.numeric(seroUnknown$sensitivityL)
seroUnknown$sensitivityH <- as.numeric(seroUnknown$sensitivityH)
for (r in c(1:nrow(seroUnknown))) {
  if (is.na(seroUnknown[[r,"sensitivityL"]])) {
    seroprevConfint <- binomial_confint(seroUnknown[[r,"nSamples"]],
                                        seroUnknown[[r, "nSeropositives"]])
    seroUnknown[[r,"sensitivityL"]] <- signif(seroprevConfint$lower*100, digits=4)
    seroUnknown[[r,"sensitivityH"]] <- signif(seroprevConfint$upper*100, digits=4)
    seroUnknown[[r,"sensitivityMean"]] <- signif(seroUnknown[[r,"nSeropositives"]]/
                                              seroUnknown[[r,"nSamples"]]*100, digits=4)
  }
}


# Relabel test times that give intervals into single times
seroUnknown$testTime[seroUnknown$testTime=="1 - 2"] <- "1.5"
seroUnknown$testTime[seroUnknown$testTime=="3 - 4"] <- "3.5"
seroUnknown$testTime[seroUnknown$testTime=="≥5"] <- "6"
# Convert test times to integer
seroUnknown$testTime <- as.numeric(seroUnknown$testTime)


### Load baseline values
baselineSensitivities <- read.csv("../data/raw_data/seroassay_sensitivity.csv",
                                  stringsAsFactors=FALSE) %>%
  dplyr::mutate(., test_name= str_replace(test_name, " $", "")) %>%
  as_tibble(.)

# Change variable names to be more friendly
newVarNames <- c("Corrected_seroreversion", "test_id", "testName", 
                 "sensitivityMean", "sensitivityL", "sensitivityH",
                 "nTrueSeropositives", "nPositivesPCR",
                 "specificityMean", "specificityL", "specificityH",
                 "nTrueSeronegatives", "nNegativePCR",
                 "citationID", "notes")
names(baselineSensitivities) <- newVarNames

baselineSensitivities <- dplyr::filter(baselineSensitivities,
                                       Corrected_seroreversion=="N/A (baseline)")

baselineNames <- c("Abbott Architect IgG Test", "DiaSorin SARS-CoV-2 S1/S2 IgG",
                   "Wantai SARS-CoV-2 Total Antibody, IgG IgM (against RBD)")
sensitivityNames <- c("Abbott Architect IgG", "Diasorin Liaison S1/S2 spike IgG",
                   "Wantai SARS-CoV-2 Total Antibody, IgG IgM anti-RBD")
for (n in c(1:length(baselineNames))) {
  baselineSensitivities$testName[baselineSensitivities$testName==baselineNames[n]] <-
    sensitivityNames[n]
}





