library(dplyr)
library(tidyr)
library(lubridate)
library(rstan)
library(bayesplot)
library(tidybayes)

source("./functions_auxiliary.R")

############
# Data of serological test and serological retest
############
sero2sero <- read.csv("../data/raw_data/serotest_to_serotest.csv",
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




