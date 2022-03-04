library(dplyr)
library(tidyr)
library(lubridate)

source("./functions_auxiliary.R")

# Load the seroassay used in each location with multiple timepoints
ifrSeroassays <- read.csv("../data/raw_data/location_seroassays.csv",
                          stringsAsFactors=FALSE)

############
# Data of known time between diagnosis and serosurvey
############
seroKnown <- read.csv("../data/raw_data/seroreversion_prior_positives_known_retest.csv",
                     stringsAsFactors=FALSE) %>%
  dplyr::mutate(., test.name=str_replace(test.name, " $", ""),
                number.of.seropositives.among.prior.positives=
                  as.integer(number.of.seropositives.among.prior.positives),
                number.of.prior.positives=as.integer(number.of.prior.positives),
                usedIFR=test.name %in% ifrSeroassays$Assay) %>%
  as_tibble(.)

# Change variable names to be more friendly
newVarNames <- c("phase_id", "country", "location", "sampleType", "startDate",
                 "endDate", "midpointDate", "testTime", "testName",
                 "citationID", "nSeropositives", "nSamples",
                 "sensitivityMean", "sensitivityL", "sensitivityH",
                 "includedTable", "notes", "usedIFR")
names(seroKnown) <- newVarNames

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

# Make the plots
seroprevKnownPlot_IFR <- dplyr::filter(seroKnown, usedIFR==1) %>%
  ggplot(aes(x=testTime, y=sensitivityMean, color=citationID, shape=sampleType)) +
  geom_pointrange(aes(ymin=sensitivityL, ymax=sensitivityH),
                  position=position_jitter(width=0.15, height=0)) +
  facet_wrap(.~testName) +
  theme_bw() +
  theme(legend.position="top",
        ) +
  guides(color=FALSE, shape=guide_legend("Sample type")) +
  xlim(0, 10) +
  xlab("Diagnosis to test (months)") +
  ylab("Sensitivity (%)")

seroprevKnownPlot_noIFR <- dplyr::filter(seroKnown, usedIFR==0) %>%
  ggplot(aes(x=testTime, y=sensitivityMean, color=citationID,
             shape=sampleType)) +
  geom_pointrange(aes(ymin=sensitivityL, ymax=sensitivityH),
                  position=position_jitter(width=0.15, height=0)) +
  facet_wrap(.~testName, ncol=4) +
  theme_bw() +
  theme(legend.position="top") +
  guides(color=FALSE, shape=guide_legend("Sample type")) +
  xlim(0, 12) +
  xlab("Diagnosis to test (months)") +
  ylab("Sensitivity (%)")

seroprevKnownPlot <- dplyr::filter(seroKnown) %>%
  ggplot(aes(x=testTime, y=sensitivityMean, color=testName)) +
  geom_pointrange(aes(ymin=sensitivityL, ymax=sensitivityH),
                  position=position_jitter(width=0.15, height=0)) +
  facet_wrap(.~sampleType) +
  theme_bw() +
  theme(legend.position="none") +
  xlim(0, 12) +
  xlab("Diagnosis to test (months)") +
  ylab("Sensitivity (%)")

ggsave("../data/figures/seroreversion_inIFR_with_diagnosis.png",
       seroprevKnownPlot_IFR, units="cm", width=25, height=16)
ggsave("../data/figures/seroreversion_notIFR_with_diagnosis.png",
       seroprevKnownPlot_noIFR, units="cm", width=32, height=38)
ggsave("../data/figures/seroreversion_all_with_diagnosis.png",
       seroprevKnownPlot, units="cm", width=30, height=18)


############
# Data of unknown time between diagnosis and serosurvey
############
seroUnknown <- read.csv("../data/raw_data/seroreversion_prior_positives_unknown_retest.csv",
                     stringsAsFactors=FALSE) %>%
  dplyr::mutate(., test.name.for.later.sampling=
                str_replace(test.name.for.later.sampling, " $", ""),
                number.of.later.seropositives.among.initial.seropositives=
                  as.integer(number.of.later.seropositives.among.initial.seropositives),
                number.of.initial.seropositives=as.integer(number.of.initial.seropositives),
                usedIFR=test.name.for.later.sampling %in% ifrSeroassays$Assay) %>%
  as_tibble(.)

# Change variable names to be more friendly
newVarNames <- c("phase_id", "country", "location", "sampleType",
                 "midpointInitial", "midpointFollowup", "testTime",
                 "testName", "citationID", "nSeropositives", "nSamples",
                 "sensitivityMean", "sensitivityL", "sensitivityH",
                 "includedTable", "notes", "usedIFR")
names(seroUnknown) <- newVarNames

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

# Make the plots
seroprevUnknownPlot_IFR <- dplyr::filter(seroUnknown, usedIFR==1) %>%
  ggplot(aes(x=testTime, y=sensitivityMean, color=citationID, shape=sampleType)) +
  geom_pointrange(aes(ymin=sensitivityL, ymax=sensitivityH),
                  position=position_jitter(width=0.15, height=0)) +
  facet_wrap(.~testName) +
  theme_bw() +
  theme(legend.position="top",
        ) +
  guides(color=FALSE, shape=guide_legend("Sample type")) +
  xlim(0, 12) +
  xlab("Test to retest (months)") +
  ylab("Relative sensitivity (%)")

seroprevUnknownPlot_noIFR <- dplyr::filter(seroUnknown, usedIFR==0) %>%
  ggplot(aes(x=testTime, y=sensitivityMean, color=citationID,
             shape=sampleType)) +
  geom_pointrange(aes(ymin=sensitivityL, ymax=sensitivityH),
                  position=position_jitter(width=0.15, height=0)) +
  facet_wrap(.~testName, ncol=3) +
  theme_bw() +
  theme(legend.position="top") +
  guides(color=FALSE, shape=guide_legend("Sample type")) +
  xlim(0, 12) +
  xlab("Test to retest (months)") +
  ylab("Relative sensitivity (%)")

seroprevUnknownPlot <- dplyr::filter(seroUnknown) %>%
  ggplot(aes(x=testTime, y=sensitivityMean, color=testName)) +
  geom_pointrange(aes(ymin=sensitivityL, ymax=sensitivityH),
                  position=position_jitter(width=0.15, height=0)) +
  facet_wrap(.~sampleType) +
  theme_bw() +
  theme(legend.position="none") +
  xlim(0, 12) +
  xlab("Test to retest (months)") +
  ylab("Relative sensitivity (%)")

ggsave("../data/figures/seroreversion_inIFR_without_diagnosis.png",
       seroprevUnknownPlot_IFR, units="cm", width=25, height=16)
ggsave("../data/figures/seroreversion_notIFR_without_diagnosis.png",
       seroprevUnknownPlot_noIFR, units="cm", width=32, height=38)
ggsave("../data/figures/seroreversion_all_without_diagnosis.png",
       seroprevUnknownPlot, units="cm", width=30, height=18)


