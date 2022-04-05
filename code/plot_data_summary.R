library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)


#########################
#########################
# Summary of serosurveys to be used for IFR estimation
#########################
#########################

allData <- read.csv("../data/raw_data/collected_data.csv", stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., MidPointCases=lubridate::date(MidPointCases),
                 EndPointOutcome=lubridate::date(EndPointOutcome),
                 Location=factor(Location))

### Plot the dates of the seroprevalence studies
countriesDates <- dplyr::select(allData, Location, MidPointCases, Period) %>%
  dplyr::distinct(.)
factorLevels <- as.character(unique(countriesDates$Location[order(countriesDates$MidPointCases)]))
countriesDates$Location <- factor(countriesDates$Location, levels=factorLevels)

deathNumbers <- dplyr::group_by(allData, Location, Period) %>%
  dplyr::filter(., !is.na(Deaths)) %>%
  dplyr::summarize(., totalDeaths=sum(Deaths))
firstDeaths <- dplyr::filter(deathNumbers, Period==1)
relativeDeaths <- deathNumbers
relativeDeaths$relativeDeaths <- NA
relativeDeaths$absoluteDeathChange <- NA

for (r in c(1:nrow(relativeDeaths))) {
  deaths1 <- firstDeaths[firstDeaths$Location==relativeDeaths$Location[r],]
  relativeDeaths$relativeDeaths[r] <- relativeDeaths$totalDeaths[r]/deaths1$totalDeaths
  relativeDeaths$absoluteDeathChange[r] <- relativeDeaths$totalDeaths[r]-deaths1$totalDeaths
}

countriesDates <- merge(countriesDates, relativeDeaths)


datesPlot <- ggplot(countriesDates, aes(x=MidPointCases, y=Location, color=Location)) +
  geom_point(aes(size=relativeDeaths)) +
  geom_line() +
  theme_bw() +
  geom_vline(xintercept=lubridate::date("2021-01-01")) +
  geom_vline(xintercept=lubridate::date("2020-08-01")) +
  scale_color_discrete(guide="none") +
  scale_size_continuous(name="Relative deaths") +
  xlab("Serosurvey median date")

ggsave("../data/figures/seroprevalence_dates.png", datesPlot,
       width=22, height=12, units="cm")

### Plot deaths in time for each location
library(COVID19)

countries <- c("Spain", "Denmark")

countriesDeaths <- covid19(country=countries, level=1) %>%
  dplyr::mutate(., Location=administrative_area_level_1) %>%
  dplyr::filter(., !is.na(deaths)) %>%
  ungroup(.) %>%
  dplyr::select(., date, deaths, Location)

franceDeaths <- read.csv("../data/downloaded_datasets/france/owid-covid-data.csv",
                         stringsAsFactors=FALSE) %>%
  dplyr::filter(., location=="France" & !is.na(total_deaths)) %>%
  dplyr::mutate(., deaths=total_deaths, Location="France") %>%
  dplyr::select(., date, deaths, Location) 

georgiaDeaths <- covid19(country="USA", level=2) %>%
  dplyr::filter(., administrative_area_level_2=="Georgia") %>%
  dplyr::mutate(., Location="Atlanta") %>%
  dplyr::filter(., !is.na(deaths)) %>%
  ungroup(.) %>%
  dplyr::select(., date, deaths, Location)

ukDeaths <- covid19(country=c("United Kingdom"), level=2)

englandDeaths <- ukDeaths %>%
  dplyr::filter(., administrative_area_level_2=="England") %>%
  dplyr::group_by(., date) %>%
  dplyr::summarize(., deaths=sum(deaths, na.rm=TRUE)) %>%
  dplyr::mutate(., Location="England") %>%
  ungroup(.) %>%
  dplyr::select(., date, deaths, Location)

scotlandDeaths <- ukDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Scotland") %>%
  dplyr::group_by(., date) %>%
  dplyr::summarize(., deaths=sum(deaths, na.rm=TRUE)) %>%
  dplyr::mutate(., Location="Scotland") %>%
  ungroup(.) %>%
  dplyr::select(., date, deaths, Location)

walesDeaths <- ukDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Wales") %>%
  dplyr::group_by(., date) %>%
  dplyr::summarize(., deaths=sum(deaths, na.rm=TRUE)) %>%
  dplyr::mutate(., Location="Wales") %>%
  ungroup(.) %>%
  dplyr::select(., date, deaths, Location)

irelandDeaths <- ukDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Northern Ireland") %>%
  dplyr::group_by(., date) %>%
  dplyr::summarize(., deaths=sum(deaths, na.rm=TRUE)) %>%
  dplyr::mutate(., Location="Ireland") %>%
  ungroup(.) %>%
  dplyr::select(., date, deaths, Location)

# Deaths in Switzerland
genevaDeaths <- read.csv("../data/downloaded_datasets/switzerland/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="GE" & altersklasse_covid19!="Unbekannt") %>%
  #dplyr::group_by(., datum) %>%
  dplyr::group_by(., datum_dboardformated) %>%
  dplyr::summarize(., deaths=sum(sumTotal)) %>%
  dplyr::mutate(., Location="Geneva")
genevaDeaths$date <- paste(sub("-", "-W", genevaDeaths$datum_dboardformated),
                          "7", sep="-")
genevaDeaths$date <- ISOweek::ISOweek2date(genevaDeaths$date)
genevaDeaths <- dplyr::select(genevaDeaths, date, deaths, Location)


neuchatelDeaths <- read.csv("../data/downloaded_datasets/switzerland/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="NE" & altersklasse_covid19!="Unbekannt") %>%
  #dplyr::group_by(., datum) %>%
  dplyr::group_by(., datum_dboardformated) %>%
  dplyr::summarize(., deaths=sum(sumTotal)) %>%
  dplyr::mutate(., Location="Neuchatel")
neuchatelDeaths$date <- paste(sub("-", "-W", neuchatelDeaths$datum_dboardformated),
                          "7", sep="-")
neuchatelDeaths$date <- ISOweek::ISOweek2date(neuchatelDeaths$date)
neuchatelDeaths <- dplyr::select(neuchatelDeaths, date, deaths, Location)

fribourgDeaths <- read.csv("../data/downloaded_datasets/switzerland/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="FR" & altersklasse_covid19!="Unbekannt") %>%
  #dplyr::group_by(., datum) %>%
  dplyr::group_by(., datum_dboardformated) %>%
  dplyr::summarize(., deaths=sum(sumTotal)) %>%
  dplyr::mutate(., Location="Fribourg")
fribourgDeaths$date <- paste(sub("-", "-W", fribourgDeaths$datum_dboardformated),
                          "7", sep="-")
fribourgDeaths$date <- ISOweek::ISOweek2date(fribourgDeaths$date)
fribourgDeaths <- dplyr::select(fribourgDeaths, date, deaths, Location)

ticinoDeaths <- read.csv("../data/downloaded_datasets/switzerland/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="TI" & altersklasse_covid19!="Unbekannt") %>%
  #dplyr::group_by(., datum) %>%
  dplyr::group_by(., datum_dboardformated) %>%
  dplyr::summarize(., deaths=sum(sumTotal)) %>%
  dplyr::mutate(., Location="Ticino")
ticinoDeaths$date <- paste(sub("-", "-W", ticinoDeaths$datum_dboardformated),
                          "7", sep="-")
ticinoDeaths$date <- ISOweek::ISOweek2date(ticinoDeaths$date)
ticinoDeaths <- dplyr::select(ticinoDeaths, date, deaths, Location)

# Germany deaths
dataGer<- read.csv("../data/downloaded_datasets/germany/RKI_COVID19.csv",
                    stringsAsFactors=FALSE) %>%
  tidyr::as_tibble(.)

freiburgDeaths <- dataGer %>%
  dplyr::filter(., Bundesland=="Baden-Württemberg" &
                Landkreis %in% c("SK Freiburg i.Breisgau", "LK Breisgau-Hochschwarzwald")) %>%
                #Landkreis %in% c("SK Freiburg i.Breisgau", "LK Emmendingen", "LK Breisgau-Hochschwarzwald")) %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum), Location="Freiburg") %>%
  dplyr::group_by(., date, Location) %>%
  dplyr::summarize(., deaths=sum(AnzahlTodesfall)) %>%
  dplyr::select(., date, deaths, Location) %>%
  ungroup(.)
freiburgDeaths$deaths <- cumsum(freiburgDeaths$deaths)

reutlingenDeaths <- dataGer %>%
  dplyr::filter(., Bundesland=="Baden-Württemberg" &
                Landkreis=="LK Reutlingen") %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum), Location="Reutlingen") %>%
  dplyr::group_by(., date, Location) %>%
  dplyr::summarize(., deaths=sum(AnzahlTodesfall)) %>%
  dplyr::select(., date, deaths, Location) %>%
  ungroup(.)
reutlingenDeaths$deaths <- cumsum(reutlingenDeaths$deaths)

osnabruckDeaths <- dataGer %>%
  dplyr::filter(., Bundesland=="Niedersachsen" &
                Landkreis %in% c("LK Osnabrück", "SK Osnabrück")) %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum), Location="Osnabruck") %>%
  dplyr::group_by(., date, Location) %>%
  dplyr::summarize(., deaths=sum(AnzahlTodesfall)) %>%
  ungroup(.)
osnabruckDeaths$deaths <- cumsum(osnabruckDeaths$deaths)

aachenDeaths <- dataGer %>%
  dplyr::filter(., Bundesland=="Nordrhein-Westfalen" &
                Landkreis=="StädteRegion Aachen") %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum), Location="Aachen") %>%
  dplyr::group_by(., date, Location) %>%
  dplyr::summarize(., deaths=sum(AnzahlTodesfall)) %>%
  ungroup(.)
aachenDeaths$deaths <- cumsum(aachenDeaths$deaths)

magdeburgDeaths <- dataGer %>%
  dplyr::filter(., Bundesland=="Sachsen-Anhalt" &
                Landkreis=="SK Magdeburg") %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum), Location="Magdeburg") %>%
  dplyr::group_by(., date, Location) %>%
  dplyr::summarize(., deaths=sum(AnzahlTodesfall)) %>%
  ungroup(.)
magdeburgDeaths$deaths <- cumsum(magdeburgDeaths$deaths)

munichDeaths <- dataGer %>%
  dplyr::filter(., Bundesland=="Bayern" &
                Landkreis=="SK München") %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum), Location="Munich") %>%
  dplyr::group_by(., date, Location) %>%
  dplyr::summarize(., deaths=sum(AnzahlTodesfall)) %>%
  ungroup(.)
munichDeaths$deaths <- cumsum(munichDeaths$deaths)

allDeaths <- rbind(aachenDeaths, georgiaDeaths, countriesDeaths,
                   englandDeaths, franceDeaths, freiburgDeaths, fribourgDeaths,
                   genevaDeaths, irelandDeaths, magdeburgDeaths, munichDeaths,
                   neuchatelDeaths, osnabruckDeaths, reutlingenDeaths,
                   scotlandDeaths, ticinoDeaths, walesDeaths)

locationList <- unique(allDeaths$Location)
deathPlots <- list()
for (loc in locationList) {
  dates <- dplyr::filter(countriesDates, Location==loc)
  deaths <- dplyr::filter(allDeaths, Location==loc)
  deathPlots[[loc]] <- ggplot(deaths, aes(x=date, y=deaths)) +
    geom_line() +
    theme_bw() +
    ggtitle(loc) +
    ylab("Deaths") +
    xlab("Date")
  for (r in c(1:nrow(dates))) {
    if (dates$MidPointCases[r] < "2021-01-01") {
      lineColor <- "blue"
    } else {
      lineColor <- "red"
    }
    deathPlots[[loc]] <- deathPlots[[loc]] +
      geom_vline(xintercept=dates$MidPointCases[r], color=lineColor)
  }
}

deathPlotsGrid <- ggpubr::ggarrange(plotlist=deathPlots, ncol=4, nrow=5)

ggsave("../data/figures/death_seroprevalence.png", deathPlotsGrid,
       width=30, height=30, units="cm")

#########################
#########################
# Plot the data of assay sensitivity over time
#########################
#########################

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


