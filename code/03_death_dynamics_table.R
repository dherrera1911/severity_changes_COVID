library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(COVID19)

source("./functions_auxiliary.R")

###################
## Get the locations that we want death curves for
###################

# Data of IFR analysis
ifrTable <- read.csv("../data/raw_data/collected_data.csv",
                     stringsAsFactors=FALSE)
locNames1 <- unique(ifrTable$Location)

# Data of seroanalysis test and seroanalysis retest
seroUnknown <- read.csv("../data/raw_data/serotest_to_serotest.csv",
                     stringsAsFactors=FALSE) %>%
  as_tibble(.)
# Change variable names to be more friendly
newVarNames <- c("phase_id", "country", "location", "sampleType",
                 "midpointInitial", "midpointFollowup", "testTime",
                 "testName", "citationID", "nSeropositives",
                 "nSamples", "sensitivityMean", "sensitivityL",
                 "sensitivityH", "includedTable", "notes")
names(seroUnknown) <- newVarNames
locNames2 <- sort(unique(paste(seroUnknown$country, seroUnknown$location)))

# Data of PCR test, seroanalysis retest, with unknown dates
seroUnknown2 <- read.csv("../data/raw_data/PCR_to_serotest_unknown_times.csv",
                     stringsAsFactors=FALSE) %>%
  as_tibble(.)
# Change variable names to be more friendly
newVarNames <- c("phase_id", "country", "location", "sampleType",
                 "startDate", "endDate", "midpointDate", 
                 "testName", "citationID", "nSeropositives",
                 "nSamples", "sensitivityMean", "sensitivityL",
                 "sensitivityH", "notes")
names(seroUnknown2) <- newVarNames
locNames3 <- sort(unique(paste(seroUnknown2$country, seroUnknown2$location)))

# Data of PCR test, seroanalysis retest, with known dates
seroKnown <- read.csv("../data/raw_data/PCR_to_serotest_known_times.csv",
                     stringsAsFactors=FALSE) %>%
  as_tibble(.)
# Change variable names to be more friendly
newVarNames <- c("phase_id", "country", "location", "sampleType", "startDate",
                 "endDate", "midpointDate", "testTime", "testName",
                 "citationID", "nSeropositives", "nSamples",
                 "sensitivityMean", "sensitivityL", "sensitivityH",
                 "includedTable", "notes", "usedIFR")
names(seroKnown) <- newVarNames
locNames4 <- sort(unique(paste(seroKnown$country, seroKnown$location)))


########################
# Tabulate the death data for this locations
########################

### Country level
countries <- c("Spain", "Denmark", "Ireland", "Germany", "Netherlands",
               "Portugal", "Canada", "United Kingdom", "Jersey", "India",
               "United States", "Andorra", "Lebanon", "Maldives", "Ethiopia",
               "Lithuania")

covid19_tidying <- function(covid19Df) {
  output <- dplyr::filter(covid19Df, (!is.na(deaths) | !is.na(confirmed))) %>%
    ungroup(.) %>%
    dplyr::select(., date, confirmed, deaths, Location) %>%
    dplyr::mutate(., date=date(date)) %>%
    dplyr::rename(., cases=confirmed)
  return(output)
}

countriesDeaths <- covid19(country=countries) %>%
  dplyr::mutate(., Location=administrative_area_level_1) %>%
  covid19_tidying(.)

countriesDeaths <- read.csv("../data/raw_data/death_dynamics.csv",
                            stringsAsFactors=FALSE) %>%
  dplyr::filter(., Location %in% countries) %>%
  dplyr::mutate(., date=date(date))

franceDeaths <- read.csv("../data/downloaded_datasets/france/owid-covid-data.csv",
                         stringsAsFactors=FALSE) %>%
  dplyr::filter(., location=="France" & (!is.na(total_deaths) | !is.na(total_cases))) %>%
  dplyr::mutate(., deaths=total_deaths, cases=total_cases, Location="France",
                date=date(date)) %>%
  dplyr::select(., date, cases, deaths, Location)

#jersey <- covid19(country="Jersey") %>%
#  dplyr::mutate(., Location=country) %>%
#  ungroup(.) %>%
#  dplyr::select(., date, deaths, Location)

### United States locations
usaDeaths <- read.csv("../data/downloaded_datasets/death_dynamics_covid19_datahub/USA.csv",
                      stringsAsFactors=FALSE)
# states
fullStates <- c("Wisconsin", "California")
usaStatesDeaths <- dplyr::filter(usaDeaths, administrative_area_level_2 %in%
                                 fullStates & administrative_area_level_3=="") %>%
  dplyr::mutate(., Location=administrative_area_level_2) %>%
  covid19_tidying(.)

# counties and cities
locations <- c("Fulton", "Salt Lake", "Hillsborough", "Hampden", "Houston",
  "Los Angeles")
states <- c("Georgia", "Utah", "Florida", "Massachusetts", "Texas",
  "California")
usaCountyDeaths <- dplyr::filter(usaDeaths, (administrative_area_level_3 %in%
  locations) & (administrative_area_level_2 %in% states) &
  !(administrative_area_level_2=="Georgia" & administrative_area_level_3=="Houston")) %>%
  dplyr::mutate(., Location=administrative_area_level_3) %>%
  covid19_tidying(.)
usaCountyDeaths$Location[usaCountyDeaths$Location=="Hampden"] <- "Holyoke"

connecticutDeaths <- dplyr::filter(usaDeaths, (administrative_area_level_2=="Connecticut") &
                                   administrative_area_level_3=="") %>%
  dplyr::mutate(., Location="Connecticut") %>%
  covid19_tidying(.)

### United Kingdom locations
# Nations
englandDeaths <- read.csv("../data/downloaded_datasets/death_dynamics_covid19_datahub/England.csv",
                          stringsAsFactors=FALSE) %>%
  dplyr::mutate(., Location=administrative_area_level_2) %>%
  covid19_tidying(.)

scotlandDeaths <- read.csv("../data/downloaded_datasets/death_dynamics_covid19_datahub/Scotland.csv",
                          stringsAsFactors=FALSE) %>%
  dplyr::mutate(., Location=administrative_area_level_2) %>%
  covid19_tidying(.)

walesDeaths <- read.csv("../data/downloaded_datasets/death_dynamics_covid19_datahub/Wales.csv",
                          stringsAsFactors=FALSE) %>%
  dplyr::mutate(., Location=administrative_area_level_2) %>%
  covid19_tidying(.)

northernIrelandDeaths <- read.csv("../data/downloaded_datasets/death_dynamics_covid19_datahub/Northern_Ireland.csv",
                          stringsAsFactors=FALSE) %>%
  dplyr::filter(., administrative_area_level_2=="Northern Ireland") %>%
  dplyr::mutate(., Location=administrative_area_level_2) %>%
  covid19_tidying(.)

### Switzerland
switzerlandDeaths <- read.csv("../data/downloaded_datasets/switzerland/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE)
switzerlandCases <- read.csv("../data/downloaded_datasets/switzerland/COVID19Cases_geoRegion.csv",
                 stringsAsFactors=FALSE)

swissCovid_tidy <- function(swissDeathDf, swissCasesDf, locCode, locationName) {
  deathsDf <- swissDeathDf %>%
    dplyr::filter(., geoRegion==locCode & altersklasse_covid19!="Unbekannt") %>%
    #dplyr::group_by(., datum) %>%
    dplyr::group_by(., datum_dboardformated) %>%
    dplyr::summarize(., deaths=sum(sumTotal)) %>%
    dplyr::mutate(., Location=locationName)
    deathsDf$date <- paste(sub("-", "-W", deathsDf$datum_dboardformated),
                              "7", sep="-")
    deathsDf$date <- ISOweek::ISOweek2date(deathsDf$date)
    deathsDf <- dplyr::select(deathsDf, date, deaths, Location)
  casesDf <- swissCasesDf %>%
    dplyr::filter(., geoRegion==locCode) %>%
    dplyr::mutate(., cases=sumTotal, Location=locationName) %>%
    dplyr::rename(., date=datum) %>%
    dplyr::select(., date, cases, Location)
  casesDf$date <- date(casesDf$date)
  output <- merge(deathsDf, casesDf)
  return(output)
}

genevaDeaths <- swissCovid_tidy(switzerlandDeaths, switzerlandCases, "GE", "Geneva")
neuchatelDeaths <- swissCovid_tidy(switzerlandDeaths, switzerlandCases, "NE", "Neuchatel")
fribourgDeaths <- swissCovid_tidy(switzerlandDeaths, switzerlandCases, "FR", "Fribourg")
ticinoDeaths <- swissCovid_tidy(switzerlandDeaths, switzerlandCases, "TI", "Ticino")

### Germany deaths
germanCovid_tidy <- function(germanyDf, locationName) {
  output <- germanyDf %>%
    dplyr::mutate(., date=lubridate::date(Meldedatum), Location=locationName) %>%
    dplyr::group_by(., date, Location) %>%
    dplyr::summarize(., deaths=sum(AnzahlTodesfall),
                     cases=sum(AnzahlFall)) %>%
    dplyr::select(., date, cases, deaths, Location) %>%
    ungroup(.)
  output$deaths <- cumsum(output$deaths)
  output$cases <- cumsum(output$cases)
  return(output)
}

dataGer<- read.csv("../data/downloaded_datasets/germany/RKI_COVID19.csv",
                    stringsAsFactors=FALSE) %>%
  tidyr::as_tibble(.)

freiburgDeaths <- dataGer %>%
  dplyr::filter(., Bundesland=="Baden-Württemberg" &
                Landkreis %in% c("SK Freiburg i.Breisgau", "LK Breisgau-Hochschwarzwald")) %>%
  germanCovid_tidy(., "Freiburg")

reutlingenDeaths <- dataGer %>%
  dplyr::filter(., Bundesland=="Baden-Württemberg" &
                Landkreis=="LK Reutlingen") %>%
  germanCovid_tidy(., "Reutlingen")

osnabruckDeaths <- dataGer %>%
  dplyr::filter(., Bundesland=="Niedersachsen" &
                Landkreis %in% c("LK Osnabrück", "SK Osnabrück")) %>%
  germanCovid_tidy(., "Osnabruck")

aachenDeaths <- dataGer %>%
  dplyr::filter(., Bundesland=="Nordrhein-Westfalen" &
                Landkreis=="StädteRegion Aachen") %>%
  germanCovid_tidy(., "Aachen")

magdeburgDeaths <- dataGer %>%
  dplyr::filter(., Bundesland=="Sachsen-Anhalt" &
                Landkreis=="SK Magdeburg") %>%
  germanCovid_tidy(., "Magdeburg")

munichDeaths <- dataGer %>%
  dplyr::filter(., Bundesland=="Bayern" &
                Landkreis=="SK München") %>%
  germanCovid_tidy(., "Munich")

gangeltDeaths <- dataGer %>%
  dplyr::filter(., Landkreis=="LK Heinsberg") %>%
  germanCovid_tidy(., "Gangelt")

mitteDeaths <- dataGer %>%
  dplyr::filter(., Landkreis=="SK Berlin Mitte") %>%
  germanCovid_tidy(., "Mitte")

kupferzellDeaths <- dataGer %>%
  dplyr::filter(., Landkreis=="LK Hohenlohekreis") %>%
  germanCovid_tidy(., "Kupferzell")

tirschenreuthDeaths <- dataGer %>%
  dplyr::filter(., Landkreis=="LK Tirschenreuth") %>%
  germanCovid_tidy(., "Tirschenreuth")

straubingDeaths <- dataGer %>%
  dplyr::filter(., Landkreis=="SK Straubing") %>%
  germanCovid_tidy(., "Straubing")

chemnitzDeaths <- dataGer %>%
  dplyr::filter(., Landkreis=="SK Chemnitz") %>%
  germanCovid_tidy(., "Chemnitz")

vorpommernDeaths <- dataGer %>%
  dplyr::filter(., Landkreis=="LK Vorpommern-Greifswald") %>%
  germanCovid_tidy(., "Vorpommern")

### Austria
dataAut <- read.csv("../data/downloaded_datasets/death_dynamics_covid19_datahub/AUT.csv",
                    stringsAsFactors=FALSE)

ischglDeaths <- dplyr::filter(dataAut, administrative_area_level_3=="Landeck") %>%
  dplyr::mutate(., Location="Ischgl") %>%
  covid19_tidying(.)

wachauDeaths <- dplyr::filter(dataAut, administrative_area_level_3 %in%
                              c("Melk", "Krems an der Donau Land")) %>%
  dplyr::group_by(., date) %>%
  dplyr::summarize(., deaths=sum(deaths), confirmed=sum(confirmed)) %>%
  dplyr::mutate(., Location="Wachau") %>%
  ungroup(.) %>%
  covid19_tidying(.)


### Italy
dataItaly <- read.csv("../data/downloaded_datasets/death_dynamics_covid19_datahub/ITA.csv",
                      stringsAsFactors=FALSE)

trentoDeaths <- dplyr::filter(dataItaly, administrative_area_level_3=="Trento") %>%
  dplyr::mutate(., Location=administrative_area_level_3) %>%
  covid19_tidying(.)
voDeaths <- dplyr::filter(dataItaly, administrative_area_level_3=="Padova") %>%
  dplyr::mutate(., Location="Vo") %>%
  covid19_tidying(.)
addaDeaths <- dplyr::filter(dataItaly, administrative_area_level_3=="Lodi") %>%
  dplyr::mutate(., Location="Adda") %>%
  covid19_tidying(.)
caldariDeaths <- dplyr::filter(dataItaly, administrative_area_level_3=="Chieti") %>%
  dplyr::mutate(., Location="Caldari") %>%
  covid19_tidying(.)

### China
wuhanDeaths <- read.csv("../data/downloaded_datasets/death_dynamics_covid19_datahub/CHN.csv",
                      stringsAsFactors=FALSE) %>%
  dplyr::mutate(., Location="Wuhan") %>%
  covid19_tidying(.)

### Russia
russiaDeaths <- read.csv("../data/downloaded_datasets/death_dynamics_covid19_datahub/RUS.csv",
                         stringsAsFactors=FALSE)

petersburgDeaths <- russiaDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Saint Petersburg") %>%
  dplyr::mutate(., Location="Saint Petersburg") %>%
  covid19_tidying(.)

tyumenDeaths <- russiaDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Tyumen Oblast") %>%
  dplyr::mutate(., Location="Tyumen") %>%
  covid19_tidying(.)

khabarovskDeaths <- russiaDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Khabarovsk Krai") %>%
  dplyr::mutate(., Location="Khabarovsk") %>%
  covid19_tidying(.)

leningradDeaths <- russiaDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Leningrad Oblast") %>%
  dplyr::mutate(., Location="Leningrad") %>%
  covid19_tidying(.)

sverdlovskDeaths <- russiaDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Sverdlovsk Oblast") %>%
  dplyr::mutate(., Location="Sverdlovsk") %>%
  covid19_tidying(.)

tatarstanDeaths <- russiaDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Tatarstan Republic") %>%
  dplyr::mutate(., Location="Tatarstan") %>%
  covid19_tidying(.)

moscowDeaths <- russiaDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Moscow") %>%
  dplyr::mutate(., Location="Moscow") %>%
  covid19_tidying(.)

chelyabinskDeaths <- russiaDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Chelyabinsk Oblast") %>%
  dplyr::mutate(., Location="Chelyabinsk") %>%
  covid19_tidying(.)

irkutskDeaths <- russiaDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Irkutsk Oblast") %>%
  dplyr::mutate(., Location="Irkutsk") %>%
  covid19_tidying(.)

saratovDeaths <- russiaDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Saratov Oblast") %>%
  dplyr::mutate(., Location="Saratov") %>%
  covid19_tidying(.)

kaliningradDeaths <- russiaDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Kaliningrad Oblast") %>%
  dplyr::mutate(., Location="Kaliningrad") %>%
  covid19_tidying(.)

murmanskDeaths <- russiaDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Murmansk Oblast") %>%
  dplyr::mutate(., Location="Murmansk") %>%
  covid19_tidying(.)

krasnoyarskDeaths <- russiaDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Krasnoyarsk Krai") %>%
  dplyr::mutate(., Location="Krasnoyarsk") %>%
  covid19_tidying(.)

novosibirskDeaths <- russiaDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Novosibirsk Oblast") %>%
  dplyr::mutate(., Location="Novosibirsk") %>%
  covid19_tidying(.)

stavropolDeaths <- russiaDeaths %>%
  dplyr::filter(., administrative_area_level_2=="Stavropol Krai") %>%
  dplyr::mutate(., Location="Stavropol") %>%
  covid19_tidying(.)


### South Africa
southAfricaDeaths <- read.csv("../data/downloaded_datasets/death_dynamics_covid19_datahub/ZAF.csv",
                              stringsAsFactors=FALSE)

joubertonDeaths <- dplyr::filter(southAfricaDeaths, administrative_area_level_2 %in%
                                 c("Mpumalanga", "North West")) %>%
  dplyr::group_by(., date) %>%
  dplyr::summarize(., deaths=sum(deaths), confirmed=sum(confirmed)) %>%
  dplyr::mutate(., Location="Jouberton Agincourt") %>%
  ungroup(.) %>%
  covid19_tidying(.)

gautengDeaths <- dplyr::filter(southAfricaDeaths, administrative_area_level_2==
                                 "Gauteng") %>%
  dplyr::mutate(., Location=administrative_area_level_2) %>%
  covid19_tidying(.)

### Argentina
argentinaDeaths <- read.csv("../data/downloaded_datasets/death_dynamics_covid19_datahub/ARG.csv",
                              stringsAsFactors=FALSE)

mugicaDeaths <- dplyr::filter(argentinaDeaths,
                              administrative_area_level_3=="Capital") %>%
  dplyr::mutate(., Location="Barrio Mugica") %>%
  covid19_tidying(.)
madrynDeaths <- dplyr::filter(argentinaDeaths,
                              administrative_area_level_3=="Biedma") %>%
  dplyr::mutate(., Location="Madryn") %>%
  covid19_tidying(.)

### Chile
chileDeaths <- read.csv("../data/downloaded_datasets/death_dynamics_covid19_datahub/CHL.csv",
                              stringsAsFactors=FALSE)

chileDeaths <- dplyr::filter(chileDeaths, administrative_area_level_3 %in%
                                c("La Serena", "Talca", "Santiago")) %>%
  dplyr::group_by(., date) %>%
  dplyr::summarize(., deaths=sum(deaths), confirmed=sum(confirmed)) %>%
  dplyr::mutate(., Location="Chile") %>%
  ungroup(.) %>%
  covid19_tidying(.)

### India
indiaDeaths <- read.csv("../data/downloaded_datasets/death_dynamics_covid19_datahub/IND.csv",
                        stringsAsFactors=FALSE)

delhiDeaths <- dplyr::filter(indiaDeaths,
                             administrative_area_level_2=="Delhi") %>%
  dplyr::mutate(., Location="Delhi") %>%
  covid19_tidying(.)
hyderabadDeaths <- dplyr::filter(indiaDeaths,
                                 administrative_area_level_2=="Telagana") %>%
  dplyr::mutate(., Location="Hyderabad") %>%
  covid19_tidying(.)
kashmirDeaths <- dplyr::filter(indiaDeaths,
                               administrative_area_level_2=="Jammu And Kashmir") %>%
  dplyr::mutate(., Location="Kashmir") %>%
  covid19_tidying(.)
odishaDeaths <- dplyr::filter(indiaDeaths,
                              administrative_area_level_2=="Orissa") %>%
  dplyr::mutate(., Location="Odisha") %>%
  covid19_tidying(.)
puducherryDeaths <- dplyr::filter(indiaDeaths,
                                  administrative_area_level_2=="Puducherry") %>%
  dplyr::mutate(., Location="Puducherry") %>%
  covid19_tidying(.)
puneDeaths <- dplyr::filter(indiaDeaths,
                            administrative_area_level_2=="Maharashtra") %>%
  dplyr::mutate(., Location="Pune") %>%
  covid19_tidying(.)
tamilDeaths <- dplyr::filter(indiaDeaths,
                             administrative_area_level_2=="Tamil Nadu") %>%
  dplyr::mutate(., Location="Tamil Nadu") %>%
  covid19_tidying(.)
bengalDeaths <- dplyr::filter(indiaDeaths,
                              administrative_area_level_2=="West Bengal") %>%
  dplyr::mutate(., Location="West Bengal") %>%
  covid19_tidying(.)
indiaNationalDeaths <- dplyr::filter(indiaDeaths, administrative_area_level_2=="") %>%
  dplyr::mutate(., Location="India") %>%
  covid19_tidying(.)

### Indonesia
indonesiaDeaths <- read.csv("../data/downloaded_datasets/death_dynamics_covid19_datahub/IDN.csv",
                        stringsAsFactors=FALSE) %>%
  dplyr::mutate(., Location="Indonesia") %>%
  covid19_tidying(.)

### Spain
spainDeaths <- read.csv("../data/downloaded_datasets/death_dynamics_covid19_datahub/ESP.csv",
                        stringsAsFactors=FALSE)
cantabriaDeaths <- dplyr::filter(spainDeaths, administrative_area_level_2=="Cantabria"
  & administrative_area_level_3=="") %>%
  dplyr::mutate(., Location="Cantabria") %>%
  covid19_tidying(.)

### Tunisia
tunisiaDeaths <- read.csv("../data/downloaded_datasets/death_dynamics_covid19_datahub/TUN.csv",
                        stringsAsFactors=FALSE) %>%
  dplyr::mutate(., Location="Tunisia") %>%
  covid19_tidying(.)

#######################
# Put data together
#######################

allDeaths <- rbind(countriesDeaths, franceDeaths, usaStatesDeaths,
  usaCountyDeaths, englandDeaths, scotlandDeaths, walesDeaths,
  northernIrelandDeaths, genevaDeaths, neuchatelDeaths, fribourgDeaths,
  ticinoDeaths, freiburgDeaths, reutlingenDeaths, osnabruckDeaths,
  aachenDeaths, magdeburgDeaths, munichDeaths, gangeltDeaths,
  mitteDeaths, kupferzellDeaths, tirschenreuthDeaths, straubingDeaths,
  chemnitzDeaths, vorpommernDeaths,
  ischglDeaths, wachauDeaths, trentoDeaths, voDeaths, addaDeaths, caldariDeaths,
  wuhanDeaths, petersburgDeaths, joubertonDeaths, gautengDeaths,
  mugicaDeaths, madrynDeaths, chileDeaths, delhiDeaths,
  hyderabadDeaths, kashmirDeaths, odishaDeaths,
  puducherryDeaths, puneDeaths, tamilDeaths, bengalDeaths,
  indonesiaDeaths, cantabriaDeaths, tunisiaDeaths, connecticutDeaths,
  petersburgDeaths, tyumenDeaths, khabarovskDeaths, leningradDeaths,
  sverdlovskDeaths, tatarstanDeaths, moscowDeaths, chelyabinskDeaths,
  irkutskDeaths, saratovDeaths, kaliningradDeaths, murmanskDeaths,
  krasnoyarskDeaths, novosibirskDeaths, stavropolDeaths)

write.csv(allDeaths, "../data/raw_data/death_dynamics.csv", row.names=FALSE)

