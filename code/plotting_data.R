library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)


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

ggsave("../data/figure/seroprevalence_dates.png", datesPlot,
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

