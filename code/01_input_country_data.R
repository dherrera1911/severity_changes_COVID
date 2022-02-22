#################################
#################################
#################################
#  
#  This script puts together the information
#  reported in several studies, and by official
#  organisms, of estimates of SARS-CoV-2 infections
#  and age-stratified hospitalizations, ICU admissions
#  and deaths for different locations. The sources are
#  indicated in the script and in the manuscript
#  
#  Nana Owusu-Boaitey and Daniel Herrera-Esposito gathered the
#  data from the cited sources. Daniel wrote the code and did
#  the preprocessing of the data.
#
#  2022
#  
#  Direct any questions to dherrera@fcien.edu.uy
#  
#  
#################################
#################################
#################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(wpp2019)
library(rriskDistributions)
source("./functions_auxiliary.R")
data(popM)
data(popF)

# Function to find IFR of each age
levinIFR <- function(age){return(10^(-3.27+0.0524*age))}
# Function to find ISR for each age
dhISR <- function(age){return(exp(-7.3+0.076*age))}

################
################
## Spain
################
################

############
### TIME 1
############
age_Spain <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70+")

### POPULATION
population_Spain <- extract_country_population(popM=popM, popF=popF,
                                               countryName="Spain",
                                               ageBins=age_Spain)

### SEROPREVALENCE
# Spain seroprevalence Pollán et al
# I take adjusted seroprevalence from Levin et al
seroprev_Spain <- c(2.4, 2.9, 3.4, 3.0, 3.6, 3.7, 3.4, 3.1)
seroprevL_Spain <- c(1.5, 2.3, 2.7, 2.4, 3.0, 3.1, 2.8, 2.5)
seroprevH_Spain <- c(3.4, 3.5, 4.1, 3.6, 4.5, 4.7, 4.5, 5.0)

### OUTCOMES
# https://cnecovid.isciii.es/covid19/#documentaci%C3%B3n-y-datos
spainDailyOutcome <- read.csv("../data/downloaded_datasets/spain/casos_hosp_uci_def_sexo_edad_provres.csv",
                    stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., fecha=lubridate::date(fecha))

spainDeaths <- dplyr::filter(spainDailyOutcome, fecha <= "2020-05-04" &
                                  grupo_edad != "NC") %>%
  group_by(., grupo_edad) %>%
  summarize(., Deaths=sum(num_def)) %>%
  dplyr::mutate(., proportion=Deaths, age=grupo_edad) %>%
  change_demography_bins(., age_Spain) %>%
  dplyr::mutate(., Deaths=proportion)

spainICU <- dplyr::filter(spainDailyOutcome, fecha <= "2020-05-04" &
                                  grupo_edad != "NC") %>%
  group_by(., grupo_edad) %>%
  summarize(., ICU=sum(num_uci)) %>%
  dplyr::mutate(., proportion=ICU, age=grupo_edad) %>%
  change_demography_bins(., age_Spain) %>%
  dplyr::mutate(., ICU=proportion)

spainHosp <- dplyr::filter(spainDailyOutcome, fecha <= "2020-05-04" &
                                  grupo_edad != "NC") %>%
  group_by(., grupo_edad) %>%
  summarize(., Hospitalized=sum(num_hosp)) %>%
  dplyr::mutate(., proportion=Hospitalized, age=grupo_edad) %>%
  change_demography_bins(., age_Spain) %>%
  dplyr::mutate(., Hospitalized=proportion)

Hospitalized_Spain <- spainHosp$Hospitalized
ICU_Spain <- spainICU$ICU
deaths_Spain <- spainDeaths$Deaths

### FULL TABLE
Spain_w1 <- data.frame(Age=age_Spain,
                    Population=population_Spain,
                    Prevalence=seroprev_Spain,
                    PrevalenceL=seroprevL_Spain,
                    PrevalenceH=seroprevH_Spain,
                    Prevalence_unvacc=seroprev_Spain,
                    PrevalenceL_unvacc=seroprevL_Spain,
                    PrevalenceH_unvacc=seroprevH_Spain,
                    Severe=Hospitalized_Spain,
                    Critical=ICU_Spain,
                    Deaths=deaths_Spain,
                    Deaths_unvacc=deaths_Spain,
                    onePlusVaccine=0,
                    twoPlusVaccine=0,
                    Type="Seroprevalence",
                    Location="Spain",
                    EndPointOutcome="2020-05-04",
                    MidPointCases="2020-05-11",
                    Period=1)

############
### TIME 2
############
age_Spain <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70+")
seroprevAges_Spain <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
                     "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                     "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90+")

### POPULATION
population_Spain <- extract_country_population(popM=popM, popF=popF,
                                               countryName="Spain",
                                               ageBins=seroprevAges_Spain)
population2_Spain <- extract_country_population(popM=popM, popF=popF,
                                               countryName="Spain",
                                               ageBins=age_Spain)

### SEROPREVALENCE
# Spain seroprevalence ENE round 4: https://www.mscbs.gob.es/gabinetePrensa/notaPrensa/pdf/15.12151220163348113.pdf
seroprev_Spain <- c(4.3, 5.9, 6.8, 6.1, 8.1, 7.6, 7.0, 6.7, 6.6, 7.3, 8.5, 7.3,
  7.8, 7.0, 7.6, 7.0, 6.9, 6.6, 6.5)
seroprevL_Spain <- c(3.1, 4.6, 5.5, 4.9, 6.6, 6.0, 5.6, 5.6, 5.7, 6.2, 7.5,
  6.3, 6.7, 5.8, 6.4, 5.6, 5.2, 4.5, 3.9)
seroprevH_Spain <- c(5.9, 7.5, 8.2, 7.6, 9.8, 9.6, 8.7, 8.0, 7.7, 8.5, 9.7,
  8.5, 9.1, 8.3, 9.1, 8.7, 9.2, 9.5, 10.8)

# redistribute bins of seroprevalence to match outcome bins
inds1 <- seq(1,14,2)
inds2 <- seq(2,14,2)
estimatedInf <- seroprev_Spain*population_Spain
seroprev_Spain <- c(estimatedInf[inds1] + estimatedInf[inds2],
                    sum(estimatedInf[15:19])) / population2_Spain
seroprev_Spain <- round(seroprev_Spain, digits=1)
estimatedLInf <- seroprevL_Spain*population_Spain
seroprevL_Spain <- c(estimatedLInf[inds1] + estimatedLInf[inds2],
                     sum(estimatedLInf[15:19])) / population2_Spain
seroprevL_Spain <- round(seroprevL_Spain, digits=1)
estimatedHInf <- seroprevH_Spain*population_Spain
seroprevH_Spain <- c(estimatedHInf[inds1] + estimatedHInf[inds2],
                     sum(estimatedHInf[15:19])) / population2_Spain
seroprevH_Spain <- round(seroprevH_Spain, digits=1)

### OUTCOMES
spainDailyOutcome <- read.csv("../data/downloaded_datasets/spain/casos_hosp_uci_def_sexo_edad_provres.csv",
                    stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., fecha=lubridate::date(fecha))

spainDeaths <- dplyr::filter(spainDailyOutcome, fecha <= "2020-11-22" &
                                  grupo_edad != "NC") %>%
  group_by(., grupo_edad) %>%
  summarize(., Deaths=sum(num_def)) %>%
  dplyr::mutate(., proportion=Deaths, age=grupo_edad) %>%
  change_demography_bins(., age_Spain) %>%
  dplyr::mutate(., Deaths=proportion)

spainICU <- dplyr::filter(spainDailyOutcome, fecha <= "2020-11-22" &
                                  grupo_edad != "NC") %>%
  group_by(., grupo_edad) %>%
  summarize(., ICU=sum(num_uci)) %>%
  dplyr::mutate(., proportion=ICU, age=grupo_edad) %>%
  change_demography_bins(., age_Spain) %>%
  dplyr::mutate(., ICU=proportion)

spainHosp <- dplyr::filter(spainDailyOutcome, fecha <= "2020-11-22" &
                                  grupo_edad != "NC") %>%
  group_by(., grupo_edad) %>%
  summarize(., Hospitalized=sum(num_hosp)) %>%
  dplyr::mutate(., proportion=Hospitalized, age=grupo_edad) %>%
  change_demography_bins(., age_Spain) %>%
  dplyr::mutate(., Hospitalized=proportion)

Hospitalized_Spain <- spainHosp$Hospitalized
ICU_Spain <- spainICU$ICU
deaths_Spain <- spainDeaths$Deaths

### FULL TABLE
Spain_w2 <- data.frame(Age=age_Spain,
                    Population=population2_Spain,
                    Prevalence=seroprev_Spain,
                    PrevalenceL=seroprevL_Spain,
                    PrevalenceH=seroprevH_Spain,
                    Prevalence_unvacc=seroprev_Spain,
                    PrevalenceL_unvacc=seroprevL_Spain,
                    PrevalenceH_unvacc=seroprevH_Spain,
                    Severe=Hospitalized_Spain,
                    Critical=ICU_Spain,
                    Deaths=deaths_Spain,
                    Deaths_unvacc=deaths_Spain,
                    onePlusVaccine=0,
                    twoPlusVaccine=0,
                    Type="Seroprevalence",
                    Location="Spain",
                    EndPointOutcome="2020-11-29",
                    MidPointCases="2020-11-22",
                    Period=2)


###########################
# Fulton, de Kalb, Georgia
###########################

#######
### TIME 1
#######

### POPULATION
# https://data.census.gov/cedsci/table?tid=ACSST5Y2019.S0101&g=1600000US1304000
atlantaPopAge <- c("0-5", "5-9", "10-14", "15-19", "20-24", "25-29",
                   "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                   "60-64", "65-69", "70-74", "75-79", "80-84", "85+")
atlantaPop <- c(26577, 27559, 22290, 32902, 48054, 55760, 46168,
                       36429, 31238, 29916, 27719, 26051, 21839, 18172,
                       14829, 9696, 6821, 6780)
# 10% of Atlanta is in DeKalb, which is not included in Fulton
# country numbers. adjust population
atlantaPop <- round(atlantaPop*0.9)
atlantaPop <- c(sum(atlantaPop[5:10])+round(atlantaPop[4]*2/5),
                sum(atlantaPop[11:13]),
                sum(atlantaPop[14:18]))

### SEROPREVALENCE
# Seroprev: Estimated Community Seroprevalence of SARS-CoV-2 Antibodies —
# Two Georgia Counties, Biggs et al
# remove ranges 0-17 which had 0 positive serology
# samples, following Levin et al
# Note: Seroprevalence estimate is from Atlanta at Fulton + Kalb counties,
# but since outcomes are only Fulton, we use the populaiton of Fulton alone to
# estimate number of infections
age_Atlanta_seroprev <- c("18-49", "50-64", "65+")
seroprev_Atlanta <- c(3.3, 4.9, 0.7)
seroprevL_Atlanta <- c(1.6, 1.8, 0.1)
seroprevH_Atlanta <- c(6.4, 12.9, 4.5)
### Correct for test characteristics
sensitivity <- 0.932
specificity <- 0.99
factor <- 1/(sensitivity + specificity - 1)
#seroprev_Atlanta <- signif((seroprev_Atlanta/100 + specificity - 1)*factor*100, digits=3)
### CHECK CORRECTION
seroprev_Atlanta <- signif(seroprev_Atlanta*factor, digits=3)
seroprevL_Atlanta <- signif(seroprevL_Atlanta*factor, digits=3)
seroprevH_Atlanta <- signif(seroprevH_Atlanta*factor, digits=3)



### OUTCOMES
# Outcome age distribution "Characteristics and Risk Factors for Hospitalization
# and Mortality among Persons with COVID-19 in Atlanta Metropolitan Area"
# Chishinga et al. Paper data are until 31 May. 
age_Atlanta_outcome <- c("0-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75+")
Hospitalized_Atlanta <- c(18, 62, 79, 115, 172, 180, 260)
ICU_Atlanta <- c(4, 13, 14, 33, 44, 57, 74)
deaths_Atlanta <- c(2, 3, 2, 16, 33, 69, 170)

# Hospital data https://www.fultoncountyga.gov/covid-19/epidemiology-reports
casesAtlantaFulton <- 1398 + round(0.425*399) # atlanta cases + portion of unknown
Hospitalized_AtlantaFulton <- round(0.179 * casesAtlantaFulton)
ICU_AtlantaFulton <- round(0.049 * casesAtlantaFulton)
deaths_AtlantaFulton <- round(0.040 * casesAtlantaFulton)

# use distribution of Chishinga with absolute numbers from governmnet report
Hospitalized_Atlanta <- round(Hospitalized_Atlanta*Hospitalized_AtlantaFulton/sum(Hospitalized_Atlanta))
ICU_Atlanta <- round(ICU_Atlanta*ICU_AtlantaFulton/sum(ICU_Atlanta))
deaths_Atlanta <- round(deaths_Atlanta*deaths_AtlantaFulton/sum(deaths_Atlanta))

# redistribute to match serology ages
Hospitalized_Atlanta <- c(sum(Hospitalized_Atlanta[2:3]),
                          sum(Hospitalized_Atlanta[4:5]),
                          sum(Hospitalized_Atlanta[6:7]))
ICU_Atlanta <- c(sum(ICU_Atlanta[2:3]),
                 sum(ICU_Atlanta[4:5]),
                 sum(ICU_Atlanta[6:7]))
deaths_Atlanta <- c(sum(deaths_Atlanta[2:3]),
                    sum(deaths_Atlanta[4:5]),
                    sum(deaths_Atlanta[6:7]))

Atlanta_w1 <- data.frame(Age=age_Atlanta_seroprev,
                      Population=atlantaPop,
                      Prevalence=seroprev_Atlanta,
                      PrevalenceL=seroprevL_Atlanta,
                      PrevalenceH=seroprevH_Atlanta,
                      Prevalence_unvacc=seroprev_Atlanta,
                      PrevalenceL_unvacc=seroprevL_Atlanta,
                      PrevalenceH_unvacc=seroprevH_Atlanta,
                      Severe=Hospitalized_Atlanta,
                      Critical=ICU_Atlanta,
                      Deaths=deaths_Atlanta,
                      Deaths_unvacc=deaths_Atlanta,
                      onePlusVaccine=0,
                      twoPlusVaccine=0,
                      Type="Seroprevalence",
                      Location="Atlanta",
                      EndPointOutcome="2020-05-04",
                      MidPointCases="2020-05-03",
                      Period=1)

#######
### TIME 2
#######

### POPULATION
atlantaPopAge <- c("0-5", "5-9", "10-14", "15-19", "20-24", "25-29",
                   "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                   "60-64", "65-69", "70-74", "75-79", "80-84", "85+")
atlantaPop <- c(26577, 27559, 22290, 32902, 48054, 55760, 46168,
                       36429, 31238, 29916, 27719, 26051, 21839, 18172,
                       14829, 9696, 6821, 6780)
# 10% of Atlanta is in DeKalb, which is not included in Fulton
# country numbers. adjust population
atlantaPop <- round(atlantaPop*0.9)
atlantaPop <- c(sum(atlantaPop[7:8]),
                sum(atlantaPop[9:10]),
                sum(atlantaPop[11:12]),
                sum(atlantaPop[13:14]),
                sum(atlantaPop[15:18]))

### SEROPREVALENCE
# Seroprevalence & pop: SARS-CoV-2 Cumulative Incidence and Wave Seroprevalence: Results From a Statewide PopulationBased Serosurvey in Atlanta
# Data is taken from supplementary Table 1, presenting for Fulton-DeKalb specifically
age_Atlanta_seroprev <- c("18-34", "35-44", "45-54", "55-64", "65+")
seroprev_Atlanta <- c(8.08, 11.48, 5.68, 7.10, 5.78)
seroprevL_Atlanta <- c(3.74, 4.83, 2.05, 2.50, 2.24)
seroprevH_Atlanta <- c(16.57, 24.91, 14.76, 18.55, 14.12)
# remove first seroprev value, to match outcome ages
seroprev_Atlanta <- seroprev_Atlanta
seroprevL_Atlanta <- seroprevL_Atlanta
seroprevH_Atlanta <- seroprevH_Atlanta

### OUTCOMES
# https://www.fultoncountyga.gov/covid-19/epidemiology-reports 17 November 2020
age_Atlanta_outcomes <- c("30-39", "40-49", "50-59", "60-69", "70+")
deaths_Atlanta <- c(14, 26, 55, 117, 435)

### FULL TABLE
Atlanta_w2 <- data.frame(Age=age_Atlanta_seroprev,
                  Population=atlantaPop,
                  Prevalence=seroprev_Atlanta,
                  PrevalenceL=seroprevL_Atlanta,
                  PrevalenceH=seroprevH_Atlanta,
                  Prevalence_unvacc=seroprev_Atlanta,
                  PrevalenceL_unvacc=seroprevL_Atlanta,
                  PrevalenceH_unvacc=seroprevH_Atlanta,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Atlanta,
                  Deaths_unvacc=deaths_Atlanta,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Atlanta",
                  EndPointOutcome="2020-11-17",
                  MidPointCases="2020-11-16",
                  Period=2)

# Data for all Georgia, instead of Fulton
# Deaths, Hospital: https://dph.georgia.gov/covid-19-daily-status-report
# https://ga-covid19.ondemand.sas.com/docs/ga_covid_data.zip
# https://agetrends.s3.amazonaws.com/AgeTrends_20210923.html#hospitalizations
# Deaths https://data.chhs.ca.gov/dataset/covid-19-time-series-metrics-by-county-and-state
#deathDataGE <- read.csv("../data/downloaded_datasets/georgia/epicurve_age_grp_symptom_date.csv",
#                       stringsAsFactors=FALSE) %>%
#  as_tibble(.) %>%
#  dplyr::filter(., (lubridate::date(symptom.date)<="2020-12-08")) %>%
#  dplyr::group_by(., age_group) %>%
#  dplyr::summarize(., deaths=sum(deaths), hospitalized=sum(confirmed_case_hospitalization)) %>%
#  dplyr::filter(., age_group!="Unknown years")
#
## Redistribute outcomes to match seroprevalence bins
#geDeaths <- deathDataGE$deaths
## relative IFR is the same if the ranges are always the same
#relIFR <- levinIFR(seq(30,39))/sum(levinIFR(seq(30,39)))
#relIFR <- c(sum(relIFR[1:4]), sum(relIFR[5:10]))
#deaths_Georgia <- c(geDeaths[5]+floor(geDeaths[6]*relIFR[1]),
#                    ceiling(geDeaths[6]*relIFR[2])+floor(geDeaths[7]*relIFR[1]),
#                    ceiling(geDeaths[7]*relIFR[2])+floor(geDeaths[8]*relIFR[1]),
#                    ceiling(geDeaths[8]*relIFR[2])+floor(geDeaths[9]*relIFR[1]),
#                    ceiling(geDeaths[9]*relIFR[2])+sum(geDeaths[10:11]))
#
#geHosp <- deathDataGE$hospitalized
## relative IFR is the same if the ranges are always the same
#relISR <- dhISR(seq(30,39))/sum(dhISR(seq(30,39)))
#relISR <- c(sum(relISR[1:4]), sum(relISR[5:10]))
#hospitalized_Georgia <- c(geHosp[5]+floor(geHosp[6]*relISR[1]),
#                    ceiling(geHosp[6]*relISR[2])+floor(geHosp[7]*relISR[1]),
#                    ceiling(geHosp[7]*relISR[2])+floor(geHosp[8]*relISR[1]),
#                    ceiling(geHosp[8]*relISR[2])+floor(geHosp[9]*relISR[1]),
#                    ceiling(geHosp[9]*relISR[2])+sum(geHosp[10:11]))

######################
# France
######################

#########
### TIME 1
#########
age_France_seroprev <- c("15-17", "18-24", "25-34", "35-44", "45-54",
                         "55-64", "65-74", "75+")

### POPULATION
popAge <- c("0-14", "15-19", "20-24", age_France_seroprev[3:8])
population_France <- extract_country_population(popM=popM, popF=popF,
                                               countryName="France",
                                               ageBins=popAge)
population_France <- c(population_France[2]*3/5,
                       population_France[2]*2/5+population_France[3],
                       population_France[4:9])

### SEROPREVALENCE
# Seroprev: Trends in social exposure to SARS-Cov-2 in France. Evidence from the
# national socio-epidemiological cohort – EPICOV
seroprev_France <- c(4.5, 4.8, 5.0, 8.3, 4.9, 4.8, 1.8, 0.7)
seroprevL_France <- c(2.2, 3.0, 3.7, 6.7, 3.9, 3.3, 1.2, 0.4)
seroprevH_France <- c(8.9, 7.6, 6.7, 10.4, 6.2, 7.1, 2.7, 1.4)

### OUTCOMES
#https://www.santepubliquefrance.fr/maladies-et-traumatismes/maladies-et-infections-respiratoires/infection-a-coronavirus/documents/bulletin-national/covid-19-point-epidemiologique-du-21-mai-2020
age_France_outcome <- c("0-14", "15-44", "45-64", "65-74", "75+")
hospital_current_France <- c(72, 840, 3545, 3897, 10536)
discharged_France <- c(726, 9085, 19378, 11755, 20338)
death_hospital_France <- c(3, 174, 1787, 3108, 12416)
deaths_outside_hospital_France <- 10245
# allocate ooh death to total deaths and hospital
deaths_France <- round(death_hospital_France +
  c(0,0,0,deaths_outside_hospital_France*0.25,
    deaths_outside_hospital_France*0.75))
hospitalized_France <- round(hospital_current_France + discharged_France +
  c(0,0,0,deaths_outside_hospital_France*0.25,
    deaths_outside_hospital_France*0.75))
# Redistribute deaths to match seroprevalence bins
relIFR1 <- levinIFR(seq(15,44))/sum(levinIFR(seq(15,44)))
relIFR1 <- c(sum(relIFR1[1:3]), sum(relIFR1[4:10]),
             sum(relIFR1[11:20]), sum(relIFR1[21:30]))
relIFR2 <- levinIFR(seq(45,64))/sum(levinIFR(seq(45,64)))
relIFR2 <- c(sum(relIFR2[1:10]), sum(relIFR2[11:20]))
deaths_France <- c(round(deaths_France[2]*relIFR1[1]),
                   round(deaths_France[2]*relIFR1[2]),
                   round(deaths_France[2]*relIFR1[3]),
                   round(deaths_France[2]*relIFR1[4]),
                   round(deaths_France[3]*relIFR2[1]),
                   round(deaths_France[3]*relIFR2[2]),
                   deaths_France[4:5])
# Redistribute hospitalized to match seroprevalence bins
relISR1 <- dhISR(seq(15,44))/sum(dhISR(seq(15,44)))
relISR1 <- c(sum(relISR1[1:3]), sum(relISR1[4:10]),
             sum(relISR1[11:20]), sum(relISR1[21:30]))
relISR2 <- dhISR(seq(45,64))/sum(dhISR(seq(45,64)))
relISR2 <- c(sum(relISR2[1:10]), sum(relISR2[11:20]))
hospitalized_France <- c(round(hospitalized_France[2]*relISR1[1]),
                   round(hospitalized_France[2]*relISR1[2]),
                   round(hospitalized_France[2]*relISR1[3]),
                   round(hospitalized_France[2]*relISR1[4]),
                   round(hospitalized_France[3]*relISR2[1]),
                   round(hospitalized_France[3]*relISR2[2]),
                   hospitalized_France[4:5])


### FULL TABLE
France_w1 <- data.frame(Age=age_France_seroprev,
                  Population=population_France,
                  Prevalence=seroprev_France,
                  PrevalenceL=seroprevL_France,
                  PrevalenceH=seroprevH_France,
                  Prevalence_unvacc=seroprev_France,
                  PrevalenceL_unvacc=seroprevL_France,
                  PrevalenceH_unvacc=seroprevH_France,
                  Severe=hospitalized_France,
                  Critical=NA,
                  Deaths=deaths_France,
                  Deaths_unvacc=deaths_France,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="France",
                  EndPointOutcome="2020-05-21",
                  MidPointCases="2020-05-21",
                  Period=1)


############
### TIME 2
############
age_France_seroprev <- c("15-17", "18-24", "25-34", "35-44", "45-54",
                         "55-64", "65-74", "75+")

### SEROPREVALENCE
# Seroprev: Trends in social exposure to SARS-Cov-2 in France. Evidence from the
# national socio-epidemiological cohort – EPICOV
seroprev_France <- c(9.8, 10.0, 7.2, 6.5, 6.5, 5.3, 4.3, 3.7)
seroprevL_France <- c(7.8, 8.6, 6.3, 5.8, 5.9, 4.8, 3.8, 2.9)
seroprevH_France <- c(12.2, 11.5, 8.3, 7.4, 7.2, 5.8, 4.9, 4.7)

### OUTCOMES
# Deaths, hospitalizations: https://www.santepubliquefrance.fr/maladies-et-traumatismes/maladies-et-infections-respiratoires/infection-a-coronavirus/documents/bulletin-national/covid-19-point-epidemiologique-du-26-novembre-2020
age_France_outcome <- c("0-14", "15-44", "45-64", "65-74", "75+")
hospital_current_France <- c(88, 1064, 4696, 6408, 18143)
discharged_France <- c(2165, 21612, 43063, 30537, 56233)
death_hospital_France <- c(3, 301, 3046, 5793, 25068)
deaths_outside_hospital_France <- 15838
# allocate ooh death to total deaths and hospital
deaths_France <- round(death_hospital_France +
  c(0,0,0,deaths_outside_hospital_France*0.25,
    deaths_outside_hospital_France*0.75))
hospitalized_France <- round(hospital_current_France + discharged_France +
  c(0,0,0,deaths_outside_hospital_France*0.25,
    deaths_outside_hospital_France*0.75))
# Redistribute deaths to match seroprevalence bins
relIFR1 <- levinIFR(seq(15,44))/sum(levinIFR(seq(15,44)))
relIFR1 <- c(sum(relIFR1[1:3]), sum(relIFR1[4:10]),
             sum(relIFR1[11:20]), sum(relIFR1[21:30]))
relIFR2 <- levinIFR(seq(45,64))/sum(levinIFR(seq(45,64)))
relIFR2 <- c(sum(relIFR2[1:10]), sum(relIFR2[11:20]))
deaths_France <- c(round(deaths_France[2]*relIFR1[1]),
                   round(deaths_France[2]*relIFR1[2]),
                   round(deaths_France[2]*relIFR1[3]),
                   round(deaths_France[2]*relIFR1[4]),
                   round(deaths_France[3]*relIFR2[1]),
                   round(deaths_France[3]*relIFR2[2]),
                   deaths_France[4:5])
# Redistribute hospitalized to match seroprevalence bins
relISR1 <- dhISR(seq(15,44))/sum(dhISR(seq(15,44)))
relISR1 <- c(sum(relISR1[1:3]), sum(relISR1[4:10]),
             sum(relISR1[11:20]), sum(relISR1[21:30]))
relISR2 <- dhISR(seq(45,64))/sum(dhISR(seq(45,64)))
relISR2 <- c(sum(relISR2[1:10]), sum(relISR2[11:20]))
hospitalized_France <- c(round(hospitalized_France[2]*relISR1[1]),
                   round(hospitalized_France[2]*relISR1[2]),
                   round(hospitalized_France[2]*relISR1[3]),
                   round(hospitalized_France[2]*relISR1[4]),
                   round(hospitalized_France[3]*relISR2[1]),
                   round(hospitalized_France[3]*relISR2[2]),
                   hospitalized_France[4:5])


### FULL TABLE
France_w2 <- data.frame(Age=age_France_seroprev,
                  Population=population_France,
                  Prevalence=seroprev_France,
                  PrevalenceL=seroprevL_France,
                  PrevalenceH=seroprevH_France,
                  Prevalence_unvacc=seroprev_France,
                  PrevalenceL_unvacc=seroprevL_France,
                  PrevalenceH_unvacc=seroprevH_France,
                  Severe=hospitalized_France,
                  Critical=NA,
                  Deaths=deaths_France,
                  Deaths_unvacc=deaths_France,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="France",
                  EndPointOutcome="2020-11-02",
                  MidPointCases="2020-11-02",
                  Period=2)


######################
# Geneva, Switzerland
######################

#######
### TIME 1
#######
age_Geneva_seroprev <- c("5-9", "10-19", "20-49", "50-64", "65+")

### POPULATION
population_Geneva <- c(26466, 53180, 219440, 98528, 83574)

### SEROPREVALENCE
# Seroprevalence data from:
# Serology-informed estimates of SARS-CoV-2 infection
# fatality risk in Geneva, Switzerland
infected_Geneva <- c(1200, 6100, 28800, 10300, 5700)
infectedL_Geneva <- c(400, 3900, 21400, 7200, 3200)
infectedH_Geneva <- c(2400, 8800, 37700, 13900, 8800)
seroprev_Geneva <- signif(infected_Geneva/population_Geneva*100, digits=3)
seroprevL_Geneva <- signif(infectedL_Geneva/population_Geneva*100, digits=3)
seroprevH_Geneva <- signif(infectedH_Geneva/population_Geneva*100, digits=3)

### OUTCOMES
# Hospital and death data from:
# https://www.covid19.admin.ch/en/weekly-report/hosp?geoView=table
hospDataGE <- read.csv("../data/downloaded_datasets/switzerland/COVID19Hosp_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="GE" & datum_dboardformated=="2020-17" &
                altersklasse_covid19!="Unbekannt")

deathDataGE <- read.csv("../data/downloaded_datasets/switzerland/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="GE" & datum_dboardformated=="2020-17" &
                altersklasse_covid19!="Unbekannt")

Hospitalized_Geneva <- hospDataGE$sumTotal
deaths_Geneva <- deathDataGE$sumTotal

# rearrange hospitalizations and deaths to match seroprevalence
Hospitalized_Geneva <- c(Hospitalized_Geneva[2],
                         sum(Hospitalized_Geneva[3:5]),
                         Hospitalized_Geneva[6]+Hospitalized_Geneva[7]/2,
                         Hospitalized_Geneva[7]/2+sum(Hospitalized_Geneva[8:9]))
deaths_Geneva <- c(deaths_Geneva[2],
                         sum(deaths_Geneva[3:5]),
                         deaths_Geneva[6]+deaths_Geneva[7]/2,
                         deaths_Geneva[7]/2+sum(deaths_Geneva[8:9]))

### FULL TABLE
Geneva_w1 <- data.frame(Age=age_Geneva_seroprev[-1],
                  Population=population_Geneva[-1],
                  Prevalence=seroprev_Geneva[-1],
                  PrevalenceL=seroprevL_Geneva[-1],
                  PrevalenceH=seroprevH_Geneva[-1],
                  Prevalence_unvacc=seroprev_Geneva[-1],
                  PrevalenceL_unvacc=seroprevL_Geneva[-1],
                  PrevalenceH_unvacc=seroprevH_Geneva[-1],
                  Severe=Hospitalized_Geneva,
                  Critical=NA,
                  Deaths=deaths_Geneva,
                  Deaths_unvacc=deaths_Geneva,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Geneva",
                  EndPointOutcome="2020-05-10",
                  MidPointCases="2020-04-22",
                  Period=1)


#######
### TIME 2
#######
age_Geneva_seroprev <- c("0-5", "6-11", "12-17", "18-24", "25-34",
                         "35-49", "50-64", "65-74", "75+")
age_Geneva_outcomes <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
                         "60-69", "70-79", "80+")

### SEROPREVALENCE
seroprev_Geneva <- c(14.9, 22.8, 23.6, 25.4, 25.9, 23.6, 21.2, 14.9, 9.3)
seroprevL_Geneva <- c(10.7, 18.7, 19.6, 19.8, 21.8, 20.8, 18.1, 10.9, 5.5)
seroprevH_Geneva <- c(19.6, 27.1, 28.0, 31.5, 30.2, 26.4, 24.6, 19.5, 13.7)

# Adjust seroprevalence to match outcomes
seroprev_Geneva <- c(mean(seroprev_Geneva[1:2]), seroprev_Geneva[3:8],
  mean(seroprev_Geneva[8:9]), seroprev_Geneva[9])
seroprevL_Geneva <- c(mean(seroprevL_Geneva[1:2]), seroprevL_Geneva[3:8],
  mean(seroprevL_Geneva[8:9]), seroprevL_Geneva[9])
seroprevH_Geneva <- c(mean(seroprevH_Geneva[1:2]), seroprevH_Geneva[3:8],
  mean(seroprevH_Geneva[8:9]), seroprevH_Geneva[9])

### OUTCOMES
# Hospital and death data from:
# https://www.covid19.admin.ch/en/weekly-report/hosp?geoView=table
hospDataGE <- read.csv("../data/downloaded_datasets/switzerland/COVID19Hosp_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="GE" & datum_dboardformated=="2020-50" &
                altersklasse_covid19!="Unbekannt")
deathDataGE <- read.csv("../data/downloaded_datasets/switzerland/nuevos_files/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="GE" & datum_dboardformated=="2020-50" &
                altersklasse_covid19!="Unbekannt")
hospitalized_Geneva <- hospDataGE$sumTotal
deaths_Geneva <- deathDataGE$sumTotal
population_Geneva <- deathDataGE$pop

### FULL TABLE

Geneva_w2 <- data.frame(Age=age_Geneva_outcomes,
                  Population=population_Geneva,
                  Prevalence=seroprev_Geneva,
                  PrevalenceL=seroprevL_Geneva,
                  PrevalenceH=seroprevH_Geneva,
                  Prevalence_unvacc=seroprev_Geneva,
                  PrevalenceL_unvacc=seroprevL_Geneva,
                  PrevalenceH_unvacc=seroprevH_Geneva,
                  Severe=hospitalized_Geneva,
                  Critical=NA,
                  Deaths=deaths_Geneva,
                  Deaths_unvacc=deaths_Geneva,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Geneva",
                  EndPointOutcome="2020-12-12",
                  MidPointCases="2020-12-08",
                  Period=2)

#######
### TIME 3
#######
age_Geneva_seroprev <- c("0-5", "6-11", "12-17", "18-24", "25-34",
                         "35-49", "50-64", "65-74", "75+")
age_Geneva_outcomes <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
                         "60-69", "70-79", "80+")

### SEROPREVALENCE
# https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2021.26.43.2100830
# Infectious origin antibodies
seroprev_Geneva <- c(20.8, 31.4, 37.7, 41.8, 31.9, 32.2, 29.8, 22.5, 16.2)
seroprevL_Geneva <- c(15.5, 26.7, 32.5, 36.3, 27.8, 28.7, 26.5, 17.0, 11.8)
seroprevH_Geneva <- c(26.7, 36.4, 43.1, 47.5, 36.2, 35.9, 33.4, 28.4, 21.1)

# All antibodies
seroprev_Geneva_vacc <- c(20.8, 31.4, 41.0, 63.6, 61.1, 64.4, 84.7, 89.2, 93.1)
seroprevL_Geneva_vacc <- c(15.5, 26.7, 35.6, 57.3, 56.1, 60.2, 81.1, 84.2, 89.6)
seroprevH_Geneva_vacc <- c(26.7, 36.4, 46.4, 69.8, 65.9, 68.6, 88.1, 93.4, 96.0)

# Adjust seroprevalence to match outcomes
seroprev_Geneva <- c(mean(seroprev_Geneva[1:2]), seroprev_Geneva[3:8],
  mean(seroprev_Geneva[8:9]), seroprev_Geneva[9])
seroprevL_Geneva <- c(mean(seroprevL_Geneva[1:2]), seroprevL_Geneva[3:8],
  mean(seroprevL_Geneva[8:9]), seroprevL_Geneva[9])
seroprevH_Geneva <- c(mean(seroprevH_Geneva[1:2]), seroprevH_Geneva[3:8],
  mean(seroprevH_Geneva[8:9]), seroprevH_Geneva[9])

### OUTCOMES
# Hospital and death data from:
# https://www.covid19.admin.ch/en/weekly-report/hosp?geoView=table
hospDataGE <- read.csv("../data/downloaded_datasets/switzerland/COVID19Hosp_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="GE" & datum_dboardformated=="2021-24" &
                altersklasse_covid19!="Unbekannt")
deathDataGE <- read.csv("../data/downloaded_datasets/switzerland/nuevos_files/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="GE" & datum_dboardformated=="2021-24" &
                altersklasse_covid19!="Unbekannt")
hospitalized_Geneva <- hospDataGE$sumTotal
deaths_Geneva <- deathDataGE$sumTotal
population_Geneva <- deathDataGE$pop

### VACCINATION
# https://www.corona-immunitas.ch/media/qrnlmrqp/report_on_status_of_vaccination_commented_122021_2021-12-10.pdf
# https://www.covid19.admin.ch/en/vaccination/doses
vaccDataGE <- read.csv("../data/downloaded_datasets/switzerland/COVID19VaccPersons_AKL10_w_v2.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="GE" & date==202124 &
                altersklasse_covid19!="Unbekannt" &
                age_group_type=="age_group_AKL10")
oneVaccGE <- vaccDataGE %>%
  dplyr::filter(., type=="COVID19AtLeastOneDosePersons")
twoVaccGE <- vaccDataGE %>%
  dplyr::filter(., type=="COVID19FullyVaccPersons")


### FULL TABLE
Geneva_w3 <- data.frame(Age=age_Geneva_outcomes,
                  Population=population_Geneva,
                  Prevalence=seroprev_Geneva,
                  PrevalenceL=seroprevL_Geneva,
                  PrevalenceH=seroprevH_Geneva,
                  Prevalence_unvacc=seroprev_Geneva,
                  PrevalenceL_unvacc=seroprevL_Geneva,
                  PrevalenceH_unvacc=seroprevH_Geneva,
                  Severe=hospitalized_Geneva,
                  Critical=NA,
                  Deaths=deaths_Geneva,
                  Deaths_unvacc=NA,
                  onePlusVaccine=oneVaccGE$per100PersonsTotal,
                  twoPlusVaccine=twoVaccGE$per100PersonsTotal,
                  Type="Seroprevalence",
                  Location="Geneva",
                  EndPointOutcome="2021-06-20",
                  MidPointCases="2021-06-19",
                  Period=3)


######################
# Neuchatel, Switzerland
######################

###########
### TIME 1
###########

### SEROPREVALENCE
# https://www.ne.ch/autorites/DFS/SCSP/medecin-cantonal/maladies-vaccinations/coronaimmunitas/Pages/R%c3%a9sultats-de-l'%c3%a9tude.aspx
age_Neuchatel_seroprev <- c("20-64", "65+")
# Phase 2, October 2020
seroprev_Neuchatel <- c(4.99, 4.27)
seroprevL_Neuchatel <- c(2.29, 1.06)
seroprevH_Neuchatel <- c(8.73, 8.66)

### OUTCOMES
# Hospital and death data from:
# https://www.covid19.admin.ch/en/weekly-report/hosp?geoView=table
hospDataNE <- read.csv("../data/downloaded_datasets/switzerland/COVID19Hosp_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="NE" & datum_dboardformated=="2020-42" &
                altersklasse_covid19!="Unbekannt")
deathDataNE <- read.csv("../data/downloaded_datasets/switzerland/nuevos_files/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="NE" & datum_dboardformated=="2020-42" &
                altersklasse_covid19!="Unbekannt")
hospitalized_Neuchatel <- hospDataNE$sumTotal
deaths_Neuchatel <- deathDataNE$sumTotal
population_Neuchatel <- deathDataNE$pop

# Change outcome bins to match seroprevalence
relIFR <- levinIFR(seq(60,69))/sum(levinIFR(seq(60,69)))
relIFR <- c(sum(relIFR[1:5]), sum(relIFR[6:10]))
relISR <- dhISR(seq(60,69))/sum(dhISR(seq(60,69)))
relISR <- c(sum(relISR[1:5]), sum(relISR[6:10]))

hospitalized_Neuchatel <- round(c(sum(hospitalized_Neuchatel[3:6])+
                            hospitalized_Neuchatel[7]*relISR[1],
                          hospitalized_Neuchatel[7]*relISR[2]+
                            sum(hospitalized_Neuchatel[8:9])))

deaths_Neuchatel <- round(c(sum(deaths_Neuchatel[3:6])+
                            deaths_Neuchatel[7]*relIFR[1],
                          deaths_Neuchatel[7]*relIFR[2]+
                            sum(deaths_Neuchatel[8:9])))

population_Neuchatel <- round(c(sum(population_Neuchatel[3:6])+
                            population_Neuchatel[7]*0.5,
                          population_Neuchatel[7]*0.5+
                            sum(population_Neuchatel[8:9])))

### FULL TABLE
Neuchatel_w1 <- data.frame(Age=age_Neuchatel_seroprev,
                  Population=population_Neuchatel,
                  Prevalence=seroprev_Neuchatel,
                  PrevalenceL=seroprevL_Neuchatel,
                  PrevalenceH=seroprevH_Neuchatel,
                  Prevalence_unvacc=seroprev_Neuchatel,
                  PrevalenceL_unvacc=seroprevL_Neuchatel,
                  PrevalenceH_unvacc=seroprevH_Neuchatel,
                  Severe=hospitalized_Neuchatel,
                  Critical=NA,
                  Deaths=deaths_Neuchatel,
                  Deaths_unvacc=NA,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Neuchatel",
                  EndPointOutcome="2020-10-18",
                  MidPointCases="2020-10-15",
                  Period=1)


###########
### TIME 2
###########

### SEROPREVALENCE
# https://www.ne.ch/autorites/DFS/SCSP/medecin-cantonal/maladies-vaccinations/coronaimmunitas/Pages/R%c3%a9sultats-de-l'%c3%a9tude.aspx
# Phase 3, March 2021. Excluding vaccinated
seroprev_Neuchatel_unvacc <- c(22.1, 18.8)
seroprevL_Neuchatel_unvacc <- c(17.1, 13.9)
seroprevH_Neuchatel_unvacc <- c(27.3, 24.1)

seroprev_Neuchatel <- c(22.4, 20.9)
seroprevL_Neuchatel <- c(17.5, 15.9)
seroprevH_Neuchatel <- c(27.7, 26.2)

### OUTCOMES
# Hospital and death data from:
# https://www.covid19.admin.ch/en/weekly-report/hosp?geoView=table
hospDataNE <- read.csv("../data/downloaded_datasets/switzerland/COVID19Hosp_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="NE" & datum_dboardformated=="2021-11" &
                altersklasse_covid19!="Unbekannt")
deathDataNE <- read.csv("../data/downloaded_datasets/switzerland/nuevos_files/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="NE" & datum_dboardformated=="2021-11" &
                altersklasse_covid19!="Unbekannt")
hospitalized_Neuchatel <- hospDataNE$sumTotal
deaths_Neuchatel <- deathDataNE$sumTotal
population_Neuchatel <- deathDataNE$pop

# Change outcome bins to match seroprevalence
relIFR <- levinIFR(seq(60,69))/sum(levinIFR(seq(60,69)))
relIFR <- c(sum(relIFR[1:5]), sum(relIFR[6:10]))
relISR <- dhISR(seq(60,69))/sum(dhISR(seq(60,69)))
relISR <- c(sum(relISR[1:5]), sum(relISR[6:10]))

hospitalized_Neuchatel <- round(c(sum(hospitalized_Neuchatel[3:6])+
                            hospitalized_Neuchatel[7]*relISR[1],
                          hospitalized_Neuchatel[7]*relISR[2]+
                            sum(hospitalized_Neuchatel[8:9])))

deaths_Neuchatel <- round(c(sum(deaths_Neuchatel[3:6])+
                            deaths_Neuchatel[7]*relIFR[1],
                          deaths_Neuchatel[7]*relIFR[2]+
                            sum(deaths_Neuchatel[8:9])))

population_Neuchatel <- round(c(sum(population_Neuchatel[3:6])+
                            population_Neuchatel[7]*0.5,
                          population_Neuchatel[7]*0.5+
                            sum(population_Neuchatel[8:9])))

### VACCINATION
# https://www.corona-immunitas.ch/media/qrnlmrqp/report_on_status_of_vaccination_commented_122021_2021-12-10.pdf
# https://www.covid19.admin.ch/en/vaccination/doses
vaccDataNE <- read.csv("../data/downloaded_datasets/switzerland/COVID19VaccPersons_AKL10_w_v2.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="NE" & date==202111 &
                altersklasse_covid19!="Unbekannt" &
                age_group_type=="age_group_AKL10")
oneVaccNE <- vaccDataNE %>%
  dplyr::filter(., type=="COVID19AtLeastOneDosePersons")
twoVaccNE <- vaccDataNE %>%
  dplyr::filter(., type=="COVID19FullyVaccPersons")
# Adjust vaccination to seroprevalence bins
oneVaccNE_age1 <- oneVaccNE[3:7,]
oneVacc1 <- sum((oneVaccNE_age1$pop/sum(oneVaccNE_age1$pop))*oneVaccNE_age1$per100PersonsTotal)
oneVaccNE_age2 <- oneVaccNE[7:9,]
oneVacc2 <- sum((oneVaccNE_age2$pop/sum(oneVaccNE_age2$pop))*oneVaccNE_age2$per100PersonsTotal)
oneVaccNE <- signif(c(oneVacc1, oneVacc2), digits=3)
twoVaccNE_age1 <- twoVaccNE[3:7,]
twoVacc1 <- sum((twoVaccNE_age1$pop/sum(twoVaccNE_age1$pop))*twoVaccNE_age1$per100PersonsTotal)
twoVaccNE_age2 <- twoVaccNE[7:9,]
twoVacc2 <- sum((twoVaccNE_age2$pop/sum(twoVaccNE_age2$pop))*twoVaccNE_age2$per100PersonsTotal)
twoVaccNE <- signif(c(twoVacc1, twoVacc2), digits=3)

### FULL TABLE
Neuchatel_w2 <- data.frame(Age=age_Neuchatel_seroprev,
                  Population=population_Neuchatel,
                  Prevalence=seroprev_Neuchatel,
                  PrevalenceL=seroprevL_Neuchatel,
                  PrevalenceH=seroprevH_Neuchatel,
                  Prevalence_unvacc=seroprev_Neuchatel_unvacc,
                  PrevalenceL_unvacc=seroprevL_Neuchatel_unvacc,
                  PrevalenceH_unvacc=seroprevH_Neuchatel_unvacc,
                  Severe=hospitalized_Neuchatel,
                  Critical=NA,
                  Deaths=deaths_Neuchatel,
                  Deaths_unvacc=NA,
                  onePlusVaccine=oneVaccNE,
                  twoPlusVaccine=twoVaccNE,
                  Type="Seroprevalence",
                  Location="Neuchatel",
                  EndPointOutcome="2021-03-18",
                  MidPointCases="2021-03-15",
                  Period=2)


###########
### TIME 3
###########

### SEROPREVALENCE
# https://www.ne.ch/autorites/DFS/SCSP/medecin-cantonal/maladies-vaccinations/coronaimmunitas/Pages/R%c3%a9sultats-de-l'%c3%a9tude.aspx
# Phase 4, August 2021. With vaccinated
seroprev_Neuchatel <- c(74.8, 93.0)
seroprevL_Neuchatel <- c(69.3, 89.2)
seroprevH_Neuchatel <- c(79.8, 96.2)

### OUTCOMES
# Hospital and death data from:
# https://www.covid19.admin.ch/en/weekly-report/hosp?geoView=table
hospDataNE <- read.csv("../data/downloaded_datasets/switzerland/COVID19Hosp_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="NE" & datum_dboardformated=="2021-30" &
                altersklasse_covid19!="Unbekannt")
deathDataNE <- read.csv("../data/downloaded_datasets/switzerland/nuevos_files/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="NE" & datum_dboardformated=="2021-30" &
                altersklasse_covid19!="Unbekannt")
hospitalized_Neuchatel <- hospDataNE$sumTotal
deaths_Neuchatel <- deathDataNE$sumTotal
population_Neuchatel <- deathDataNE$pop

# Change outcome bins to match seroprevalence
relIFR <- levinIFR(seq(60,69))/sum(levinIFR(seq(60,69)))
relIFR <- c(sum(relIFR[1:5]), sum(relIFR[6:10]))
relISR <- dhISR(seq(60,69))/sum(dhISR(seq(60,69)))
relISR <- c(sum(relISR[1:5]), sum(relISR[6:10]))

hospitalized_Neuchatel <- round(c(sum(hospitalized_Neuchatel[3:6])+
                            hospitalized_Neuchatel[7]*relISR[1],
                          hospitalized_Neuchatel[7]*relISR[2]+
                            sum(hospitalized_Neuchatel[8:9])))

deaths_Neuchatel <- round(c(sum(deaths_Neuchatel[3:6])+
                            deaths_Neuchatel[7]*relIFR[1],
                          deaths_Neuchatel[7]*relIFR[2]+
                            sum(deaths_Neuchatel[8:9])))

population_Neuchatel <- round(c(sum(population_Neuchatel[3:6])+
                            population_Neuchatel[7]*0.5,
                          population_Neuchatel[7]*0.5+
                            sum(population_Neuchatel[8:9])))

### VACCINATION
# https://www.corona-immunitas.ch/media/qrnlmrqp/report_on_status_of_vaccination_commented_122021_2021-12-10.pdf
# https://www.covid19.admin.ch/en/vaccination/doses
vaccDataNE <- read.csv("../data/downloaded_datasets/switzerland/COVID19VaccPersons_AKL10_w_v2.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="NE" & date==202130 &
                altersklasse_covid19!="Unbekannt" &
                age_group_type=="age_group_AKL10")
oneVaccNE <- vaccDataNE %>%
  dplyr::filter(., type=="COVID19AtLeastOneDosePersons")
twoVaccNE <- vaccDataNE %>%
  dplyr::filter(., type=="COVID19FullyVaccPersons")
# Adjust vaccination to seroprevalence bins
oneVaccNE_age1 <- oneVaccNE[3:7,]
oneVacc1 <- sum((oneVaccNE_age1$pop/sum(oneVaccNE_age1$pop))*oneVaccNE_age1$per100PersonsTotal)
oneVaccNE_age2 <- oneVaccNE[7:9,]
oneVacc2 <- sum((oneVaccNE_age2$pop/sum(oneVaccNE_age2$pop))*oneVaccNE_age2$per100PersonsTotal)
oneVaccNE <- signif(c(oneVacc1, oneVacc2), digits=3)
twoVaccNE_age1 <- twoVaccNE[3:7,]
twoVacc1 <- sum((twoVaccNE_age1$pop/sum(twoVaccNE_age1$pop))*twoVaccNE_age1$per100PersonsTotal)
twoVaccNE_age2 <- twoVaccNE[7:9,]
twoVacc2 <- sum((twoVaccNE_age2$pop/sum(twoVaccNE_age2$pop))*twoVaccNE_age2$per100PersonsTotal)
twoVaccNE <- signif(c(twoVacc1, twoVacc2), digits=3)


### FULL TABLE
Neuchatel_w3 <- data.frame(Age=age_Neuchatel_seroprev,
                  Population=population_Neuchatel,
                  Prevalence=seroprev_Neuchatel,
                  PrevalenceL=seroprevL_Neuchatel,
                  PrevalenceH=seroprevH_Neuchatel,
                  Prevalence_unvacc=NA,
                  PrevalenceL_unvacc=NA,
                  PrevalenceH_unvacc=NA,
                  Severe=hospitalized_Neuchatel,
                  Critical=NA,
                  Deaths=deaths_Neuchatel,
                  Deaths_unvacc=NA,
                  onePlusVaccine=oneVaccNE,
                  twoPlusVaccine=twoVaccNE,
                  Type="Seroprevalence",
                  Location="Neuchatel",
                  EndPointOutcome="2021-08-01",
                  MidPointCases="2021-07-28",
                  Period=3)


######################
# Fribourg, Switzerland
######################

# There is some info on vaccination status on reports, although
# seroprevalence is not age stratified and vaccine status segregated

########
### TIME 1
########

### SEROPREVALENCE
age_Fribourg_seroprev <- c("20-64", "65+")
tested_Fribourg <- c(227, 191)
positive_Fribourg <- c(20, 13)
rawPrevalenceCI <- binomial_confint(tested_Fribourg, positive_Fribourg)
seroprev_Fribourg <- signif(positive_Fribourg/tested_Fribourg*100, digits=3)
seroprevL_Fribourg <- signif(rawPrevalenceCI$lower*100, digits=3)
seroprevH_Fribourg <- signif(rawPrevalenceCI$upper*100, digits=3)

### OUTCOMES
# Hospital and death data from:
# https://www.covid19.admin.ch/en/weekly-report/hosp?geoView=table
hospDataFR <- read.csv("../data/downloaded_datasets/switzerland/COVID19Hosp_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="FR" & datum_dboardformated=="2020-38" &
                altersklasse_covid19!="Unbekannt")

deathDataFR <- read.csv("../data/downloaded_datasets/switzerland/nuevos_files/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="FR" & datum_dboardformated=="2020-38" &
                altersklasse_covid19!="Unbekannt")
hospitalized_Fribourg <- hospDataFR$sumTotal
deaths_Fribourg <- deathDataFR$sumTotal
population_Fribourg <- deathDataFR$pop
relIFR <- levinIFR(seq(60,69))/sum(levinIFR(seq(60,69)))
relIFR <- c(sum(relIFR[1:5]), sum(relIFR[6:10]))
relISR <- dhISR(seq(60,69))/sum(dhISR(seq(60,69)))
relISR <- c(sum(relISR[1:5]), sum(relISR[6:10]))

hospitalized_Fribourg <- round(c(sum(hospitalized_Fribourg[3:6])+
                            hospitalized_Fribourg[7]*relISR[1],
                          hospitalized_Fribourg[7]*relISR[2]+
                            sum(hospitalized_Fribourg[8:9])))
deaths_Fribourg <- round(c(sum(deaths_Fribourg[3:6])+
                            deaths_Fribourg[7]*relIFR[1],
                          deaths_Fribourg[7]*relIFR[2]+
                            sum(deaths_Fribourg[8:9])))
population_Fribourg <- round(c(sum(population_Fribourg[3:6])+
                            population_Fribourg[7]*0.5,
                          population_Fribourg[7]*0.5+
                            sum(population_Fribourg[8:9])))

Fribourg_w1 <- data.frame(Age=age_Fribourg_seroprev,
                  Population=population_Fribourg,
                  Prevalence=seroprev_Fribourg,
                  PrevalenceL=seroprevL_Fribourg,
                  PrevalenceH=seroprevH_Fribourg,
                  Prevalence_unvacc=seroprev_Fribourg,
                  PrevalenceL_unvacc=seroprevL_Fribourg,
                  PrevalenceH_unvacc=seroprevH_Fribourg,
                  Severe=hospitalized_Fribourg,
                  Critical=NA,
                  Deaths=deaths_Fribourg,
                  Deaths_unvacc=deaths_Fribourg,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Fribourg",
                  EndPointOutcome="2020-09-20",
                  MidPointCases="2020-09-20",
                  Period=1)

#########
### TIME 2
#########

### SEROPREVALENCE
age_Fribourg_seroprev <- c("20-64", "65+")
tested_Fribourg <- c(302, 147)
positive_Fribourg <- c(57, 35)
rawPrevalenceCI <- binomial_confint(tested_Fribourg, positive_Fribourg)
seroprev_Fribourg <- signif(positive_Fribourg/tested_Fribourg*100, digits=3)
seroprevL_Fribourg <- signif(rawPrevalenceCI$lower*100, digits=3)
seroprevH_Fribourg <- signif(rawPrevalenceCI$upper*100, digits=3)

### OUTCOMES
# Hospital and death data from:
# https://www.covid19.admin.ch/en/weekly-report/hosp?geoView=table
hospDataFR <- read.csv("../data/downloaded_datasets/switzerland/COVID19Hosp_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="FR" & datum_dboardformated=="2020-52" &
                altersklasse_covid19!="Unbekannt")

deathDataFR <- read.csv("../data/downloaded_datasets/switzerland/nuevos_files/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="FR" & datum_dboardformated=="2020-52" &
                altersklasse_covid19!="Unbekannt")
hospitalized_Fribourg <- hospDataFR$sumTotal
deaths_Fribourg <- deathDataFR$sumTotal
population_Fribourg <- deathDataFR$pop
relIFR <- levinIFR(seq(60,69))/sum(levinIFR(seq(60,69)))
relIFR <- c(sum(relIFR[1:5]), sum(relIFR[6:10]))
relISR <- dhISR(seq(60,69))/sum(dhISR(seq(60,69)))
relISR <- c(sum(relISR[1:5]), sum(relISR[6:10]))

hospitalized_Fribourg <- round(c(sum(hospitalized_Fribourg[3:6])+
                            hospitalized_Fribourg[7]*relISR[1],
                          hospitalized_Fribourg[7]*relISR[2]+
                            sum(hospitalized_Fribourg[8:9])))
deaths_Fribourg <- round(c(sum(deaths_Fribourg[3:6])+
                            deaths_Fribourg[7]*relIFR[1],
                          deaths_Fribourg[7]*relIFR[2]+
                            sum(deaths_Fribourg[8:9])))
population_Fribourg <- round(c(sum(population_Fribourg[3:6])+
                            population_Fribourg[7]*0.5,
                          population_Fribourg[7]*0.5+
                            sum(population_Fribourg[8:9])))


### FULL TABLE
Fribourg_w2 <- data.frame(Age=age_Fribourg_seroprev,
                  Population=population_Fribourg,
                  Prevalence=seroprev_Fribourg,
                  PrevalenceL=seroprevL_Fribourg,
                  PrevalenceH=seroprevH_Fribourg,
                  Prevalence_unvacc=seroprev_Fribourg,
                  PrevalenceL_unvacc=seroprevL_Fribourg,
                  PrevalenceH_unvacc=seroprevH_Fribourg,
                  Severe=hospitalized_Fribourg,
                  Critical=NA,
                  Deaths=deaths_Fribourg,
                  Deaths_unvacc=deaths_Fribourg,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Fribourg",
                  EndPointOutcome="2020-12-23",
                  MidPointCases="2020-12-23",
                  Period=2)


#########
### TIME 3
#########

### SEROPREVALENCE
age_Fribourg_seroprev <- c("20-64", "65+")
tested_Fribourg <- c(261, 243)
positive_Fribourg <- c(176, 218)
rawPrevalenceCI <- binomial_confint(tested_Fribourg, positive_Fribourg)
seroprev_Fribourg <- signif(positive_Fribourg/tested_Fribourg*100, digits=3)
seroprevL_Fribourg <- signif(rawPrevalenceCI$lower*100, digits=3)
seroprevH_Fribourg <- signif(rawPrevalenceCI$upper*100, digits=3)

### OUTCOMES
# Hospital and death data from:
# https://www.covid19.admin.ch/en/weekly-report/hosp?geoView=table
hospDataFR <- read.csv("../data/downloaded_datasets/switzerland/COVID19Hosp_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="FR" & datum_dboardformated=="2021-26" &
                altersklasse_covid19!="Unbekannt")
deathDataFR <- read.csv("../data/downloaded_datasets/switzerland/nuevos_files/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="FR" & datum_dboardformated=="2021-28" &
                altersklasse_covid19!="Unbekannt")
hospitalized_Fribourg <- hospDataFR$sumTotal
deaths_Fribourg <- deathDataFR$sumTotal
population_Fribourg <- deathDataFR$pop
relIFR <- levinIFR(seq(60,69))/sum(levinIFR(seq(60,69)))
relIFR <- c(sum(relIFR[1:5]), sum(relIFR[6:10]))
relISR <- dhISR(seq(60,69))/sum(dhISR(seq(60,69)))
relISR <- c(sum(relISR[1:5]), sum(relISR[6:10]))

hospitalized_Fribourg <- round(c(sum(hospitalized_Fribourg[3:6])+
                            hospitalized_Fribourg[7]*relISR[1],
                          hospitalized_Fribourg[7]*relISR[2]+
                            sum(hospitalized_Fribourg[8:9])))
deaths_Fribourg <- round(c(sum(deaths_Fribourg[3:6])+
                            deaths_Fribourg[7]*relIFR[1],
                          deaths_Fribourg[7]*relIFR[2]+
                            sum(deaths_Fribourg[8:9])))
population_Fribourg <- round(c(sum(population_Fribourg[3:6])+
                            population_Fribourg[7]*0.5,
                          population_Fribourg[7]*0.5+
                            sum(population_Fribourg[8:9])))

### VACCINATION
vaccDataFR <- read.csv("../data/downloaded_datasets/switzerland/COVID19VaccPersons_AKL10_w_v2.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="FR" & date==202128 &
                altersklasse_covid19!="Unbekannt" &
                age_group_type=="age_group_AKL10")
oneVaccFR <- vaccDataFR %>%
  dplyr::filter(., type=="COVID19AtLeastOneDosePersons")
twoVaccFR <- vaccDataFR %>%
  dplyr::filter(., type=="COVID19FullyVaccPersons")
# Adjust vaccination to seroprevalence bins
oneVaccFR_age1 <- oneVaccFR[3:7,]
oneVacc1 <- sum((oneVaccFR_age1$pop/sum(oneVaccFR_age1$pop))*oneVaccFR_age1$per100PersonsTotal)
oneVaccFR_age2 <- oneVaccFR[7:9,]
oneVacc2 <- sum((oneVaccFR_age2$pop/sum(oneVaccFR_age2$pop))*oneVaccFR_age2$per100PersonsTotal)
oneVaccFR <- signif(c(oneVacc1, oneVacc2), digits=3)
twoVaccFR_age1 <- twoVaccFR[3:7,]
twoVacc1 <- sum((twoVaccFR_age1$pop/sum(twoVaccFR_age1$pop))*twoVaccFR_age1$per100PersonsTotal)
twoVaccFR_age2 <- twoVaccFR[7:9,]
twoVacc2 <- sum((twoVaccFR_age2$pop/sum(twoVaccFR_age2$pop))*twoVaccFR_age2$per100PersonsTotal)
twoVaccFR <- signif(c(twoVacc1, twoVacc2), digits=3)


### FULL TABLE
Fribourg_w3 <- data.frame(Age=age_Fribourg_seroprev,
                  Population=population_Fribourg,
                  Prevalence=seroprev_Fribourg,
                  PrevalenceL=seroprevL_Fribourg,
                  PrevalenceH=seroprevH_Fribourg,
                  Prevalence_unvacc=NA,
                  PrevalenceL_unvacc=NA,
                  PrevalenceH_unvacc=NA,
                  Severe=hospitalized_Fribourg,
                  Critical=NA,
                  Deaths=deaths_Fribourg,
                  Deaths_unvacc=NA,
                  onePlusVaccine=oneVaccFR,
                  twoPlusVaccine=twoVaccFR,
                  Type="Seroprevalence",
                  Location="Fribourg",
                  EndPointOutcome="2021-07-18",
                  MidPointCases="2021-07-18",
                  Period=3)


######################
# Ticino, Switzerland
######################

#### Ticino test captures only infection antibodies

########
### TIME 1
########

### SEROPREVALENCE
age_Ticino_seroprev <- c("5-9", "10-19", "20-49", "50-64", "65+")
seroprev_Ticino <- c(3.7, 9.2, 8.1, 10.8, 9.7)
seroprevL_Ticino <- c(0.0, 4.2, 5.3, 6.8, 5.8)
seroprevH_Ticino <- c(15.0, 17.6, 11.7, 15.9, 14.8)

### OUTCOMES
# Hospital and death data from:
# https://www.covid19.admin.ch/en/weekly-report/hosp?geoView=table
hospDataTI <- read.csv("../data/downloaded_datasets/switzerland/COVID19Hosp_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="TI" & datum_dboardformated=="2020-23" &
                altersklasse_covid19!="Unbekannt")

deathDataTI <- read.csv("../data/downloaded_datasets/switzerland/nuevos_files/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="TI" & datum_dboardformated=="2020-23" &
                altersklasse_covid19!="Unbekannt")
hospitalized_Ticino <- hospDataTI$sumTotal
deaths_Ticino <- deathDataTI$sumTotal
population_Ticino <- deathDataTI$pop

relIFR <- levinIFR(seq(60,69))/sum(levinIFR(seq(60,69)))
relIFR <- c(sum(relIFR[1:5]), sum(relIFR[6:10]))
relISR <- dhISR(seq(60,69))/sum(dhISR(seq(60,69)))
relISR <- c(sum(relISR[1:5]), sum(relISR[6:10]))

hospitalized_Ticino <- c(round(hospitalized_Ticino[1]/2),
  hospitalized_Ticino[2],
  sum(hospitalized_Ticino[3:5]),
  hospitalized_Ticino[6]+round(hospitalized_Ticino[7]*relISR[1]),
  round(hospitalized_Ticino[7]*relISR[2])+sum(hospitalized_Ticino[8:9]))

deaths_Ticino <- c(round(deaths_Ticino[1]/2),
  deaths_Ticino[2],
  sum(deaths_Ticino[3:5]),
  deaths_Ticino[6]+round(deaths_Ticino[7]*relIFR[1]),
  round(deaths_Ticino[7]*relIFR[2])+sum(deaths_Ticino[8:9]))

population_Ticino <- c(round(population_Ticino[1]/2),
  population_Ticino[2],
  sum(population_Ticino[3:5]),
  population_Ticino[6]+round(population_Ticino[7]*0.5),
  round(population_Ticino[7]*0.5)+sum(population_Ticino[8:9]))


Ticino_w1 <- data.frame(Age=age_Ticino_seroprev,
                  Population=population_Ticino,
                  Prevalence=seroprev_Ticino,
                  PrevalenceL=seroprevL_Ticino,
                  PrevalenceH=seroprevH_Ticino,
                  Prevalence_unvacc=seroprev_Ticino,
                  PrevalenceL_unvacc=seroprevL_Ticino,
                  PrevalenceH_unvacc=seroprevH_Ticino,
                  Severe=hospitalized_Ticino,
                  Critical=NA,
                  Deaths=deaths_Ticino,
                  Deaths_unvacc=deaths_Ticino,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Ticino",
                  EndPointOutcome="2020-06-03",
                  MidPointCases="2020-06-07",
                  Period=1)


########
### TIME 2
########

### SEROPREVALENCE
age_Ticino_seroprev <- c("5-9", "10-19", "20-49", "50-64", "65+")
seroprev_Ticino <- c(3.9, 5.8, 8.4, 9.8, 9.0)
seroprevL_Ticino <- c(0.0, 1.9, 5.5, 5.8, 5.2)
seroprevH_Ticino <- c(15.8, 13.1, 12.1, 14.8, 13.9)

### OUTCOMES
# Hospital and death data from:
# https://www.covid19.admin.ch/en/weekly-report/hosp?geoView=table
hospDataTI <- read.csv("../data/downloaded_datasets/switzerland/COVID19Hosp_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="TI" & datum_dboardformated=="2020-35" &
                altersklasse_covid19!="Unbekannt")

deathDataTI <- read.csv("../data/downloaded_datasets/switzerland/nuevos_files/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="TI" & datum_dboardformated=="2020-35" &
                altersklasse_covid19!="Unbekannt")
hospitalized_Ticino <- hospDataTI$sumTotal
deaths_Ticino <- deathDataTI$sumTotal
population_Ticino <- deathDataTI$pop

relIFR <- levinIFR(seq(60,69))/sum(levinIFR(seq(60,69)))
relIFR <- c(sum(relIFR[1:5]), sum(relIFR[6:10]))
relISR <- dhISR(seq(60,69))/sum(dhISR(seq(60,69)))
relISR <- c(sum(relISR[1:5]), sum(relISR[6:10]))

hospitalized_Ticino <- c(round(hospitalized_Ticino[1]/2),
  hospitalized_Ticino[2],
  sum(hospitalized_Ticino[3:5]),
  hospitalized_Ticino[6]+round(hospitalized_Ticino[7]*relISR[1]),
  round(hospitalized_Ticino[7]*relISR[2])+sum(hospitalized_Ticino[8:9]))

deaths_Ticino <- c(round(deaths_Ticino[1]/2),
  deaths_Ticino[2],
  sum(deaths_Ticino[3:5]),
  deaths_Ticino[6]+round(deaths_Ticino[7]*relIFR[1]),
  round(deaths_Ticino[7]*relIFR[2])+sum(deaths_Ticino[8:9]))

population_Ticino <- c(round(population_Ticino[1]/2),
  population_Ticino[2],
  sum(population_Ticino[3:5]),
  population_Ticino[6]+round(population_Ticino[7]*0.5),
  round(population_Ticino[7]*0.5)+sum(population_Ticino[8:9]))


Ticino_w2 <- data.frame(Age=age_Ticino_seroprev,
                  Population=population_Ticino,
                  Prevalence=seroprev_Ticino,
                  PrevalenceL=seroprevL_Ticino,
                  PrevalenceH=seroprevH_Ticino,
                  Prevalence_unvacc=seroprev_Ticino,
                  PrevalenceL_unvacc=seroprevL_Ticino,
                  PrevalenceH_unvacc=seroprevH_Ticino,
                  Severe=hospitalized_Ticino,
                  Critical=NA,
                  Deaths=deaths_Ticino,
                  Deaths_unvacc=deaths_Ticino,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Ticino",
                  EndPointOutcome="2020-08-30",
                  MidPointCases="2020-08-27",
                  Period=2)



########
### TIME 3
########

### SEROPREVALENCE
age_Ticino_seroprev <- c("5-9", "10-19", "20-49", "50-64", "65+")
seroprev_Ticino <- c(11.8, 12.4, 14.8, 16.0, 12.3)
seroprevL_Ticino <- c(4.1, 6.2, 11.0, 11.0, 8.0)
seroprevH_Ticino <- c(25.9, 21.4, 19.3, 22.1, 17.9)

### OUTCOMES
# Hospital and death data from:
# https://www.covid19.admin.ch/en/weekly-report/hosp?geoView=table
hospDataTI <- read.csv("../data/downloaded_datasets/switzerland/COVID19Hosp_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="TI" & datum_dboardformated=="2020-48" &
                altersklasse_covid19!="Unbekannt")

deathDataTI <- read.csv("../data/downloaded_datasets/switzerland/nuevos_files/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="TI" & datum_dboardformated=="2020-48" &
                altersklasse_covid19!="Unbekannt")
hospitalized_Ticino <- hospDataTI$sumTotal
deaths_Ticino <- deathDataTI$sumTotal
population_Ticino <- deathDataTI$pop

relIFR <- levinIFR(seq(60,69))/sum(levinIFR(seq(60,69)))
relIFR <- c(sum(relIFR[1:5]), sum(relIFR[6:10]))
relISR <- dhISR(seq(60,69))/sum(dhISR(seq(60,69)))
relISR <- c(sum(relISR[1:5]), sum(relISR[6:10]))

hospitalized_Ticino <- c(round(hospitalized_Ticino[1]/2),
  hospitalized_Ticino[2],
  sum(hospitalized_Ticino[3:5]),
  hospitalized_Ticino[6]+round(hospitalized_Ticino[7]*relISR[1]),
  round(hospitalized_Ticino[7]*relISR[2])+sum(hospitalized_Ticino[8:9]))

deaths_Ticino <- c(round(deaths_Ticino[1]/2),
  deaths_Ticino[2],
  sum(deaths_Ticino[3:5]),
  deaths_Ticino[6]+round(deaths_Ticino[7]*relIFR[1]),
  round(deaths_Ticino[7]*relIFR[2])+sum(deaths_Ticino[8:9]))

population_Ticino <- c(round(population_Ticino[1]/2),
  population_Ticino[2],
  sum(population_Ticino[3:5]),
  population_Ticino[6]+round(population_Ticino[7]*0.5),
  round(population_Ticino[7]*0.5)+sum(population_Ticino[8:9]))


Ticino_w3 <- data.frame(Age=age_Ticino_seroprev,
                  Population=population_Ticino,
                  Prevalence=seroprev_Ticino,
                  PrevalenceL=seroprevL_Ticino,
                  PrevalenceH=seroprevH_Ticino,
                  Prevalence_unvacc=seroprev_Ticino,
                  PrevalenceL_unvacc=seroprevL_Ticino,
                  PrevalenceH_unvacc=seroprevH_Ticino,
                  Severe=hospitalized_Ticino,
                  Critical=NA,
                  Deaths=deaths_Ticino,
                  Deaths_unvacc=deaths_Ticino,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Ticino",
                  EndPointOutcome="2020-11-29",
                  MidPointCases="2020-11-29",
                  Period=3)

########
### TIME 4
########

### SEROPREVALENCE
age_Ticino_seroprev <- c("5-9", "10-19", "20-49", "50-64", "65+")
seroprev_Ticino <- c(27.0, 25.8, 20.7, 22.9, 22.0)
seroprevL_Ticino <- c(13.7, 15.9, 16.0, 16.7, 15.9)
seroprevH_Ticino <- c(45.6, 37.6, 26.2, 30.4, 29.4)

### OUTCOMES
# Hospital and death data from:
# https://www.covid19.admin.ch/en/weekly-report/hosp?geoView=table
hospDataTI <- read.csv("../data/downloaded_datasets/switzerland/COVID19Hosp_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="TI" & datum_dboardformated=="2021-20" &
                altersklasse_covid19!="Unbekannt")

deathDataTI <- read.csv("../data/downloaded_datasets/switzerland/nuevos_files/COVID19Death_geoRegion_AKL10_w.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="TI" & datum_dboardformated=="2021-20" &
                altersklasse_covid19!="Unbekannt")
hospitalized_Ticino <- hospDataTI$sumTotal
deaths_Ticino <- deathDataTI$sumTotal
population_Ticino <- deathDataTI$pop

relIFR <- levinIFR(seq(60,69))/sum(levinIFR(seq(60,69)))
relIFR <- c(sum(relIFR[1:5]), sum(relIFR[6:10]))
relISR <- dhISR(seq(60,69))/sum(dhISR(seq(60,69)))
relISR <- c(sum(relISR[1:5]), sum(relISR[6:10]))

hospitalized_Ticino <- c(round(hospitalized_Ticino[1]/2),
  hospitalized_Ticino[2],
  sum(hospitalized_Ticino[3:5]),
  hospitalized_Ticino[6]+round(hospitalized_Ticino[7]*relISR[1]),
  round(hospitalized_Ticino[7]*relISR[2])+sum(hospitalized_Ticino[8:9]))

deaths_Ticino <- c(round(deaths_Ticino[1]/2),
  deaths_Ticino[2],
  sum(deaths_Ticino[3:5]),
  deaths_Ticino[6]+round(deaths_Ticino[7]*relIFR[1]),
  round(deaths_Ticino[7]*relIFR[2])+sum(deaths_Ticino[8:9]))

population_Ticino <- c(round(population_Ticino[1]/2),
  population_Ticino[2],
  sum(population_Ticino[3:5]),
  population_Ticino[6]+round(population_Ticino[7]*0.5),
  round(population_Ticino[7]*0.5)+sum(population_Ticino[8:9]))


### VACCINATION
# https://www.corona-immunitas.ch/media/qrnlmrqp/report_on_status_of_vaccination_commented_122021_2021-12-10.pdf
# https://www.covid19.admin.ch/en/vaccination/doses
vaccDataTI <- read.csv("../data/downloaded_datasets/switzerland/COVID19VaccPersons_AKL10_w_v2.csv",
                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., geoRegion=="TI" & date==202120 &
                altersklasse_covid19!="Unbekannt" &
                age_group_type=="age_group_AKL10")
oneVaccTI <- vaccDataTI %>%
  dplyr::filter(., type=="COVID19AtLeastOneDosePersons")
twoVaccTI <- vaccDataTI %>%
  dplyr::filter(., type=="COVID19FullyVaccPersons")
# Adjust vaccination to seroprevalence bins
oneVaccTI_age1 <- oneVaccTI[3:5,]
oneVacc1 <- sum((oneVaccTI_age1$pop/sum(oneVaccTI_age1$pop))*oneVaccTI_age1$per100PersonsTotal)
oneVaccTI_age2 <- oneVaccTI[6:7,]
oneVacc2 <- sum((oneVaccTI_age2$pop/sum(oneVaccTI_age2$pop))*oneVaccTI_age2$per100PersonsTotal)
oneVaccTI_age3 <- oneVaccTI[7:9,]
oneVacc3 <- sum((oneVaccTI_age3$pop/sum(oneVaccTI_age3$pop))*oneVaccTI_age3$per100PersonsTotal)
oneVaccTI <- signif(c(oneVaccTI$per100PersonsTota[1:2], oneVacc1, oneVacc2, oneVacc3), digits=3)

twoVaccTI_age1 <- twoVaccTI[3:5,]
twoVacc1 <- sum((twoVaccTI_age1$pop/sum(twoVaccTI_age1$pop))*twoVaccTI_age1$per100PersonsTotal)
twoVaccTI_age2 <- twoVaccTI[6:7,]
twoVacc2 <- sum((twoVaccTI_age2$pop/sum(twoVaccTI_age2$pop))*twoVaccTI_age2$per100PersonsTotal)
twoVaccTI_age3 <- twoVaccTI[7:9,]
twoVacc3 <- sum((twoVaccTI_age3$pop/sum(twoVaccTI_age3$pop))*twoVaccTI_age3$per100PersonsTotal)
twoVaccTI <- signif(c(twoVaccTI$per100PersonsTotal[1:2],
                    twoVacc1, twoVacc2, twoVacc3), digits=3)


### FULL TABLE
Ticino_w4 <- data.frame(Age=age_Ticino_seroprev,
                  Population=population_Ticino,
                  Prevalence=seroprev_Ticino,
                  PrevalenceL=seroprevL_Ticino,
                  PrevalenceH=seroprevH_Ticino,
                  Prevalence_unvacc=NA,
                  PrevalenceL_unvacc=NA,
                  PrevalenceH_unvacc=NA,
                  Severe=hospitalized_Ticino,
                  Critical=NA,
                  Deaths=deaths_Ticino,
                  Deaths_unvacc=NA,
                  onePlusVaccine=oneVaccTI,
                  twoPlusVaccine=twoVaccTI,
                  Type="Seroprevalence",
                  Location="Ticino",
                  EndPointOutcome="2020-11-29",
                  MidPointCases="2021-05-23",
                  Period=4)



######################
# Freiburg, Germany
######################

# Hospitalizations?
# https://zenodo.org/record/5945065#.Yfq1dPvQ-sE

#############
### TIME 1
#############

### POPULATION
# From https://citypopulation.de/en/germany/admin/baden_w%C3%BCrttemberg/08311__freiburg_im_breisgau/
#age_Freiburg_pop <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
#                      "60-69", "70-79", "80-89", "90+")
#population_Freiburg <- c(21552, 20267, 43197, 37397, 26465, 30075, 23300,
#                         15883, 10549, 2255)
# Source: Seroprevalence paper, Table 6/7
population_Freiburg <- c(119900, 164717, 97329, 30183)

### SEROPREVALENCE
# SARS-CoV-2 Seroprevalence in Germany - A Population Based Sequential Study in Five Regions
age_Freiburg_seroprev <- c("18-34", "35-59", "60-79", "80+")
seroprev_Freiburg <- c(1.6, 1.7, 1.5, 0.0)
seroprevL_Freiburg <- c(0.9, 1.1, 0.8, 0.0)
seroprevH_Freiburg <- c(2.7, 2.6, 2.7, 3.0)

### OUTCOMES
# Deaths 1: from https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv
age_Freiburg_deaths <- c("15-34", "35-59", "60-79", "80+")
dataGer<- read.csv("../data/downloaded_datasets/germany/RKI_COVID19.csv",
                    stringsAsFactors=FALSE) %>%
  tidyr::as_tibble(.)

dataFreiburg <- dataGer %>%
  dplyr::filter(., Bundesland=="Baden-Württemberg" &
                Landkreis %in% c("SK Freiburg i.Breisgau", "LK Breisgau-Hochschwarzwald")) %>%
                #Landkreis %in% c("SK Freiburg i.Breisgau", "LK Emmendingen", "LK Breisgau-Hochschwarzwald")) %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum)) %>%
  dplyr::filter(., date <= lubridate::date("2020/08/18")) %>%
  dplyr::group_by(., Altersgruppe) %>%
  dplyr::summarize(., totalDeaths=sum(AnzahlTodesfall))
deaths_Freiburg <- dataFreiburg$totalDeaths[3:6]

# Deaths 2: From the paper Table 6/7 (deaths at study start)
#deaths_Freiburg <- c(0, 4, 43, 104)

### FULL TABLE
Freiburg_w1 <- data.frame(Age=age_Freiburg_deaths,
                  Population=population_Freiburg,
                  Prevalence=seroprev_Freiburg,
                  PrevalenceL=seroprevL_Freiburg,
                  PrevalenceH=seroprevH_Freiburg,
                  Prevalence_unvacc=seroprev_Freiburg,
                  PrevalenceL_unvacc=seroprevL_Freiburg,
                  PrevalenceH_unvacc=seroprevH_Freiburg,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Freiburg,
                  Deaths_unvacc=deaths_Freiburg,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Freiburg",
                  EndPointOutcome="2020-08-18",
                  MidPointCases="2020-08-18",
                  Period=1)


#############
### TIME 2
#############

### SEROPREVALENCE
## Seroprev 2
## Pop from same study as seroprevalence (some numbers in original table are wrong,
## can be estimated from infections reported in table)
#age_Freiburg_seroprev <- c("18-34", "35-59", "60-79", "80+")
seroprev_Freiburg <- c(2.3, 2.6, 1.8, 4.4)
seroprevL_Freiburg <- c(1.3, 1.6, 0.9, 1.1)
seroprevH_Freiburg <- c(3.4, 4.1, 3.6, 15.8)

### OUTCOMES
# Deaths 2: from https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv
age_Freiburg_deaths <- c("15-34", "35-59", "60-79", "80+")
dataFreiburg <- dataGer %>%
  dplyr::filter(., Bundesland=="Baden-Württemberg" &
                Landkreis %in% c("SK Freiburg i.Breisgau", "LK Breisgau-Hochschwarzwald")) %>%
                #Landkreis %in% c("SK Freiburg i.Breisgau", "LK Emmendingen", "LK Breisgau-Hochschwarzwald")) %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum)) %>%
  dplyr::filter(., date <= lubridate::date("2020/12/01")) %>%
  dplyr::group_by(., Altersgruppe) %>%
  dplyr::summarize(., totalDeaths=sum(AnzahlTodesfall))
deaths_Freiburg <- dataFreiburg$totalDeaths[3:6]

### FULL TABLE
Freiburg_w2 <- data.frame(Age=age_Freiburg_deaths,
                  Population=population_Freiburg,
                  Prevalence=seroprev_Freiburg,
                  PrevalenceL=seroprevL_Freiburg,
                  PrevalenceH=seroprevH_Freiburg,
                  Prevalence_unvacc=seroprev_Freiburg,
                  PrevalenceL_unvacc=seroprevL_Freiburg,
                  PrevalenceH_unvacc=seroprevH_Freiburg,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Freiburg,
                  Deaths_unvacc=deaths_Freiburg,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Freiburg",
                  EndPointOutcome="2020-12-01",
                  MidPointCases="2020-12-01",
                  Period=2)


######################
# Reutlingen, Germany
######################

##########
### TIME 1
##########

### POPULATION
# Pop from same study as seroprevalence 
population_Reutlingen <- c(57802, 100381, 59367, 19604)

### SEROPREVALENCE
# Seroprevalence: SARS-CoV-2 seroprevalence in Germany -
# a population based sequential study in five regions
# Some numbers are wrong in Table 6. Here we recalculated them from reported
# infected in table
age_Reutlingen_seroprev <- c("18-34", "35-59", "60-79", "80+")
seroprev_Reutlingen <- c(2.6, 2.6, 2.3, 1.1)
seroprevL_Reutlingen <- c(1.5, 1.7, 1.4, 0.2)
seroprevH_Reutlingen <- c(4.6, 3.7, 3.8, 7.5)

### OUTCOMES
# Deaths 1: from seroprev study
#age_Reutlingen_deaths <- c("18-34", "35-59", "60-79", "80+")
#deaths_Reutlingen <- c(0, 2, 17, 62)

# Deaths 2: from https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv
dataReutlingen <- dataGer %>%
  dplyr::filter(., Bundesland=="Baden-Württemberg" &
                Landkreis=="LK Reutlingen") %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum)) %>%
  dplyr::filter(., date <= lubridate::date("2020/07/14")) %>%
  dplyr::group_by(., Altersgruppe) %>%
  dplyr::summarize(., totalDeaths=sum(AnzahlTodesfall))
deaths_Reutlingen <- dataReutlingen$totalDeaths[3:6]

### FULL TABLE
Reutlingen_w1 <- data.frame(Age=age_Reutlingen_seroprev,
                  Population=population_Reutlingen,
                  Prevalence=seroprev_Reutlingen,
                  PrevalenceL=seroprevL_Reutlingen,
                  PrevalenceH=seroprevH_Reutlingen,
                  Prevalence_unvacc=seroprev_Reutlingen,
                  PrevalenceL_unvacc=seroprevL_Reutlingen,
                  PrevalenceH_unvacc=seroprevH_Reutlingen,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Reutlingen,
                  Deaths_unvacc=deaths_Reutlingen,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Reutlingen",
                  EndPointOutcome="2020-07-14",
                  MidPointCases="2020-07-14",
                  Period=1)

##########
### TIME 2
##########

### SEROPREVALENCE
# Seroprevalence: SARS-CoV-2 seroprevalence in Germany -
# a population based sequential study in five regions
# Pop from same study as seroprevalence (some numbers in original table are wrong,
# can be estimated from infections reported in table)
age_Reutlingen_seroprev <- c("18-34", "35-59", "60-79", "80+")
seroprev_Reutlingen <- c(2.0, 3.4, 1.8, 5.7)
seroprevL_Reutlingen <- c(1.0, 2.4, 1.1, 1.9)
seroprevH_Reutlingen <- c(4.0, 4.6, 3.1, 16.1)

### OUTCOMES
# Deaths 1: from seroprev study
#age_Reutlingen_deaths <- c("18-34", "35-59", "60-79", "80+")
#deaths_Reutlingen <- c(0, 2, 17, 62)
# Deaths 2: from https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv
dataReutlingen <- dataGer %>%
  dplyr::filter(., Bundesland=="Baden-Württemberg" &
                Landkreis=="LK Reutlingen") %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum)) %>%
  dplyr::filter(., date <= lubridate::date("2020/10/27")) %>%
  dplyr::group_by(., Altersgruppe) %>%
  dplyr::summarize(., totalDeaths=sum(AnzahlTodesfall))
deaths_Reutlingen <- dataReutlingen$totalDeaths[3:6]

### FULL TABLE
Reutlingen_w2 <- data.frame(Age=age_Reutlingen_seroprev,
                  Population=population_Reutlingen,
                  Prevalence=seroprev_Reutlingen,
                  PrevalenceL=seroprevL_Reutlingen,
                  PrevalenceH=seroprevH_Reutlingen,
                  Prevalence_unvacc=seroprev_Reutlingen,
                  PrevalenceL_unvacc=seroprevL_Reutlingen,
                  PrevalenceH_unvacc=seroprevH_Reutlingen,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Reutlingen,
                  Deaths_unvacc=deaths_Reutlingen,
                  Type="Seroprevalence",
                  Location="Reutlingen",
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  EndPointOutcome="2020-10-27",
                  MidPointCases="2020-10-27",
                  Period=2)


######################
# Osnabruck, Germany
######################

############
### TIME 1
############

### POPULATION
age_Osnabruck_seroprev <- c("18-34", "35-59", "60-79", "80+")
population_Osnabruck <- c(116113, 177847, 105957, 34650)

### SEROPREVALENCE
# Pop from same study as seroprevalence (some numbers in original table are wrong,
# can be estimated from infections reported in table)
seroprev_Osnabruck <- c(1.5, 1.4, 1.5, 0)
seroprevL_Osnabruck <- c(1.0, 0.7, 0.8, 0)
seroprevH_Osnabruck <- c(1.9, 3.0, 2.6, 0)

### OUTCOMES
# Deaths 1: from seroprev study
#age_Osnabruck_deaths <- c("18-34", "35-59", "60-79", "80+")
#deaths_Osnabruck <- c(0, 8, 33, 40)
# Deaths 2: from https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv
dataOsnabruck <- dataGer %>%
  dplyr::filter(., Bundesland=="Niedersachsen" &
                Landkreis %in% c("LK Osnabrück", "SK Osnabrück")) %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum)) %>%
  dplyr::filter(., date <= lubridate::date("2020/10/15")) %>%
  dplyr::group_by(., Altersgruppe) %>%
  dplyr::summarize(., totalDeaths=sum(AnzahlTodesfall))
deaths_Osnabruck <- dataOsnabruck$totalDeaths[3:6]

### FULL TABLE
Osnabruck_w1 <- data.frame(Age=age_Osnabruck_seroprev,
                  Population=population_Osnabruck,
                  Prevalence=seroprev_Osnabruck,
                  PrevalenceL=seroprevL_Osnabruck,
                  PrevalenceH=seroprevH_Osnabruck,
                  Prevalence_unvacc=seroprev_Osnabruck,
                  PrevalenceL_unvacc=seroprevL_Osnabruck,
                  PrevalenceH_unvacc=seroprevH_Osnabruck,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Osnabruck,
                  Deaths_unvacc=deaths_Osnabruck,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Osnabruck",
                  EndPointOutcome="2020-10-30",
                  MidPointCases="2020-10-30",
                  Period=1)



############
### TIME 2
############

### POPULATION
age_Osnabruck_seroprev <- c("18-34", "35-59", "59+")
population_Osnabruck <- c(116113, 177847, 140607)

### SEROPREVALENCE
# Pop from same study as seroprevalence (some numbers in original table are wrong,
# can be estimated from infections reported in table)
seroprev_Osnabruck_unvacc <- c(4.2, 4.9, 3.0)
seroprevL_Osnabruck_unvacc <- c(2.4, 3.6, 1.8)
seroprevH_Osnabruck_unvacc <- c(7.2, 6.6, 5.2)

### OUTCOMES
# Deaths 1: from seroprev study
#age_Osnabruck_deaths <- c("18-34", "35-59", "60-79", "80+")
#deaths_Osnabruck <- c(0, 8, 33, 40)
# Deaths 2: from https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv
dataOsnabruck <- dataGer %>%
  dplyr::filter(., Bundesland=="Niedersachsen" &
                Landkreis %in% c("LK Osnabrück", "SK Osnabrück")) %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum)) %>%
  dplyr::filter(., date <= lubridate::date("2021/03/16")) %>%
  dplyr::group_by(., Altersgruppe) %>%
  dplyr::summarize(., totalDeaths=sum(AnzahlTodesfall))
deaths_Osnabruck <- c(dataOsnabruck$totalDeaths[3:4], sum(dataOsnabruck$totalDeaths[5:6]))

### VACCINE
# https://raw.githubusercontent.com/robert-koch-institut/COVID-19-Impfungen_in_Deutschland/master/Archiv/2021-07-21_Deutschland_Landkreise_COVID-19-Impfungen.csv
# https://github.com/robert-koch-institut/COVID-19-Impfungen_in_Deutschland
vaccineOsnabruck <- read.csv("../data/downloaded_datasets/germany/2021-07-21_Deutschland_Landkreise_COVID-19-Impfungen.csv",
                             stringsAsFactors=FALSE) %>%
  tidyr::as_tibble(.) %>%
  dplyr::mutate(., date=lubridate::date(Impfdatum)) %>%
  dplyr::filter(., LandkreisId_Impfort=="03459" & Impfdatum<="2021-03-16") %>%
  dplyr::group_by(., Altersgruppe, Impfschutz) %>%
  dplyr::summarize(., totalVacc=sum(Anzahl))
oneDoseOsnabruck <- vaccineOsnabruck[vaccineOsnabruck$Impfschutz==1,]$totalVacc
oneDoseOsnabruck <- c(round(oneDoseOsnabruck[2]/3), round(oneDoseOsnabruck[2]*2/3),
                      oneDoseOsnabruck[3])
twoDoseOsnabruck <- vaccineOsnabruck[vaccineOsnabruck$Impfschutz==2,]$totalVacc
twoDoseOsnabruck <- c(round(twoDoseOsnabruck[2]/3), round(twoDoseOsnabruck[2]*2/3),
                      twoDoseOsnabruck[3])

### FULL TABLE
Osnabruck_w2 <- data.frame(Age=age_Osnabruck_seroprev,
                  Population=population_Osnabruck,
                  Prevalence=NA,
                  PrevalenceL=NA,
                  PrevalenceH=NA,
                  Prevalence_unvacc=seroprev_Osnabruck_unvacc,
                  PrevalenceL_unvacc=seroprevL_Osnabruck_unvacc,
                  PrevalenceH_unvacc=seroprevH_Osnabruck_unvacc,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Osnabruck,
                  Deaths_unvacc=NA,
                  onePlusVaccine=signif(oneDoseOsnabruck/population_Osnabruck*100, digits=3),
                  twoPlusVaccine=signif(twoDoseOsnabruck/population_Osnabruck*100, digits=3),
                  Type="Seroprevalence",
                  Location="Osnabruck",
                  EndPointOutcome="2021-03-16",
                  MidPointCases="2021-03-16",
                  Period=2)


######################
# Aachen, Germany
######################

#############
### TIME 1
#############

### POPULATION
age_Aachen_seroprev <- c("18-34", "35-59", "60-79", "80+")
population_Aachen <- c(142227, 178758, 114037, 35763)

### SEROPREVALENCE
# Pop from same study as seroprevalence (some numbers in original table are wrong,
# can be estimated from infections reported in table)
age_Aachen_seroprev <- c("18-34", "35-59", "60-79", "80+")
seroprev_Aachen <- c(2.2, 1.5, 2.5, 5.4)
seroprevL_Aachen <- c(1.2, 0.9, 1.5, 2.2)
seroprevH_Aachen <- c(4.1, 2.6, 4.2, 12.3)

### OUTCOMES
# Deaths 1: from seroprev study
#age_Aachen_deaths <- c("18-34", "35-59", "60-79", "80+")
#deaths_Aachen <- c(0, 8, 30, 65)
# Deaths 2: from https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv
dataAachen <- dataGer %>%
  dplyr::filter(., Bundesland=="Nordrhein-Westfalen" &
                Landkreis=="StädteRegion Aachen") %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum)) %>%
  dplyr::filter(., date <= lubridate::date("2020/09/24")) %>%
  dplyr::group_by(., Altersgruppe) %>%
  dplyr::summarize(., totalDeaths=sum(AnzahlTodesfall))
deaths_Aachen <- dataAachen$totalDeaths[3:6]

### FULL TABLE
Aachen_w1 <- data.frame(Age=age_Aachen_seroprev,
                  Population=population_Aachen,
                  Prevalence=seroprev_Aachen,
                  PrevalenceL=seroprevL_Aachen,
                  PrevalenceH=seroprevH_Aachen,
                  Prevalence_unvacc=seroprev_Aachen,
                  PrevalenceL_unvacc=seroprevL_Aachen,
                  PrevalenceH_unvacc=seroprevH_Aachen,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Aachen,
                  Deaths_unvacc=deaths_Aachen,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Aachen",
                  EndPointOutcome="2020-09-24",
                  MidPointCases="2020-09-24",
                  Period=1)



#############
### TIME 2
#############

### SEE IF THERE IS VACCINE SEGREGATED DEATH DATA

### SEROPREVALENCE
# Pop from same study as seroprevalence (some numbers in original table are wrong,
# can be estimated from infections reported in table)
age_Aachen_seroprev <- c("18-34", "35-59", "60-79", "80+")
seroprev_Aachen <- c(6.0, 7.6, 6.6, 8.8)
seroprevL_Aachen <- c(4.2, 5.9, 4.8, 4.2)
seroprevH_Aachen <- c(8.6, 9.6, 8.9, 17.6)

# Without vaccinated individuals
seroprev_Aachen_unvacc <- c(3.9, 5.5, 5.9, 9.4)
seroprevL_Aachen_unvacc <- c(2.4, 4.1, 4.2, 4.4)
seroprevH_Aachen_unvacc <- c(6.2, 7.4, 8.2, 18.7)

### OUTCOMES
# Deaths 1: from seroprev study
age_Aachen_deaths <- c("18-34", "35-59", "60-79", "80+")
deaths_Aachen <- c(0, 15, 109, 259)

# Deaths 2: from https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv
dataAachen <- dataGer %>%
  dplyr::filter(., Bundesland=="Nordrhein-Westfalen" &
                Landkreis=="StädteRegion Aachen") %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum)) %>%
  dplyr::filter(., date <= lubridate::date("2021/02/15")) %>%
  dplyr::group_by(., Altersgruppe) %>%
  dplyr::summarize(., totalDeaths=sum(AnzahlTodesfall))
deaths_Aachen <- dataAachen$totalDeaths[3:6]

### VACCINE
# https://raw.githubusercontent.com/robert-koch-institut/COVID-19-Impfungen_in_Deutschland/master/Archiv/2021-07-21_Deutschland_Landkreise_COVID-19-Impfungen.csv
# https://github.com/robert-koch-institut/COVID-19-Impfungen_in_Deutschland
vaccineAachen <- read.csv("../data/downloaded_datasets/germany/2021-07-21_Deutschland_Landkreise_COVID-19-Impfungen.csv",
                             stringsAsFactors=FALSE) %>%
  tidyr::as_tibble(.) %>%
  dplyr::mutate(., date=lubridate::date(Impfdatum)) %>%
  dplyr::filter(., LandkreisId_Impfort=="05334" & Impfdatum<="2021-02-15") %>%
  dplyr::group_by(., Altersgruppe, Impfschutz) %>%
  dplyr::summarize(., totalVacc=sum(Anzahl))

# Redistribute vaccines to match seroprevalence bins
population_Aachen2 <- c(population_Aachen[1:2], sum(population_Aachen[3:4]))
oneDoseAachen <- vaccineAachen[vaccineAachen$Impfschutz==1,]$totalVacc
oneDoseAachen <- c(round(oneDoseAachen[2]/3), round(oneDoseAachen[2]*2/3),
                      oneDoseAachen[3])/population_Aachen2*100
oneDoseAachen <- signif(c(oneDoseAachen[1:2], rep(oneDoseAachen[3],2)), digits=3)
twoDoseAachen <- vaccineAachen[vaccineAachen$Impfschutz==2,]$totalVacc
twoDoseAachen <- c(round(twoDoseAachen[2]/3), round(twoDoseAachen[2]*2/3),
                      twoDoseAachen[3])/population_Aachen2*100
twoDoseAachen <- signif(c(twoDoseAachen[1:2], rep(twoDoseAachen[3],2)), digits=3)


### FULL TABLE
Aachen_w2 <- data.frame(Age=age_Aachen_seroprev,
                  Population=population_Aachen,
                  Prevalence=seroprev_Aachen,
                  PrevalenceL=seroprevL_Aachen,
                  PrevalenceH=seroprevH_Aachen,
                  Prevalence_unvacc=seroprev_Aachen_unvacc,
                  PrevalenceL_unvacc=seroprevL_Aachen_unvacc,
                  PrevalenceH_unvacc=seroprevH_Aachen_unvacc,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Aachen,
                  Deaths_unvacc=NA,
                  onePlusVaccine=oneDoseAachen,
                  twoPlusVaccine=twoDoseAachen,
                  Type="Seroprevalence",
                  Location="Aachen",
                  EndPointOutcome="2021-02-15",
                  MidPointCases="2021-02-15",
                  Period=2)


######################
# Magdeburg, Germany
######################

############
### TIME 1
############

### POPULATION
age_Magdeburg_seroprev <- c("18-34", "35-59", "60-79", "80+")
population_Magdeburg <- c(54425, 74127, 54468, 18576)

### SEROPREVALENCE
# Pop from same study as seroprevalence (some numbers in original table are wrong,
# can be estimated from infections reported in table)
seroprev_Magdeburg <- c(2.2, 2.4, 1.8, 5.1)
seroprevL_Magdeburg <- c(1.2, 1.6, 1.1, 2.7)
seroprevH_Magdeburg <- c(3.8, 3.5, 2.9, 9.3)

### OUTCOMES
# Deaths 1: from seroprev study
age_Magdeburg_deaths <- c("18-34", "35-59", "60-79", "80+")
deaths_Magdeburg <- c(0, 2, 3, 10)
# Deaths 2: from https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv
dataMagdeburg <- dataGer %>%
  dplyr::filter(., Bundesland=="Sachsen-Anhalt" &
                Landkreis=="SK Magdeburg") %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum)) %>%
  dplyr::filter(., date <= lubridate::date("2020/12/01")) %>%
  dplyr::group_by(., Altersgruppe) %>%
  dplyr::summarize(., totalDeaths=sum(AnzahlTodesfall))
deaths_Magdeburg <- dataMagdeburg$totalDeaths[3:6]

### FULL TABLE
Magdeburg_w1 <- data.frame(Age=age_Magdeburg_seroprev,
                  Population=population_Magdeburg,
                  Prevalence=seroprev_Magdeburg,
                  PrevalenceL=seroprevL_Magdeburg,
                  PrevalenceH=seroprevH_Magdeburg,
                  Prevalence_unvacc=seroprev_Magdeburg,
                  PrevalenceL_unvacc=seroprevL_Magdeburg,
                  PrevalenceH_unvacc=seroprevH_Magdeburg,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Magdeburg,
                  Deaths_unvacc=deaths_Magdeburg,
                  Type="Seroprevalence",
                  Location="Magdeburg",
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  EndPointOutcome="2020-12-01",
                  MidPointCases="2020-12-01",
                  Period=1)


############
### TIME 2
############

### POPULATION
age_Magdeburg_seroprev <- c("18-34", "35-59", "60+")
population_Magdeburg <- c(54425, 74127, 73004)

### SEROPREVALENCE
# Pop from same study as seroprevalence (some numbers in original table are wrong,
# can be estimated from infections reported in table)
seroprev_Magdeburg_unvacc <- c(8.2, 11.1, 8.7)
seroprevL_Magdeburg_unvacc <- c(6.1, 8.1, 7.1)
seroprevH_Magdeburg_unvacc <- c(10.9, 15.1, 10.6)

### OUTCOMES
# Deaths 1: from https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv
dataMagdeburg <- dataGer %>%
  dplyr::filter(., Bundesland=="Sachsen-Anhalt" &
                Landkreis=="SK Magdeburg") %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum)) %>%
  dplyr::filter(., date <= lubridate::date("2021/04/21")) %>%
  dplyr::group_by(., Altersgruppe) %>%
  dplyr::summarize(., totalDeaths=sum(AnzahlTodesfall))
deaths_Magdeburg <- c(dataMagdeburg$totalDeaths[3:4],
                      sum(dataMagdeburg$totalDeaths[5:6]))


### VACCINE
# https://raw.githubusercontent.com/robert-koch-institut/COVID-19-Impfungen_in_Deutschland/master/Archiv/2021-07-21_Deutschland_Landkreise_COVID-19-Impfungen.csv
# https://github.com/robert-koch-institut/COVID-19-Impfungen_in_Deutschland
vaccineMagdeburg <- read.csv("../data/downloaded_datasets/germany/2021-07-21_Deutschland_Landkreise_COVID-19-Impfungen.csv",
                             stringsAsFactors=FALSE) %>%
  tidyr::as_tibble(.) %>%
  dplyr::mutate(., date=lubridate::date(Impfdatum)) %>%
  dplyr::filter(., LandkreisId_Impfort=="15003" & Impfdatum<="2021-04-21") %>%
  dplyr::group_by(., Altersgruppe, Impfschutz) %>%
  dplyr::summarize(., totalVacc=sum(Anzahl))

oneDoseMagdeburg <- vaccineMagdeburg[vaccineMagdeburg$Impfschutz==1,]$totalVacc
oneDoseMagdeburg <- c(round(oneDoseMagdeburg[2]/3), round(oneDoseMagdeburg[2]*2/3),
                      oneDoseMagdeburg[3])
twoDoseMagdeburg <- vaccineMagdeburg[vaccineMagdeburg$Impfschutz==2,]$totalVacc
twoDoseMagdeburg <- c(round(twoDoseMagdeburg[2]/3), round(twoDoseMagdeburg[2]*2/3),
                      twoDoseMagdeburg[3])

### FULL TABLE
Magdeburg_w2 <- data.frame(Age=age_Magdeburg_seroprev,
                  Population=population_Magdeburg,
                  Prevalence=NA,
                  PrevalenceL=NA,
                  PrevalenceH=NA,
                  Prevalence_unvacc=seroprev_Magdeburg_unvacc,
                  PrevalenceL_unvacc=seroprevL_Magdeburg_unvacc,
                  PrevalenceH_unvacc=seroprevH_Magdeburg_unvacc,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Magdeburg,
                  Deaths_unvacc=NA,
                  Type="Seroprevalence",
                  Location="Magdeburg",
                  onePlusVaccine=signif(oneDoseMagdeburg/population_Magdeburg*100, digits=3),
                  twoPlusVaccine=signif(twoDoseMagdeburg/population_Magdeburg*100, digits=3),
                  EndPointOutcome="2021-04-21",
                  MidPointCases="2021-04-21",
                  Period=2)


######################
# Munich, Germany
######################

###########
### TIME 1
###########

### POPULATION
age_Munich_pop <- c("0-19", "20-34", "35-49", "50-64", "65-79", "80+")
population_Munich <- c(263053, 390382, 348651, 291562, 184764, 83308)

### SEROPREVALENCE
# From Prevalence and Risk Factors of Infection in the Representative
# COVID-19 Cohort Munich
age_Munich_seroprev <- c("0-19", "20-34", "35-49", "50-64", "65-79", "80+")
tested_Munich <- c(261, 1322, 1505, 1233, 504, 77)
positive_Munich <- c(4, 28, 29, 19, 12, 1)
rawPrevalenceCI <- binomial_confint(tested_Munich, positive_Munich)
seroprev_Munich <- signif(positive_Munich/tested_Munich*100, digits=2)
seroprevL_Munich <- signif(rawPrevalenceCI$lower*100, digits=2)
seroprevH_Munich <- signif(rawPrevalenceCI$upper*100, digits=2)
# adjust for test characteristics
sensitivity <- 0.855
specificity <- 0.998
factor <- 1/(sensitivity + specificity - 1)
seroprev_Munich <- signif((seroprev_Munich/100 + specificity - 1)*factor*100, digits=2)
seroprevL_Munich <- signif(seroprevL_Munich*factor, digits=2)
seroprevH_Munich <- signif(seroprevH_Munich*factor, digits=2)

### OUTCOMES
# Deaths: from https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv
dataMunich <- dataGer %>%
  dplyr::filter(., Bundesland=="Bayern" &
                Landkreis=="SK München") %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum)) %>%
  dplyr::filter(., date <= lubridate::date("2020/05/09")) %>%
  dplyr::group_by(., Altersgruppe) %>%
  dplyr::summarize(., totalDeaths=sum(AnzahlTodesfall))
muDeaths <- dataMunich$totalDeaths

relIFR <- levinIFR(seq(15,34))/sum(levinIFR(seq(15,34)))
relIFR <- c(sum(relIFR[1:5]), sum(relIFR[6:20]))
relIFR2 <- levinIFR(seq(35,59))/sum(levinIFR(seq(35,59)))
relIFR2 <- c(sum(relIFR2[1:15]), sum(relIFR2[16:25]))
relIFR3 <- levinIFR(seq(60,79))/sum(levinIFR(seq(60,79)))
relIFR3 <- c(sum(relIFR3[1:5]), sum(relIFR3[6:20]))
deaths_Munich <- round(c(sum(muDeaths[1:3]*c(1,1,relIFR[1])),
                   ceiling(muDeaths[3]*relIFR[2]),
                   floor(muDeaths[4]*relIFR2[1]),
                   ceiling(muDeaths[4]*relIFR2[2])+floor(muDeaths[5]*relIFR3[1]),
                   ceiling(muDeaths[5]*relIFR3[2]),
                   muDeaths[6]))

Munich_w1 <- data.frame(Age=age_Munich_seroprev,
                  Population=population_Munich,
                  Prevalence=seroprev_Munich,
                  PrevalenceL=seroprevL_Munich,
                  PrevalenceH=seroprevH_Munich,
                  Prevalence_unvacc=seroprev_Munich,
                  PrevalenceL_unvacc=seroprevL_Munich,
                  PrevalenceH_unvacc=seroprevH_Munich,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Munich,
                  Deaths_unvacc=deaths_Munich,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Munich",
                  EndPointOutcome="2020-05-09",
                  MidPointCases="2020-05-09",
                  Period=1)


###########
### TIME 2
###########

### SEROPREVALENCE
# Seroprevalence: From first to second wave: follow-up
# of the prospective COVID-19 cohort (KoCo19)
# in Munich (Germany) (Supplementary materials)
age_Munich_seroprev <- c("0-19", "20-34", "35-49", "50-64", "65-79", "80+")
tested_Munich <- c(212, 1040, 1271, 1166, 599, 145)
ever_positive_Munich <- c(9, 34, 44, 34, 18, 2)
sero_remission_Munich <- c(0, 0, 3, 2, 1, 0)
positive_Munich <- ever_positive_Munich - sero_remission_Munich

rawPrevalenceCI <- binomial_confint(tested_Munich, positive_Munich)
seroprev_Munich <- signif(positive_Munich/tested_Munich*100, digits=2)
seroprevL_Munich <- signif(rawPrevalenceCI$lower*100, digits=2)
seroprevH_Munich <- signif(rawPrevalenceCI$upper*100, digits=2)

### OUTCOMES
# Deaths: from https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv
dataMunich <- dataGer %>%
  dplyr::filter(., Bundesland=="Bayern" &
                Landkreis=="SK München") %>%
  dplyr::mutate(., date=lubridate::date(Meldedatum)) %>%
  dplyr::filter(., date <= lubridate::date("2020/12/17")) %>%
  dplyr::group_by(., Altersgruppe) %>%
  dplyr::summarize(., totalDeaths=sum(AnzahlTodesfall))

muDeaths <- dataMunich$totalDeaths

relIFR <- levinIFR(seq(15,34))/sum(levinIFR(seq(15,34)))
relIFR <- c(sum(relIFR[1:5]), sum(relIFR[6:20]))
relIFR2 <- levinIFR(seq(35,59))/sum(levinIFR(seq(35,59)))
relIFR2 <- c(sum(relIFR2[1:15]), sum(relIFR2[16:25]))
relIFR3 <- levinIFR(seq(60,79))/sum(levinIFR(seq(60,79)))
relIFR3 <- c(sum(relIFR3[1:5]), sum(relIFR3[6:20]))
deaths_Munich <- round(c(sum(muDeaths[1:3]*c(1,1,relIFR[1])),
                   ceiling(muDeaths[3]*relIFR[2]),
                   floor(muDeaths[4]*relIFR2[1]),
                   ceiling(muDeaths[4]*relIFR2[2])+floor(muDeaths[5]*relIFR3[1]),
                   ceiling(muDeaths[5]*relIFR3[2]),
                   muDeaths[6]))

### FULL TABLE
Munich_w2 <- data.frame(Age=age_Munich_seroprev,
                  Population=population_Munich,
                  Prevalence=seroprev_Munich,
                  PrevalenceL=seroprevL_Munich,
                  PrevalenceH=seroprevH_Munich,
                  Prevalence_unvacc=seroprev_Munich,
                  PrevalenceL_unvacc=seroprevL_Munich,
                  PrevalenceH_unvacc=seroprevH_Munich,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Munich,
                  Deaths_unvacc=deaths_Munich,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Munich",
                  EndPointOutcome="2020-12-17",
                  MidPointCases="2020-12-17",
                  Period=2)

# Interesting: https://www.covid19.statistik.uni-muenchen.de/pdfs/codag_bericht_19.pdf
# https://www.muenchen.de/rathaus/Stadtinfos/Coronavirus-Fallzahlen.html


#################
# Denmark
#################

##########
### Time 1
##########

### POPULATION
age_Denmark_seroprev <- c("18-39", "40-64", "65+")
population_Denmark <- c(1584906, 1916667, 1173913)

### SEROPREVALENCE
# Prevalence of SARS‑CoV‑2 antibodies in Denmark: nationwide,
# population‑based seroepidemiological study.
# First stage
seroprev_Denmark <- c(2.3, 0.6, 0.6)
seroprevL_Denmark <- c(1.1, 0.0, 0.0)
seroprevH_Denmark <- c(3.7, 1.5, 2.1)

### OUTCOMES
# https://covid19.ssi.dk/overvagningsdata/download-fil-med-overvaagningdata
# https://experience.arcgis.com/experience/242ec2acc014456295189631586f1d26
age_Denmark_deaths <- c("0-59", "60-69", "70-79", "80-89", "90+")
deaths_Denmark <- c(16, 53, 156, 208, 115)
age_Denmark_ICU <- c("0-39","40-49", "50-59", "60-69", "70-79", "80+")
ICU_Denmark <- c(16, 21, 54, 82, 120, 40)
age_Denmark_Hospitalized <- c("0-9", "10-19", "20-29", "30-39", "40-49",
                              "50-59", "60-69", "70-79", "80-89", "90+")
hospitalized_Denmark <- c(13, 13, 51, 87, 189, 326, 361, 563, 459, 125)

# Change outcome bins to match seroprevalence
relIFR1 <- levinIFR(seq(1,59))/sum(levinIFR(seq(1,59)))
relIFR1 <- c(sum(relIFR1[1:39]), sum(relIFR1[40:59]))
relIFR2 <- levinIFR(seq(60,69))/sum(levinIFR(seq(60,69)))
relIFR2 <- c(sum(relIFR2[1:5]), sum(relIFR2[6:10]))
deaths_Denmark <- c(round(deaths_Denmark[1]*relIFR1[1]),
                    round(deaths_Denmark[1]*relIFR1[2])+
                    round(deaths_Denmark[2]*relIFR2[1]),
                    round(deaths_Denmark[2]*relIFR2[2])+
                      sum(deaths_Denmark[3:5]))

relISR1 <- dhISR(seq(60, 69))/sum(dhISR(seq(60,69)))
relISR1 <- c(sum(relISR1[1:5]), sum(relISR1[6:10]))
hospitalized_Denmark <- c(sum(hospitalized_Denmark[3:4]),
                          sum(hospitalized_Denmark[5:6])+
                            round(hospitalized_Denmark[7]*relISR1[1]),
                          round(hospitalized_Denmark[7]*relISR1[2])+
                            sum(hospitalized_Denmark[8:10]))

### FULL TABLE
Denmark_w1 <- data.frame(Age=age_Denmark_seroprev,
                  Population=population_Denmark,
                  Prevalence=seroprev_Denmark,
                  PrevalenceL=seroprevL_Denmark,
                  PrevalenceH=seroprevH_Denmark,
                  Prevalence_unvacc=seroprev_Denmark,
                  PrevalenceL_unvacc=seroprevL_Denmark,
                  PrevalenceH_unvacc=seroprevH_Denmark,
                  Severe=hospitalized_Denmark,
                  Critical=NA,
                  Deaths=deaths_Denmark,
                  Deaths_unvacc=deaths_Denmark,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Denmark",
                  EndPointOutcome="2020-05-18",
                  MidPointCases="2020-05-18",
                  Period=1)


##########
### Time 2
##########

# All data from Espenhain et al.
# Prevalence of SARS‑CoV‑2 antibodies in Denmark: nationwide,
# population‑based seroepidemiological study
# also useful:
# https://experience.arcgis.com/experience/aa41b29149f24e20a4007a0c4e13db1d

### POPULATION
age_Denmark_seroprev <- c("12-17", "18-39", "40-64", "65+")
population_Denmark <- c(400000, 1584906, 1916667, 1173913)

### SEROPREVALENCE
# December, last phase figure
seroprev_Denmark <- c(6.5, 5.3, 3.6, 2.3)
seroprevL_Denmark <- c(3.8, 3.8, 2.5, 1.0)
seroprevH_Denmark <- c(10, 6.8, 4.7, 4.0)

### OUTCOMES
hospitalized_Denmark <- c(58, 642, 2033, 3164)
deaths_Denmark <- c(0, 3, 56, 1022)

### FULL TABLE
Denmark_w2 <- data.frame(Age=age_Denmark_seroprev,
                  Population=population_Denmark,
                  Prevalence=seroprev_Denmark,
                  PrevalenceL=seroprevL_Denmark,
                  PrevalenceH=seroprevH_Denmark,
                  Prevalence_unvacc=seroprev_Denmark,
                  PrevalenceL_unvacc=seroprevL_Denmark,
                  PrevalenceH_unvacc=seroprevH_Denmark,
                  Severe=hospitalized_Denmark,
                  Critical=NA,
                  Deaths=deaths_Denmark,
                  Deaths_unvacc=deaths_Denmark,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Denmark",
                  EndPointOutcome="2020-12-12",
                  MidPointCases="2020-12-02",
                  Period=2)




######################
# England
######################

################
### Time 1
################

### POPULATION
# population data from Levin et al
age_England_Pop <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
                     "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                     "60-64", "65-69", "70-74", "75-79", "80-84", "85+")
population_England <- c(3496750, 3135711, 3258677, 2079229, 5267401,
                        3836609, 3683915, 3732161, 4099089, 4100526,
                        3601694, 3183915, 3377162, 2674161, 2178672,
                        1777547, 1338005, 1254688)

population_seroprev_England <- c(round(population_England[4]*2/5)+population_England[5],
                                 sum(population_England[6:7]),
                                 sum(population_England[8:9]),
                                 sum(population_England[10:11]),
                                 sum(population_England[12:13]),
                                 sum(population_England[14:15]),
                                 sum(population_England[16:18]))


### SEROPREVALENCE
# Seroprevalence from:
# SARS-CoV-2 antibody prevalence in England following the first peak of the pandemic
# (supplementary)
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-21237-w/MediaObjects/41467_2021_21237_MOESM1_ESM.pdf
age_England_seroprev <- c("18-24", "25-34", "35-44", "45-54", "55-64",
                       "65-74", "75+")
seroprev_England <- c(7.9, 7.8, 6.1, 6.4, 5.9, 3.2, 3.3)
seroprevL_England <- c(7.3, 7.4, 5.7, 6.0, 5.5, 2.8, 2.9)
seroprevH_England <- c(8.5, 8.3, 6.6, 6.9, 6.4, 3.6, 3.8)

### OUTCOMES
# DEATHS
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/weeklyprovisionalfiguresondeathsregisteredinenglandandwales
englandWalesDeaths <- read.csv("../data/downloaded_datasets/england/weekly_deaths_2020_ONS.csv",
                    stringsAsFactors=FALSE) %>%
  tidyr::pivot_longer(., cols=c(2:54), names_to="week", values_to="deaths") %>%
  dplyr::mutate(., week=as.integer(gsub("X", "", week))) %>%
  dplyr::filter(., week<=27) %>%
  dplyr::group_by(., Age) %>%
  dplyr::summarize(., deaths=sum(deaths))

# vector indicating to what seroprev bin each death group belongs
ageGroups <- c(0, 0, 0, 0, 1, 2, 2, 3, 3, 4, 0, 4, 5, 5, 6, 6, 7, 7, 7, 7)
englandWalesDeaths$bin <- ageGroups
englandWalesDeaths <- dplyr::group_by(englandWalesDeaths, bin) %>%
  dplyr::summarize(., deaths=sum(deaths))

# Total deaths segregated between England and Wales
totalEnglandWalesDeaths <- read.csv("../data/downloaded_datasets/england/wales_england_deaths_total.csv",
                                    stringsAsFactors=FALSE)
totalEnglandWalesDeaths$Date <- lubridate::dmy(totalEnglandWalesDeaths$Date)
totalEnglandWalesDeaths <- dplyr::filter(totalEnglandWalesDeaths, !is.na(Date)
                                         & Date<="2020-07-05")
propDeathsEngland <- with(totalEnglandWalesDeaths, sum(England)/(sum(England)+sum(Wales)))
# Adjust age-stratified deaths
deaths_England <- round(englandWalesDeaths$deaths[2:8]*propDeathsEngland)

# HOSPITALIZATIONS
# https://coronavirus.data.gov.uk/details/download
# Hospital 2 https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-hospital-activity/
hospitalized_England <- read.csv("../data/downloaded_datasets/england/hospital_data.csv",
                                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., date=="2020-07-01" & age %in% c("18_to_64", "65_to_84", "85+"))

age_England_Hosp2 <- c("18-64", "65+")
youngerInd <- which(hospitalized_England$age=="18_to_64")
hospitalized_England2 <- c(hospitalized_England$value[youngerInd],
                           sum(hospitalized_England$value[-youngerInd]))

# Estimate seroprevalence for hospitalization age bins
cases_England <- round(population_seroprev_England*seroprev_England/100)
casesL_England <- round(population_seroprev_England*seroprevL_England/100)
casesH_England <- round(population_seroprev_England*seroprevH_England/100)

cases_England <- c(sum(cases_England[1:5]), sum(cases_England[6:7]))
casesL_England <- c(sum(casesL_England[1:5]), sum(casesL_England[6:7]))
casesH_England <- c(sum(casesH_England[1:5]), sum(casesH_England[6:7]))

population_Hosp_England <- c(sum(population_seroprev_England[1:5]),
                             sum(population_seroprev_England[6:7]))

# seroprevalence for hospitalization bins
seroprev_England2 <- signif(cases_England/population_Hosp_England*100, digits=2)
seroprevL_England2 <- signif(casesL_England/population_Hosp_England*100, digits=2)
seroprevH_England2 <- signif(casesH_England/population_Hosp_England*100, digits=2)


### FULL TABLE

England_deaths_w1 <- data.frame(Age=age_England_seroprev,
                       Population=population_seroprev_England,
                       Prevalence=seroprev_England,
                       PrevalenceL=seroprevL_England,
                       PrevalenceH=seroprevH_England,
                       Prevalence_unvacc=seroprev_England,
                       PrevalenceL_unvacc=seroprevL_England,
                       PrevalenceH_unvacc=seroprevH_England,
                       Severe=NA,
                       Critical=NA,
                       Deaths=deaths_England,
                       Deaths_unvacc=deaths_England,
                       onePlusVaccine=0,
                       twoPlusVaccine=0,
                       Type="Seroprevalence",
                       Location="England",
                       EndPointOutcome="2020-07-05",
                       MidPointCases="2020-07-01",
                       Period=1)

England_hospital_w1 <- data.frame(Age=age_England_Hosp2,
                       Population=population_Hosp_England,
                       Prevalence=seroprev_England2,
                       PrevalenceL=seroprevL_England2,
                       PrevalenceH=seroprevH_England2,
                       Prevalence_unvacc=seroprev_England2,
                       PrevalenceL_unvacc=seroprevL_England2,
                       PrevalenceH_unvacc=seroprevH_England2,
                       Severe=hospitalized_England2,
                       Critical=NA,
                       Deaths=NA,
                       Deaths_unvacc=NA,
                       onePlusVaccine=0,
                       twoPlusVaccine=0,
                       Type="Seroprevalence",
                       Location="England",
                       EndPointOutcome="2020-07-01",
                       MidPointCases="2020-07-01",
                       Period=1)


################
### Time 2
################

### POPULATION
age_England_Pop <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
                     "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                     "60-64", "65-69", "70-74", "75-79", "80-84", "85+")
population_England <- c(3496750, 3135711, 3258677, 2079229, 5267401,
                        3836609, 3683915, 3732161, 4099089, 4100526,
                        3601694, 3183915, 3377162, 2674161, 2178672,
                        1777547, 1338005, 1254688)

population_seroprev_England <- c(round(population_England[4]*2/5)+population_England[5],
                                 sum(population_England[6:7]),
                                 sum(population_England[8:9]),
                                 sum(population_England[10:11]),
                                 sum(population_England[12:13]),
                                 sum(population_England[14:15]),
                                 sum(population_England[16:18]))


### SEROPREVALENCE
# Seroprevalence from:
# Increasing SARS-CoV-2 antibody prevalence in England at the start of the second
# wave: REACT-2 Round 4 cross-sectional study in 160,000 adults
# (supplementary)
age_England_seroprev <- c("18-24", "25-34", "35-44", "45-54", "55-64",
                       "65-74", "75+")
seroprev_England <- c(9.85, 6.60, 5.28, 5.48, 5.53, 3.35, 2.63)
seroprevL_England <- c(9.33, 6.25, 4.94, 5.16, 5.17, 3.03, 2.30)
seroprevH_England <- c(10.38, 6.96, 5.63, 5.82, 5.89, 3.70, 2.98)

### OUTCOMES
# DEATHS
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/weeklyprovisionalfiguresondeathsregisteredinenglandandwales
englandWalesDeaths <- read.csv("../data/downloaded_datasets/england/weekly_deaths_2020_ONS.csv",
                    stringsAsFactors=FALSE) %>%
  tidyr::pivot_longer(., cols=c(2:54), names_to="week", values_to="deaths") %>%
  dplyr::mutate(., week=as.integer(gsub("X", "", week))) %>%
  dplyr::filter(., week<=46) %>%
  dplyr::group_by(., Age) %>%
  dplyr::summarize(., deaths=sum(deaths))
# vector indicating to what seroprev bin each death group belongs
ageGroups <- c(0, 0, 0, 0, 1, 2, 2, 3, 3, 4, 0, 4, 5, 5, 6, 6, 7, 7, 7, 7)

englandWalesDeaths$bin <- ageGroups
englandWalesDeaths <- dplyr::group_by(englandWalesDeaths, bin) %>%
  dplyr::summarize(., deaths=sum(deaths))

# Total deaths segregated between England and Wales
totalEnglandWalesDeaths <- read.csv("../data/downloaded_datasets/england/wales_england_deaths_total.csv",
                                    stringsAsFactors=FALSE)
totalEnglandWalesDeaths$Date <- lubridate::dmy(totalEnglandWalesDeaths$Date)
totalEnglandWalesDeaths <- dplyr::filter(totalEnglandWalesDeaths, !is.na(Date)
                                         & Date<="2020-11-03")
propDeathsEngland <- with(totalEnglandWalesDeaths, sum(England)/(sum(England)+sum(Wales)))
# Adjust age-stratified deaths
deaths_England <- round(englandWalesDeaths$deaths[2:8]*propDeathsEngland)

# HOSPITALIZATIONS
# https://coronavirus.data.gov.uk/details/download
# Hospital 2 https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-hospital-activity/
hospitalized_England <- read.csv("../data/downloaded_datasets/england/hospital_data.csv",
                                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., date=="2020-11-03" & age %in% c("18_to_64", "65_to_84", "85+"))

age_England_Hosp2 <- c("18-64", "65+")
youngerInd <- which(hospitalized_England$age=="18_to_64")
hospitalized_England2 <- c(hospitalized_England$value[youngerInd],
                           sum(hospitalized_England$value[-youngerInd]))

# Estimate seroprevalence for hospitalization age bins
cases_England <- round(population_seroprev_England*seroprev_England/100)
casesL_England <- round(population_seroprev_England*seroprevL_England/100)
casesH_England <- round(population_seroprev_England*seroprevH_England/100)

cases_England <- c(sum(cases_England[1:5]), sum(cases_England[6:7]))
casesL_England <- c(sum(casesL_England[1:5]), sum(casesL_England[6:7]))
casesH_England <- c(sum(casesH_England[1:5]), sum(casesH_England[6:7]))

population_Hosp_England <- c(sum(population_seroprev_England[1:5]),
                             sum(population_seroprev_England[6:7]))

# seroprevalence for hospitalization bins
seroprev_England2 <- signif(cases_England/population_Hosp_England*100, digits=2)
seroprevL_England2 <- signif(casesL_England/population_Hosp_England*100, digits=2)
seroprevH_England2 <- signif(casesH_England/population_Hosp_England*100, digits=2)


### FULL TABLE
England_deaths_w2 <- data.frame(Age=age_England_seroprev,
                       Population=population_seroprev_England,
                       Prevalence=seroprev_England,
                       PrevalenceL=seroprevL_England,
                       PrevalenceH=seroprevH_England,
                       Prevalence_unvacc=seroprev_England,
                       PrevalenceL_unvacc=seroprevL_England,
                       PrevalenceH_unvacc=seroprevH_England,
                       Severe=NA,
                       Critical=NA,
                       Deaths=deaths_England,
                       Deaths_unvacc=deaths_England,
                       onePlusVaccine=0,
                       twoPlusVaccine=0,
                       Type="Seroprevalence",
                       Location="England",
                       EndPointOutcome="2020-11-03",
                       MidPointCases="2020-11-03",
                       Period=2)

England_hospital_w2 <- data.frame(Age=age_England_Hosp2,
                       Population=population_Hosp_England,
                       Prevalence=seroprev_England2,
                       PrevalenceL=seroprevL_England2,
                       PrevalenceH=seroprevH_England2,
                       Prevalence_unvacc=seroprev_England2,
                       PrevalenceL_unvacc=seroprevL_England2,
                       PrevalenceH_unvacc=seroprevH_England2,
                       Severe=hospitalized_England2,
                       Critical=NA,
                       Deaths=NA,
                       Deaths_unvacc=NA,
                       onePlusVaccine=0,
                       twoPlusVaccine=0,
                       Type="Seroprevalence",
                       Location="England",
                       EndPointOutcome="2020-11-03",
                       MidPointCases="2020-11-03",
                       Period=2)


################
### Time 3
################

### POPULATION
age_England_Pop <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
                     "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                     "60-64", "65-69", "70-74", "75-79", "80-84", "85+")
population_England <- c(3496750, 3135711, 3258677, 2079229, 5267401,
                        3836609, 3683915, 3732161, 4099089, 4100526,
                        3601694, 3183915, 3377162, 2674161, 2178672,
                        1777547, 1338005, 1254688)

population_seroprev_England <- c(round(population_England[4]*2/5)+population_England[5],
                                 sum(population_England[6:7]),
                                 sum(population_England[8:9]),
                                 sum(population_England[10:11]),
                                 sum(population_England[12:13]),
                                 sum(population_England[14:15]),
                                 sum(population_England[16:18]))


### SEROPREVALENCE
# Seroprevalence from:
# REACT-2 Round 5: increasing prevalence of SARS-CoV-2 antibodies
# demonstrate impact of the second wave and of vaccine roll-out in England

### Unvaccinated seroprev may be obtainable from the data in the paper
# Combine Table 4 with numbers of vaccinated positives with the raw
# test results
age_England_seroprev <- c("18-29", "30-39", "40-49", "50-59", "60-69",
                       "70-79", "80+")
seroprev_England <- c(14.5, 9.9, 9.6, 9.1, 6.9, 4.5, 10.3)
seroprevL_England <- c(14.0, 9.5, 9.2, 8.7, 6.5, 4.1, 7.4)
seroprevH_England <- c(15.0, 10.4, 10.1, 9.6, 7.3, 5.1, 14.1)

### OUTCOMES
# DEATHS
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/weeklyprovisionalfiguresondeathsregisteredinenglandandwales
englandWalesDeaths2020 <- read.csv("../data/downloaded_datasets/england/weekly_deaths_2020_ONS.csv",
                    stringsAsFactors=FALSE) %>%
  tidyr::pivot_longer(., cols=c(2:54), names_to="week", values_to="deaths") %>%
  dplyr::mutate(., week=as.integer(gsub("X", "", week))) %>%
  dplyr::group_by(., Age) %>%
  dplyr::summarize(., deaths=sum(deaths))
englandWalesDeaths2021 <- read.csv("../data/downloaded_datasets/england/weekly_deaths_2021_ONS.csv",
                    stringsAsFactors=FALSE) %>%
  tidyr::pivot_longer(., cols=c(2:54), names_to="week", values_to="deaths") %>%
  dplyr::mutate(., week=as.integer(gsub("X", "", week))) %>%
  dplyr::filter(., week<=04) %>%
  dplyr::group_by(., Age) %>%
  dplyr::summarize(., deaths=sum(deaths))
# add 2020 and 2021 deaths
englandWalesDeaths2021$deaths <- englandWalesDeaths2020$deaths + englandWalesDeaths2021$deaths

# vector indicating to what seroprev bin each death group belongs
ageGroups <- c(0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 0, 4, 4, 5, 5, 6, 6, 7, 7, 7)
englandWalesDeaths2021$bin <- ageGroups
englandWalesDeaths2021 <- dplyr::group_by(englandWalesDeaths2021, bin) %>%
  dplyr::summarize(., deaths=sum(deaths))

# Total deaths segregated between England and Wales
totalEnglandWalesDeaths <- read.csv("../data/downloaded_datasets/england/wales_england_deaths_total.csv",
                                    stringsAsFactors=FALSE)
totalEnglandWalesDeaths$Date <- lubridate::dmy(totalEnglandWalesDeaths$Date)
totalEnglandWalesDeaths <- dplyr::filter(totalEnglandWalesDeaths, !is.na(Date)
                                         & Date<="2021-02-01")
propDeathsEngland <- with(totalEnglandWalesDeaths, sum(England)/(sum(England)+sum(Wales)))
# Adjust age-stratified deaths
deaths_England <- round(englandWalesDeaths2021$deaths[2:8]*propDeathsEngland)

# DEATHS BY VACCINATION STATUS
# Source: https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/deathsbyvaccinationstatusengland
# September file
deathsVaccine <- read.csv("../data/downloaded_datasets/england/deaths_by_vaccination_2021_ONS.csv",
                    stringsAsFactors=FALSE)

unvaccinatedStrings <- c("Unvaccinated", "Within 21 days of first dose")
unvaccinatedDeathsEngland <- deathsVaccine %>%
  dplyr::filter(., Vaccination.status %in% unvaccinatedStrings & Week.number<=4) %>%
  dplyr::group_by(., Age.group) %>%
  dplyr::summarize(., deaths=sum(Number.of.deaths))

vaccinatedDeathsEngland <- deathsVaccine %>%
  dplyr::filter(., !(Vaccination.status %in% unvaccinatedStrings) & Week.number<=4) %>%
  dplyr::group_by(., Age.group) %>%
  dplyr::summarize(., deaths=sum(Number.of.deaths))
# Match seroprevalence bins
relIFR1 <- levinIFR(seq(10, 59))/sum(levinIFR(seq(10, 59)))
relIFR1 <- c(sum(relIFR1[9:20]), sum(relIFR1[21:30]),
             sum(relIFR1[31:40]), sum(relIFR1[41:50]))
vaccinatedDeathsEngland <- c(round(vaccinatedDeathsEngland$deaths[1]*relIFR1),
                             vaccinatedDeathsEngland$deaths[2:4])


### VACCINATION
# https://coronavirus.data.gov.uk/details/download
# Check also: https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19antibodydatafortheuk/2021

vaccinationEngland <- read.csv("../data/downloaded_datasets/england/vaccination_by_age_england.csv",
                    stringsAsFactors=FALSE) %>%
  dplyr::filter(., date=="2021-02-01")
ageGroups <- c(0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7)
ageGroupsString <- c("0-17", age_England_seroprev)[ageGroups+1]
vaccinationEngland$ageGroups <- ageGroupsString


administeredVaccinesEngland <- dplyr::group_by(vaccinationEngland, ageGroups) %>%
  dplyr::summarize(., firstDose=sum(cumPeopleVaccinatedFirstDoseByVaccinationDate),
                    secondDose=sum(cumPeopleVaccinatedSecondDoseByVaccinationDate),
                    thirdDose=sum(cumPeopleVaccinatedThirdInjectionByVaccinationDate))
populationVaccEngland <- dplyr::group_by(vaccinationEngland, ageGroups) %>%
  dplyr::summarize(., vaccPop=sum(VaccineRegisterPopulationByVaccinationDate))
vaccPercentageEngland <- signif(administeredVaccinesEngland[,2:4]/
                                populationVaccEngland$vaccPop*100, digits=3)

### check: https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1039677/Vaccine_surveillance_report_-_week_49.pdf
### Check vaccine surveillance reports. Need to add hospitalizations from al reports?


# HOSPITALIZATIONS
# https://coronavirus.data.gov.uk/details/download
# Hospital 2 https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-hospital-activity/
hospitalized_England <- read.csv("../data/downloaded_datasets/england/hospital_data.csv",
                                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., date=="2021-02-01" & age %in% c("18_to_64", "65_to_84", "85+"))

age_England_Hosp2 <- c("18-64", "65+")
youngerInd <- which(hospitalized_England$age=="18_to_64")
hospitalized_England2 <- c(hospitalized_England$value[youngerInd],
                           sum(hospitalized_England$value[-youngerInd]))

# Estimate seroprevalence for hospitalization age bins
cases_England <- round(population_seroprev_England*seroprev_England/100)
casesL_England <- round(population_seroprev_England*seroprevL_England/100)
casesH_England <- round(population_seroprev_England*seroprevH_England/100)

cases_England <- c(sum(cases_England[1:5]), sum(cases_England[6:7]))
casesL_England <- c(sum(casesL_England[1:5]), sum(casesL_England[6:7]))
casesH_England <- c(sum(casesH_England[1:5]), sum(casesH_England[6:7]))

population_Hosp_England <- c(sum(population_seroprev_England[1:5]),
                             sum(population_seroprev_England[6:7]))

# seroprevalence for hospitalization bins
seroprev_England2 <- signif(cases_England/population_Hosp_England*100, digits=2)
seroprevL_England2 <- signif(casesL_England/population_Hosp_England*100, digits=2)
seroprevH_England2 <- signif(casesH_England/population_Hosp_England*100, digits=2)

populationVaccEngland2 <- c(sum(populationVaccEngland$vaccPop[2:5])+
                         round(populationVaccEngland$vaccPop[6]/2),
                         round(populationVaccEngland$vaccPop[6]/2+
                             sum(populationVaccEngland$vaccPop[7:8])))

# Adjust vaccine stratification to hospital data
administeredVaccinesEngland_firstDose2 <- signif(c(sum(administeredVaccinesEngland$firstDose[2:5])+
                         round(administeredVaccinesEngland$firstDose[6]/2),
                         round(administeredVaccinesEngland$firstDose[6]/2+
                               sum(administeredVaccinesEngland$firstDose[7:8])))/
                         populationVaccEngland2*100, digits=3)
administeredVaccinesEngland_secondDose2 <- signif(c(sum(administeredVaccinesEngland$secondDose[2:5])+
                         round(administeredVaccinesEngland$secondDose[6]/2),
                         round(administeredVaccinesEngland$secondDose[6]/2+
                               sum(administeredVaccinesEngland$secondDose[7:8])))/
                         populationVaccEngland2*100, digits=3)

### FULL TABLE
England_deaths_w3 <- data.frame(Age=age_England_seroprev,
                       Population=population_seroprev_England,
                       Prevalence=seroprev_England,
                       PrevalenceL=seroprevL_England,
                       PrevalenceH=seroprevH_England,
                       Prevalence_unvacc=NA,
                       PrevalenceL_unvacc=NA,
                       PrevalenceH_unvacc=NA,
                       Severe=NA,
                       Critical=NA,
                       Deaths=deaths_England,
                       Deaths_unvacc=deaths_England-vaccinatedDeathsEngland,
                       onePlusVaccine=vaccPercentageEngland$firstDose[-1],
                       twoPlusVaccine=vaccPercentageEngland$secondDose[-1],
                       Type="Seroprevalence",
                       Location="England",
                       EndPointOutcome="2021-02-01",
                       MidPointCases="2021-02-01",
                       Period=3)

England_hospital_w3 <- data.frame(Age=age_England_Hosp2,
                       Population=population_Hosp_England,
                       Prevalence=seroprev_England2,
                       PrevalenceL=seroprevL_England2,
                       PrevalenceH=seroprevH_England2,
                       Prevalence_unvacc=NA,
                       PrevalenceL_unvacc=NA,
                       PrevalenceH_unvacc=NA,
                       Severe=hospitalized_England2,
                       Critical=NA,
                       Deaths=NA,
                       Deaths_unvacc=NA,
                       onePlusVaccine=administeredVaccinesEngland_firstDose2,
                       twoPlusVaccine=administeredVaccinesEngland_secondDose2,
                       Type="Seroprevalence",
                       Location="England",
                       EndPointOutcome="2021-02-01",
                       MidPointCases="2021-02-01",
                       Period=3)



################
### Time 4
################

### POPULATION
age_England_Pop <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
                     "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                     "60-64", "65-69", "70-74", "75-79", "80-84", "85+")
population_England <- c(3496750, 3135711, 3258677, 2079229, 5267401,
                        3836609, 3683915, 3732161, 4099089, 4100526,
                        3601694, 3183915, 3377162, 2674161, 2178672,
                        1777547, 1338005, 1254688)

population_seroprev_England <- c(round(population_England[4]*2/5)+population_England[5],
                                 sum(population_England[6:7]),
                                 sum(population_England[8:9]),
                                 sum(population_England[10:11]),
                                 sum(population_England[12:13]),
                                 sum(population_England[14:15]),
                                 sum(population_England[16:18]))


### SEROPREVALENCE
# Seroprevalence from:
# Vaccine uptake and SARS-CoV-2 antibody prevalence among
# 207,337 adults during May 2021 in England: REACT-2 
## Unvaccinated
age_England_seroprev <- c("18-24", "25-34", "35-44", "45-54", "55-64",
                       "65-74", "75+")
seroprev_England_unvacc <- c(17.0, 12.0, 14.5, 27.8, 29.4, 24.8, 41.0)
seroprevL_England_unvacc <- c(16.3, 11.5, 13.8, 25.4, 25.4, 18.8, 32.9)
seroprevH_England_unvacc <- c(17.6, 12.5, 15.3, 30.3, 33.8, 31.9, 50.0)

# All
age_England_seroprev <- c("18-24", "25-34", "35-44", "45-54", "55-64",
                       "65-74", "75+")
seroprev_England <- c(35.8, 37.5, 48.4, 66.6, 66.2, 90.2, 95.3)
seroprevL_England <- c(35.1, 36.9, 47.8, 66.0, 65.5, 89.6, 94.6)
seroprevH_England <- c(36.5, 38.0, 49.0, 67.2, 66.9, 90.9, 95.9)


### OUTCOMES
# DEATHS
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/weeklyprovisionalfiguresondeathsregisteredinenglandandwales
englandWalesDeaths2020 <- read.csv("../data/downloaded_datasets/england/weekly_deaths_2020_ONS.csv",
                    stringsAsFactors=FALSE) %>%
  tidyr::pivot_longer(., cols=c(2:54), names_to="week", values_to="deaths") %>%
  dplyr::mutate(., week=as.integer(gsub("X", "", week))) %>%
  dplyr::group_by(., Age) %>%
  dplyr::summarize(., deaths=sum(deaths))
englandWalesDeaths2021 <- read.csv("../data/downloaded_datasets/england/weekly_deaths_2021_ONS.csv",
                    stringsAsFactors=FALSE) %>%
  tidyr::pivot_longer(., cols=c(2:54), names_to="week", values_to="deaths") %>%
  dplyr::mutate(., week=as.integer(gsub("X", "", week))) %>%
  dplyr::filter(., week<=19) %>%
  dplyr::group_by(., Age) %>%
  dplyr::summarize(., deaths=sum(deaths))
englandWalesDeaths2021$deaths <- englandWalesDeaths2020$deaths + englandWalesDeaths2021$deaths
# Total deaths segregated between England and Wales
totalEnglandWalesDeaths <- read.csv("../data/downloaded_datasets/england/wales_england_deaths_total.csv",
                                    stringsAsFactors=FALSE)
totalEnglandWalesDeaths$Date <- lubridate::dmy(totalEnglandWalesDeaths$Date)
totalEnglandWalesDeaths <- dplyr::filter(totalEnglandWalesDeaths, !is.na(Date)
                                         & Date<="2021-02-01")
propDeathsEngland <- with(totalEnglandWalesDeaths, sum(England)/(sum(England)+sum(Wales)))
# Adjust age-stratified deaths
deaths_England <- round(englandWalesDeaths2021$deaths[2:8]*propDeathsEngland)


# vector indicating to what seroprev bin each death group belongs
ageGroups <- c(0, 0, 0, 0, 1, 2, 2, 3, 3, 4, 0, 4, 5, 5, 6, 6, 7, 7, 7, 7)
englandWalesDeaths2021$bin <- ageGroups
englandWalesDeaths2021 <- dplyr::group_by(englandWalesDeaths2021, bin) %>%
  dplyr::summarize(., deaths=sum(deaths))
deaths_England <- englandWalesDeaths2021$deaths[2:8]
#### HAVE TO CORRECT FOR WALES DEATHS

# DEATHS BY VACCINATION STATUS
# Source: https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/deathsbyvaccinationstatusengland
# September file
deathsVaccine <- read.csv("../data/downloaded_datasets/england/deaths_by_vaccination_2021_ONS.csv",
                    stringsAsFactors=FALSE)
unvaccinatedStrings <- c("Unvaccinated", "Within 21 days of first dose")
unvaccinatedDeathsEngland <- deathsVaccine %>%
  dplyr::filter(., Vaccination.status %in% unvaccinatedStrings & Week.number<=19) %>%
  dplyr::group_by(., Age.group) %>%
  dplyr::summarize(., deaths=sum(Number.of.deaths))
vaccinatedDeathsEngland <- deathsVaccine %>%
  dplyr::filter(., !(Vaccination.status %in% unvaccinatedStrings) & Week.number<=19) %>%
  dplyr::group_by(., Age.group) %>%
  dplyr::summarize(., deaths=sum(Number.of.deaths))
# Match seroprevalence bins
relIFR1 <- levinIFR(seq(10, 59))/sum(levinIFR(seq(10, 59)))
relIFR1 <- c(sum(relIFR1[9:20]), sum(relIFR1[21:30]),
             sum(relIFR1[31:40]), sum(relIFR1[41:50]))
vaccinatedDeathsEngland <- c(round(vaccinatedDeathsEngland$deaths[1]*relIFR1),
                             vaccinatedDeathsEngland$deaths[2:4])


### VACCINATION
# Source: https://coronavirus.data.gov.uk/details/download
# Check also: https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19antibodydatafortheuk/2021
vaccinationEngland <- read.csv("../data/downloaded_datasets/england/vaccination_by_age_england.csv",
                    stringsAsFactors=FALSE) %>%
  dplyr::filter(., date=="2021-05-18")
ageGroups <- c(0, 0, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 7)
ageGroupsString <- c("0-17", age_England_seroprev)[ageGroups+1]
vaccinationEngland$ageGroups <- ageGroupsString

administeredVaccinesEngland <- dplyr::group_by(vaccinationEngland, ageGroups) %>%
  dplyr::summarize(., firstDose=sum(cumPeopleVaccinatedFirstDoseByVaccinationDate),
                    secondDose=sum(cumPeopleVaccinatedSecondDoseByVaccinationDate),
                    thirdDose=sum(cumPeopleVaccinatedThirdInjectionByVaccinationDate))
populationVaccEngland <- dplyr::group_by(vaccinationEngland, ageGroups) %>%
  dplyr::summarize(., vaccPop=sum(VaccineRegisterPopulationByVaccinationDate))
vaccPercentageEngland <- signif(administeredVaccinesEngland[,2:4]/
                                populationVaccEngland$vaccPop*100, digits=3)



# HOSPITALIZATIONS
# https://coronavirus.data.gov.uk/details/download
# Hospital 2 https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-hospital-activity/
hospitalized_England <- read.csv("../data/downloaded_datasets/england/hospital_data.csv",
                                 stringsAsFactors=FALSE) %>%
  dplyr::filter(., date=="2021-05-18" & age %in% c("18_to_64", "65_to_84", "85+"))

age_England_Hosp2 <- c("18-64", "65+")
youngerInd <- which(hospitalized_England$age=="18_to_64")
hospitalized_England2 <- c(hospitalized_England$value[youngerInd],
                           sum(hospitalized_England$value[-youngerInd]))

# Adjust vaccine stratification to hospital data
populationVaccEngland2 <- c(sum(populationVaccEngland$vaccPop[2:5])+
                         round(populationVaccEngland$vaccPop[6]/2),
                         round(populationVaccEngland$vaccPop[6]/2+
                             sum(populationVaccEngland$vaccPop[7:8])))

administeredVaccinesEngland_firstDose2 <- signif(c(sum(administeredVaccinesEngland$firstDose[2:5])+
                         round(administeredVaccinesEngland$firstDose[6]/2),
                         round(administeredVaccinesEngland$firstDose[6]/2+
                               sum(administeredVaccinesEngland$firstDose[7:8])))/
                         populationVaccEngland2*100, digits=3)
administeredVaccinesEngland_secondDose2 <- signif(c(sum(administeredVaccinesEngland$secondDose[2:5])+
                         round(administeredVaccinesEngland$secondDose[6]/2),
                         round(administeredVaccinesEngland$secondDose[6]/2+
                               sum(administeredVaccinesEngland$secondDose[7:8])))/
                         populationVaccEngland2*100, digits=3)


# Estimate seroprevalence for hospitalization age bins
cases_England <- round(population_seroprev_England*seroprev_England/100)
casesL_England <- round(population_seroprev_England*seroprevL_England/100)
casesH_England <- round(population_seroprev_England*seroprevH_England/100)
cases_England <- c(sum(cases_England[1:5]), sum(cases_England[6:7]))
casesL_England <- c(sum(casesL_England[1:5]), sum(casesL_England[6:7]))
casesH_England <- c(sum(casesH_England[1:5]), sum(casesH_England[6:7]))

unvaccPop <- population_seroprev_England-administeredVaccinesEngland$firstDose[-1]
cases_England_unvacc <- round(unvaccPop*seroprev_England_unvacc/100)
casesL_England_unvacc <- round(unvaccPop*seroprevL_England_unvacc/100)
casesH_England_unvacc <- round(unvaccPop*seroprevH_England_unvacc/100)
cases_England_unvacc <- c(sum(cases_England_unvacc[1:5]), sum(cases_England_unvacc[6:7]))
casesL_England_unvacc <- c(sum(casesL_England_unvacc[1:5]), sum(casesL_England_unvacc[6:7]))
casesH_England_unvacc <- c(sum(casesH_England_unvacc[1:5]), sum(casesH_England_unvacc[6:7]))
# seroprevalence for hospitalization bins
seroprev_England2 <- signif(cases_England/population_Hosp_England*100, digits=3)
seroprevL_England2 <- signif(casesL_England/population_Hosp_England*100, digits=3)
seroprevH_England2 <- signif(casesH_England/population_Hosp_England*100, digits=3)

population_Hosp_England <- c(sum(population_seroprev_England[1:5]),
                             sum(population_seroprev_England[6:7]))
# unvacc seroprevalence
unvaccPopHosp <- population_Hosp_England - administeredVaccinesEngland_firstDose2
seroprev_England2_unvacc <- signif(cases_England_unvacc/unvaccPopHosp*100, digits=3)
seroprevL_England2_unvacc <- signif(casesL_England_unvacc/unvaccPopHosp*100, digits=3)
seroprevH_England2_unvacc <- signif(casesH_England_unvacc/unvaccPopHosp*100, digits=3)


### FULL TABLE
England_deaths_w4 <- data.frame(Age=age_England_seroprev,
                       Population=population_seroprev_England,
                       Prevalence=seroprev_England,
                       PrevalenceL=seroprevL_England,
                       PrevalenceH=seroprevH_England,
                       Prevalence_unvacc=seroprev_England_unvacc,
                       PrevalenceL_unvacc=seroprevL_England_unvacc,
                       PrevalenceH_unvacc=seroprevH_England_unvacc,
                       Severe=NA,
                       Critical=NA,
                       Deaths=deaths_England,
                       Deaths_unvacc=deaths_England-vaccinatedDeathsEngland,
                       onePlusVaccine=vaccPercentageEngland$firstDose[-1],
                       twoPlusVaccine=vaccPercentageEngland$secondDose[-1],
                       Type="Seroprevalence",
                       Location="England",
                       EndPointOutcome="2021-05-18",
                       MidPointCases="2021-05-18",
                       Period=4)

England_hospital_w4 <- data.frame(Age=age_England_Hosp2,
                       Population=population_Hosp_England,
                       Prevalence=seroprev_England2,
                       PrevalenceL=seroprevL_England2,
                       PrevalenceH=seroprevH_England2,
                       Prevalence_unvacc=seroprev_England2_unvacc,
                       PrevalenceL_unvacc=seroprevL_England2_unvacc,
                       PrevalenceH_unvacc=seroprevH_England2_unvacc,
                       Severe=hospitalized_England2,
                       Critical=NA,
                       Deaths=NA,
                       Deaths_unvacc=NA,
                       onePlusVaccine=administeredVaccinesEngland_firstDose2,
                       twoPlusVaccine=administeredVaccinesEngland_secondDose2,
                       Type="Seroprevalence",
                       Location="England",
                       EndPointOutcome="2021-05-18",
                       MidPointCases="2021-05-18",
                       Period=4)


################
### Scotland
################

# deaths https://www.nrscotland.gov.uk/statistics-and-data/statistics/statistics-by-theme/vital-events/general-publications/weekly-and-monthly-data-on-births-and-deaths/deaths-involving-coronavirus-covid-19-in-scotland
# deaths https://statistics.gov.scot/slice?dataset=http%3A%2F%2Fstatistics.gov.scot%2Fdata%2Fdeaths-involving-coronavirus-covid-19&http%3A%2F%2Fstatistics.gov.scot%2Fdef%2Fdimension%2Fage=http%3A%2F%2Fstatistics.gov.scot%2Fdef%2Fconcept%2Fage%2F1-14-years

#############
### TIME 1
#############
### POPULATION
scotlandPopAge <- c("0-5", "5-9", "10-14", "15-19", "20-24", "25-29",
                   "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                   "60-64", "65-69", "70-74", "75-79", "80-84", "85-89",
                    "90+")
population_Scotland <- c(263806, 297903, 298081, 282120, 341755, 377204,
                         374069, 355666, 324366, 349924, 393113, 399344,
                         352569, 300433, 285830, 198210, 143296, 84562, 13927)

population_Scotland2 <- c(sum(population_Scotland[4:5]),
                         sum(population_Scotland[6:7]), sum(population_Scotland[8:10]),
                         sum(population_Scotland[11:13]), sum(population_Scotland[14:15]),
                         sum(population_Scotland[16:19]))

### SEROPREVALENCE
# https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19antibodydatafortheuk/2021
age_Scotland_seroprev <- c("16-24", "25-34", "35-49", "50-59", "60-64", "65-69",
                         "70-74", "75-79", "80+")
seroprev_Scotland <- c(16.0, 9.7, 6.4, 6.3, 5.5, 5.8, 3.7, 2.8, 2.0)
seroprevL_Scotland <- c(11.6, 7.0, 4.7, 4.6, 3.9, 4.1, 2.6, 1.8, 1.3)
seroprevH_Scotland <- c(21.9, 13.2, 8.8, 8.8, 7.7, 8.2, 5.4, 4.1, 3.1)
# change some seroprevalence bins to match death data. Make new
# seroprevalence bins by weighting originals by bin-population size
newBinsAges <- c("50-64", "65-74", "75+")
relPop1 <- c(sum(population_Scotland[11:12]), population_Scotland[13])/
  sum(population_Scotland[11:13])
relPop2 <- c(population_Scotland[14:15])/sum(population_Scotland[14:15])
relPop3 <- c(population_Scotland[16], sum(population_Scotland[17:19]))/
  sum(population_Scotland[16:19])
seroprevNew <- c(sum(seroprev_Scotland[4:5]*relPop1),
  sum(seroprev_Scotland[6:7]*relPop2), sum(seroprev_Scotland[8:9]*relPop3)) 
seroprevLNew <- c(sum(seroprevL_Scotland[4:5]*relPop1),
  sum(seroprevL_Scotland[6:7]*relPop2), sum(seroprevL_Scotland[8:9]*relPop3)) 
seroprevHNew <- c(sum(seroprevH_Scotland[4:5]*relPop1),
  sum(seroprevH_Scotland[6:7]*relPop2), sum(seroprevH_Scotland[8:9]*relPop3)) 

# New seroprevalence bins
age_Scotland_seroprev <- c("16-24", "25-34", "35-49", "50-64", "65-74", "75+")
seroprev_Scotland <- signif(c(seroprev_Scotland[1:3], seroprevNew), digits=3)
seroprevL_Scotland <- signif(c(seroprevL_Scotland[1:3], seroprevLNew), digits=3)
seroprevH_Scotland <- signif(c(seroprevH_Scotland[1:3], seroprevHNew), digits=3)

### OUTCOMES
# https://www.nrscotland.gov.uk/statistics-and-data/statistics/statistics-by-theme/vital-events/general-publications/weekly-and-monthly-data-on-births-and-deaths/deaths-involving-coronavirus-covid-19-in-scotland/archive
# https://www.nrscotland.gov.uk/files//statistics/covid19/covid-deaths-report-week-53.pdf
dataScotland <- read.csv("../data/downloaded_datasets/scotland/deaths2020.csv",
                         stringsAsFactors=FALSE) %>%
  tidyr::pivot_longer(., cols=c(2:54), names_to="week", values_to="deaths") %>%
  dplyr::mutate(., week=as.integer(gsub("X", "", week))) %>%
  dplyr::filter(., week<=46) %>%
  dplyr::group_by(., Age) %>%
  dplyr::summarize(., deaths=sum(deaths))
# Adjust death bins to match seroprevalence
relIFR1 <- levinIFR(seq(15,44))/sum(levinIFR(seq(15,44)))
relIFR1 <- c(sum(relIFR1[1:10]), sum(relIFR1[11:20]), sum(relIFR1[21:30]))
relIFR2 <- levinIFR(seq(45,64))/sum(levinIFR(seq(45,64)))
relIFR2 <- c(sum(relIFR2[1:5]), sum(relIFR2[6:20]))

deaths_Scotland <- c(round(dataScotland$deaths[2]*relIFR1[1:2]),
                     round(dataScotland$deaths[2]*relIFR1[3]+
                           dataScotland$deaths[3]*relIFR2[1]),
                     round(dataScotland$deaths[3]*relIFR2[2]),
                     dataScotland$deaths[4],
                     sum(dataScotland$deaths[5:6]))

### FULL TABLE
Scotland_w1 <- data.frame(Age=age_Scotland_seroprev,
                  Population=population_Scotland2,
                  Prevalence=seroprev_Scotland,
                  PrevalenceL=seroprevL_Scotland,
                  PrevalenceH=seroprevH_Scotland,
                  Prevalence_unvacc=seroprev_Scotland,
                  PrevalenceL_unvacc=seroprevL_Scotland,
                  PrevalenceH_unvacc=seroprevH_Scotland,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Scotland,
                  Deaths_unvacc=deaths_Scotland,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Scotland",
                  EndPointOutcome="2020-12-20",
                  MidPointCases="2020-12-20",
                  Period=1)


#############
### TIME 2
#############

### POPULATION
scotlandPopAge <- c("0-5", "5-9", "10-14", "15-19", "20-24", "25-29",
                   "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                   "60-64", "65-69", "70-74", "75-79", "80-84", "85-89",
                    "90+")
population_Scotland <- c(263806, 297903, 298081, 282120, 341755, 377204,
                         374069, 355666, 324366, 349924, 393113, 399344,
                         352569, 300433, 285830, 198210, 143296, 84562, 13927)

population_Scotland2 <- c(sum(population_Scotland[4:5]),
                         sum(population_Scotland[6:7]), sum(population_Scotland[8:10]),
                         sum(population_Scotland[11:13]), sum(population_Scotland[14:15]),
                         sum(population_Scotland[16:19]))

### SEROPREVALENCE
# https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19antibodydatafortheuk/2021
age_Scotland_seroprev <- c("16-24", "25-34", "35-49", "50-59", "60-64", "65-69",
                         "70-74", "75-79", "80+")
seroprev_Scotland <- c(30.7, 31.3, 39.9, 81.3, 78.5, 85.6, 90.8, 90.2, 87.5)
seroprevL_Scotland <- c(25.0, 25.8, 34.0, 76.9, 73.2, 81.6, 87.8, 86.6, 82.8)
seroprevH_Scotland <- c(37.5, 37.3, 46.1, 85.2, 82.9, 88.8, 93.1, 92.8, 91.1)

# change some seroprevalence bins to match death data. Make new
# seroprevalence bins by weighting originals by bin-population size
newBinsAges <- c("50-64", "65-74", "75+")
relPop1 <- c(sum(population_Scotland[11:12]), population_Scotland[13])/
  sum(population_Scotland[11:13])
relPop2 <- c(population_Scotland[14:15])/sum(population_Scotland[14:15])
relPop3 <- c(population_Scotland[16], sum(population_Scotland[17:19]))/
  sum(population_Scotland[16:19])
seroprevNew <- c(sum(seroprev_Scotland[4:5]*relPop1),
  sum(seroprev_Scotland[6:7]*relPop2), sum(seroprev_Scotland[8:9]*relPop3)) 
seroprevLNew <- c(sum(seroprevL_Scotland[4:5]*relPop1),
  sum(seroprevL_Scotland[6:7]*relPop2), sum(seroprevL_Scotland[8:9]*relPop3)) 
seroprevHNew <- c(sum(seroprevH_Scotland[4:5]*relPop1),
  sum(seroprevH_Scotland[6:7]*relPop2), sum(seroprevH_Scotland[8:9]*relPop3)) 

# New seroprevalence bins
age_Scotland_seroprev <- c("16-24", "25-34", "35-49", "50-64", "65-74", "75+")
seroprev_Scotland <- signif(c(seroprev_Scotland[1:3], seroprevNew), digits=3)
seroprevL_Scotland <- signif(c(seroprevL_Scotland[1:3], seroprevLNew), digits=3)
seroprevH_Scotland <- signif(c(seroprevH_Scotland[1:3], seroprevHNew), digits=3)

### OUTCOMES
# https://www.nrscotland.gov.uk/statistics-and-data/statistics/statistics-by-theme/vital-events/general-publications/weekly-and-monthly-data-on-births-and-deaths/deaths-involving-coronavirus-covid-19-in-scotland/archive
# https://www.nrscotland.gov.uk/files//statistics/covid19/covid-deaths-report-week-53.pdf
dataScotland2020 <- read.csv("../data/downloaded_datasets/scotland/deaths2020.csv",
                         stringsAsFactors=FALSE) %>%
  tidyr::pivot_longer(., cols=c(2:54), names_to="week", values_to="deaths") %>%
  dplyr::mutate(., week=as.integer(gsub("X", "", week))) %>%
  dplyr::filter(., week<=46) %>%
  dplyr::group_by(., Age) %>%
  dplyr::summarize(., deaths=sum(deaths))
# Adjust death bins to match seroprevalence
relIFR1 <- levinIFR(seq(15,44))/sum(levinIFR(seq(15,44)))
relIFR1 <- c(sum(relIFR1[1:10]), sum(relIFR1[11:20]), sum(relIFR1[21:30]))
relIFR2 <- levinIFR(seq(45,64))/sum(levinIFR(seq(45,64)))
relIFR2 <- c(sum(relIFR2[1:5]), sum(relIFR2[6:20]))

deaths_Scotland2020 <- c(round(dataScotland2020$deaths[2]*relIFR1[1:2]),
                     round(dataScotland2020$deaths[2]*relIFR1[3]+
                           dataScotland2020$deaths[3]*relIFR2[1]),
                     round(dataScotland2020$deaths[3]*relIFR2[2]),
                     dataScotland2020$deaths[4],
                     sum(dataScotland2020$deaths[5:6]))

dataScotland2021 <- read.csv("../data/downloaded_datasets/scotland/deaths2021.csv",
                             stringsAsFactors=FALSE) %>%
  tidyr::pivot_longer(., cols=starts_with("X"), names_to="Age", values_to="Deaths") %>%
  dplyr::filter(., Week.number<=21 & Registration.year==2021) %>%
  dplyr::group_by(., Age) %>%
  dplyr::summarize(., totalDeaths=sum(Deaths))

deaths_Scotland2021 <- c(round(dataScotland2021$totalDeaths[2]*relIFR1[1:2]),
                     round(dataScotland2021$totalDeaths[2]*relIFR1[3]+
                           dataScotland2021$totalDeaths[3]*relIFR2[1]),
                     round(dataScotland2021$totalDeaths[3]*relIFR2[2]),
                     dataScotland2021$totalDeaths[4],
                     sum(dataScotland2021$totalDeaths[5:6]))

# deaths by vaccination status
# https://www.publichealthscotland.scot/publications/covid-19-statistical-report/covid-19-statistical-report-21-july-2021/
ageDeaths_vacc <- c("0-39", "40-49", "50-59", "60-69", "70-79", "80+")
deathsUnvacc <- c(21, 55, 184, 412, 765, 1525)
deathsOneVacc <- c(1, 1, 5, 13, 43, 194)
deathsTwoVacc <- c(0, 1, 1, 5, 26, 31)
ratioVaccDeaths <- (deathsOneVacc+deathsTwoVacc)/(deathsUnvacc+deathsOneVacc+deathsTwoVacc)
deaths_ScotlandVacc2021 <- round(ratioVaccDeaths*deaths_Scotland2021)
deaths_ScotlandUnvacc2021 <- deaths_Scotland2021 - deaths_ScotlandVacc2021

### VACCINATION
ageGroups <- c("12 to 15", "16 to 17", "18 to 29", "30 to 39", "40 to 49",
               "50 to 54", "55 to 59", "60 to 64", "65 to 69", "70 to 74",
               "75 to 79", "80 years and over")
vaccScotlandOneDose <- read.csv("../data/downloaded_datasets/scotland/daily_vacc_age_sex_20220201.csv",
                                stringsAsFactors=FALSE) %>%
  dplyr::filter(., Date==20210509 & Sex=="Total" & Dose=="Dose 1" &
  AgeGroup %in% ageGroups)
vaccinationScotlandOneDose <- signif(pmin(100,
                                c(vaccScotlandOneDose$CumulativePercentCoverage[3:5],
                                mean(vaccScotlandOneDose$CumulativePercentCoverage[6:8]),
                                mean(vaccScotlandOneDose$CumulativePercentCoverage[9:10]),
                                mean(vaccScotlandOneDose$CumulativePercentCoverage[11:12]))),
                                     digits=3)

vaccScotlandTwoDose <- read.csv("../data/downloaded_datasets/scotland/daily_vacc_age_sex_20220201.csv",
                                stringsAsFactors=FALSE) %>%
  dplyr::filter(., Date==20210509 & Sex=="Total" & Dose=="Dose 2" &
  AgeGroup %in% ageGroups)
vaccinationScotlandTwoDose <- signif(pmin(100,
                                c(vaccScotlandTwoDose$CumulativePercentCoverage[3:5],
                                mean(vaccScotlandTwoDose$CumulativePercentCoverage[6:8]),
                                mean(vaccScotlandTwoDose$CumulativePercentCoverage[9:10]),
                                mean(vaccScotlandTwoDose$CumulativePercentCoverage[11:12]))),
                                     digits=3)


### FULL TABLE
Scotland_w2 <- data.frame(Age=age_Scotland_seroprev,
                  Population=population_Scotland2,
                  Prevalence=seroprev_Scotland,
                  PrevalenceL=seroprevL_Scotland,
                  PrevalenceH=seroprevH_Scotland,
                  Prevalence_unvacc=NA,
                  PrevalenceL_unvacc=NA,
                  PrevalenceH_unvacc=NA,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Scotland2020+deaths_Scotland2021,
                  Deaths_unvacc=deaths_Scotland2020+deaths_ScotlandUnvacc2021,
                  onePlusVaccine=vaccinationScotlandOneDose,
                  twoPlusVaccine=vaccinationScotlandTwoDose,
                  Type="Seroprevalence",
                  Location="Scotland",
                  EndPointOutcome="2021-05-09",
                  MidPointCases="2021-05-06",
                  Period=2)

#############
### Wales
#############

# Wales hospital by vaccine data
# https://www2.nphs.wales.nhs.uk/CommunitySurveillanceDocs.nsf/3dc04669c9e1eaa880257062003b246b/a4f536f72da3962b8025875a0031b3c8/$FILE/Survey%20of%20vaccine%20status%20in%20cases%20and%20hospital%20inpatients.pdf


#############
### TIME 1
#############
### POPULATION
walesPopAge <- c("0-5", "5-9", "10-14", "15-19", "20-24", "25-29",
                   "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                   "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90+")
population_Wales <- c(161341, 182178, 184437, 173821, 206557, 208114,
                         196672, 185873, 172930, 192465, 216970, 222222,
                         197416, 179955, 181886, 131023, 90566, 54400, 43749)

population_Wales2 <- c(sum(population_Wales[4:5]),
                         sum(population_Wales[6:7]), sum(population_Wales[8:10]),
                         sum(population_Wales[11:12]), population_Wales[13:16], 
                         sum(population_Wales[17:19]))

### SEROPREVALENCE
# https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19antibodydatafortheuk/2021
age_Wales_seroprev <- c("16-24", "25-34", "35-49", "50-59", "60-64", "65-69",
                         "70-74", "75-79", "80+")
seroprev_Wales <- c(16.9, 12.5, 8.1, 5.5, 4.7, 3.7, 3.7, 3.1, 3.1)
seroprevL_Wales <- c(12.0, 9.0, 5.8, 3.9, 3.3, 2.6, 2.5, 2.0, 2.0)
seroprevH_Wales <- c(23.0, 17.3, 11.1, 7.7, 6.8, 5.4, 5.4, 4.7, 4.7)


### OUTCOMES
# DEATHS
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/weeklyprovisionalfiguresondeathsregisteredinenglandandwales
# Take total England + Wales stratified deaths, and scale by the proportion of
# deaths that occur in Wales
englandWalesDeaths <- read.csv("../data/downloaded_datasets/england/weekly_deaths_2020_ONS.csv",
                    stringsAsFactors=FALSE) %>%
  tidyr::pivot_longer(., cols=c(2:54), names_to="week", values_to="deaths") %>%
  dplyr::mutate(., week=as.integer(gsub("X", "", week))) %>%
  dplyr::filter(., week<=52) %>%
  dplyr::group_by(., Age) %>%
  dplyr::summarize(., deaths=sum(deaths))

# vector indicating to what seroprev bin each death group belongs
ageGroups <- c(0, 0, 0, 1, 1, 2, 2, 3, 3, 3, 0, 4, 4, 5, 6, 7, 8, 9, 9, 9)
englandWalesDeaths$bin <- ageGroups
englandWalesDeaths <- dplyr::group_by(englandWalesDeaths, bin) %>%
  dplyr::summarize(., deaths=sum(deaths))

# Total deaths segregated between England and Wales
totalEnglandWalesDeaths <- read.csv("../data/downloaded_datasets/england/wales_england_deaths_total.csv",
                                    stringsAsFactors=FALSE)
totalEnglandWalesDeaths$Date <- lubridate::dmy(totalEnglandWalesDeaths$Date)
totalEnglandWalesDeaths <- dplyr::filter(totalEnglandWalesDeaths, !is.na(Date)
                                         & Date<="2020-12-20")
propDeathsWales <- with(totalEnglandWalesDeaths, sum(Wales)/(sum(England)+sum(Wales)))

# Adjust age-stratified deaths
deaths_Wales <- round(englandWalesDeaths$deaths[2:10]*propDeathsWales)


### FULL TABLE
Wales_w1 <- data.frame(Age=age_Wales_seroprev,
                  Population=population_Wales2,
                  Prevalence=seroprev_Wales,
                  PrevalenceL=seroprevL_Wales,
                  PrevalenceH=seroprevH_Wales,
                  Prevalence_unvacc=seroprev_Wales,
                  PrevalenceL_unvacc=seroprevL_Wales,
                  PrevalenceH_unvacc=seroprevH_Wales,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Wales,
                  Deaths_unvacc=deaths_Wales,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Wales",
                  EndPointOutcome="2020-12-20",
                  MidPointCases="2020-12-20",
                  Period=1)



#############
### TIME 2
#############
### POPULATION
walesPopAge <- c("0-5", "5-9", "10-14", "15-19", "20-24", "25-29",
                   "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                   "60-64", "65-69", "70-74", "75-79", "80-84", "85-89",
                    "90+")
population_Wales <- c(161341, 182178, 184437, 173821, 206557, 208114,
                         196672, 185873, 172930, 192465, 216970, 222222,
                         197416, 179955, 181886, 131023, 90566, 54400, 43749)

population_Wales2 <- c(sum(population_Wales[4:5]),
                         sum(population_Wales[6:7]), sum(population_Wales[8:10]),
                         sum(population_Wales[11:12]), population_Wales[13:16], 
                         sum(population_Wales[17:19]))

### SEROPREVALENCE
# https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19antibodydatafortheuk/2021
age_Wales_seroprev <- c("16-24", "25-34", "35-49", "50-59", "60-64", "65-69",
                         "70-74", "75-79", "80+")
seroprev_Wales <- c(41.0, 45.9, 55.1, 84.9, 82.8, 84.6, 93.5, 94.0, 91.2)
seroprevL_Wales <- c(34.3, 39.4, 48.6, 81.1, 78.3, 80.4, 91.3, 91.7, 87.6) 
seroprevH_Wales <- c(48.1, 52.7, 61.3, 88.2, 86.5, 88.1, 95.2, 95.7, 93.8)


### OUTCOMES
# DEATHS
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/weeklyprovisionalfiguresondeathsregisteredinenglandandwales
englandWalesDeaths2020 <- read.csv("../data/downloaded_datasets/england/weekly_deaths_2020_ONS.csv",
                    stringsAsFactors=FALSE) %>%
  tidyr::pivot_longer(., cols=c(2:54), names_to="week", values_to="deaths") %>%
  dplyr::mutate(., week=as.integer(gsub("X", "", week))) %>%
  dplyr::group_by(., Age) %>%
  dplyr::summarize(., deaths=sum(deaths))
englandWalesDeaths2021 <- read.csv("../data/downloaded_datasets/england/weekly_deaths_2021_ONS.csv",
                    stringsAsFactors=FALSE) %>%
  tidyr::pivot_longer(., cols=c(2:54), names_to="week", values_to="deaths") %>%
  dplyr::mutate(., week=as.integer(gsub("X", "", week))) %>%
  dplyr::filter(., week<=18) %>%
  dplyr::group_by(., Age) %>%
  dplyr::summarize(., deaths=sum(deaths))
englandWalesDeaths2021$deaths <- englandWalesDeaths2020$deaths + englandWalesDeaths2021$deaths

# vector indicating to what seroprev bin each death group belongs
ageGroups <- c(0, 0, 0, 1, 1, 2, 2, 3, 3, 3, 0, 4, 4, 5, 6, 7, 8, 9, 9, 9)
englandWalesDeaths2021$bin <- ageGroups
englandWalesDeaths2021 <- dplyr::group_by(englandWalesDeaths2021, bin) %>%
  dplyr::summarize(., deaths=sum(deaths))

# Total deaths segregated between England and Wales
totalEnglandWalesDeaths <- read.csv("../data/downloaded_datasets/england/wales_england_deaths_total.csv",
                                    stringsAsFactors=FALSE)
totalEnglandWalesDeaths$Date <- lubridate::dmy(totalEnglandWalesDeaths$Date)
totalEnglandWalesDeaths <- dplyr::filter(totalEnglandWalesDeaths, !is.na(Date)
                                         & Date<="2021-05-09")
propDeathsWales <- with(totalEnglandWalesDeaths, sum(Wales)/(sum(England)+sum(Wales)))

# Adjust age-stratified deaths
deaths_Wales <- round(englandWalesDeaths2021$deaths[2:10]*propDeathsWales)

### VACCINES
# Source: 
# https://hduhb.nhs.wales/about-us/governance-arrangements/board-committees/audit-and-risk-assurance-committee-arac/arac/audit-and-risk-assurance-committee-meeting-22-june-2021/item-4-4-covid-19-vaccination-roll-out/
age_Wales_vaccine <- c("0-49", "50-59", "60-69", "70-79", "80+")
oneDoseWales <- c(NA, 91.3, 94.4, 96.6, 97.2)
# Match seroprevalence bins
oneDoseWales <- c(rep(oneDoseWales[1],3), oneDoseWales[2],
                  rep(oneDoseWales[3],2), rep(oneDoseWales[4],2), oneDoseWales[5])

### FULL TABLE
Wales_w2 <- data.frame(Age=age_Wales_seroprev,
                  Population=population_Wales2,
                  Prevalence=seroprev_Wales,
                  PrevalenceL=seroprevL_Wales,
                  PrevalenceH=seroprevH_Wales,
                  Prevalence_unvacc=NA,
                  PrevalenceL_unvacc=NA,
                  PrevalenceH_unvacc=NA,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Wales,
                  Deaths_unvacc=NA,
                  onePlusVaccine=oneDoseWales,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Wales",
                  EndPointOutcome="2021-05-09",
                  MidPointCases="2021-05-05",
                  Period=2)


#############
### Northern Ireland
#############

#############
### TIME 1
#############
### POPULATION
irelandPopAge <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
  "60-69", "70-79", "80-89", "90+")
population_Ireland <- c(245610, 238885, 232143, 250934, 240883, 258086,
                        199713, 146596, 68733, 13927)

population_Ireland2 <- c(population_Ireland[2]/3+population_Ireland[3]/2,
                         population_Ireland[3]/2+population_Ireland[4]/2,
                         population_Ireland[4]/2+population_Ireland[5],
                         sum(population_Ireland[6:7]),
                         sum(population_Ireland[8:10]))

### SEROPREVALENCE
# https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19antibodydatafortheuk/2021
age_Ireland_seroprev <- c("16-24", "25-34", "35-49", "50-69", "70+")
seroprev_Ireland <- c(7.7, 16.3, 7.6, 4.6, 2.5)
seroprevL_Ireland <- c(2.5, 8.2, 3.6, 2.3, 0.9)
seroprevH_Ireland <- c(20.2, 30.3, 15.1, 8.6, 6.6)

### OUTCOMES
# DEATHS
# Source: https://www.nisra.gov.uk/publications/excess-mortality-covid-19-related-deaths-december-2020
age_Ireland_deaths <- c("<1", "01-04", "05-09", "10-14", "15-19", "20-24",
                        "25-29", "30-34", "35-39", "40-44", "45-49", "50-54",
                        "55-59", "60-64", "65-69", "70-74", "75-79", "80-84",
                        "85-89", "90+")
deaths_Ireland <- c(0, 0, 0, 0, 1, 0, 0, 2, 3, 8, 5, 27, 44, 63, 94, 175, 268,
                   385, 406, 422)
# match serology bins
deaths_Ireland <- c(sum(deaths_Ireland[5:6]), sum(deaths_Ireland[7:8]),
                 sum(deaths_Ireland[9:11]), sum(deaths_Ireland[12:15]),
                 sum(deaths_Ireland[16:20]))

### FULL TABLE
Ireland_w1 <- data.frame(Age=age_Ireland_seroprev,
                  Population=population_Ireland2,
                  Prevalence=seroprev_Ireland,
                  PrevalenceL=seroprevL_Ireland,
                  PrevalenceH=seroprevH_Ireland,
                  Prevalence_unvacc=seroprev_Ireland,
                  PrevalenceL_unvacc=seroprevL_Ireland,
                  PrevalenceH_unvacc=seroprevH_Ireland,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Ireland,
                  Deaths_unvacc=deaths_Ireland,
                  onePlusVaccine=0,
                  twoPlusVaccine=0,
                  Type="Seroprevalence",
                  Location="Ireland",
                  EndPointOutcome="2020-12-30",
                  MidPointCases="2020-12-30",
                  Period=1)


#############
### TIME 2
#############
### POPULATION
irelandPopAge <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
  "60-69", "70-79", "80-89", "90+")
population_Ireland <- c(245610, 238885, 232143, 250934, 240883, 258086,
                        199713, 146596, 68733, 13927)

population_Ireland2 <- c(population_Ireland[2]/3+population_Ireland[3]/2,
                         population_Ireland[3]/2+population_Ireland[4]/2,
                         population_Ireland[4]/2+population_Ireland[5],
                         sum(population_Ireland[6:7]),
                         sum(population_Ireland[8:10]))

### SEROPREVALENCE
# https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19antibodydatafortheuk/2021
age_Ireland_seroprev <- c("16-24", "25-34", "35-49", "50-69", "70+")
seroprev_Ireland <- c(29.8, 44.8, 73.5, 89.7, 90.3)
seroprevL_Ireland <- c(20.9, 34.6, 64.5, 85.1, 85.7)
seroprevH_Ireland <- c(40.6, 55.9, 81.2, 93.0, 93.7)

### OUTCOMES
# DEATHS
# Source: https://www.nisra.gov.uk/publications/weekly-deaths-week-ending-07-may-2021
age_Ireland_deaths <- c("1-14", "15-44", "45-64", "65-74", "75-85", "85+")
deaths_Ireland <- c(0 , 25 , 254 , 429 , 1006 , 1246)

relIFR1 <- levinIFR(seq(15,44))/sum(levinIFR(seq(15,44)))
relIFR1 <- c(sum(relIFR1[1:10]), sum(relIFR1[11:20]), sum(relIFR1[21:30]))
relIFR2 <- levinIFR(seq(45,64))/sum(levinIFR(seq(45,64)))
relIFR2 <- c(sum(relIFR2[1:5]), sum(relIFR2[6:20]))
relIFR3 <- levinIFR(seq(65,74))/sum(levinIFR(seq(65,74)))
relIFR3 <- c(sum(relIFR3[1:5]), sum(relIFR3[6:10]))

deaths_Ireland <- c(round(deaths_Ireland[2]*relIFR1[1]),
                   round(deaths_Ireland[2]*relIFR1[2]),
                   round(deaths_Ireland[2]*relIFR1[3]+deaths_Ireland[3]*relIFR2[1]),
                   round(deaths_Ireland[3]*relIFR2[2]+deaths_Ireland[4]*relIFR3[1]),
                   round(deaths_Ireland[4]*relIFR3[2]+sum(deaths_Ireland[5:6])))

### VACCINATION
# Vaccine deaths: https://www.health-ni.gov.uk/publications/vaccination-status-deaths-and-hospitalisations-previous-publications
age_Ireland_vacc <- c("16-24", "25-34", "35-49", "50-69", "70+")
oneVaccIreland <- c(17.8, 41.1, 85.1, 95.5, 98.9)
twoVaccIreland <- c(9.5, 18.0, 22.1, 49.4, 87.7)


### FULL TABLE
Ireland_w2 <- data.frame(Age=age_Ireland_seroprev,
                  Population=population_Ireland2,
                  Prevalence=seroprev_Ireland,
                  PrevalenceL=seroprevL_Ireland,
                  PrevalenceH=seroprevH_Ireland,
                  Prevalence_unvacc=NA,
                  PrevalenceL_unvacc=NA,
                  PrevalenceH_unvacc=NA,
                  Severe=NA,
                  Critical=NA,
                  Deaths=deaths_Ireland,
                  Deaths_unvacc=NA,
                  onePlusVaccine=oneVaccIreland,
                  twoPlusVaccine=twoVaccIreland,
                  Type="Seroprevalence",
                  Location="Ireland",
                  EndPointOutcome="2021-05-09",
                  MidPointCases="2021-05-05",
                  Period=2)



countriesDf <- rbind(Spain_w1, Spain_w2, Atlanta_w1, Atlanta_w2,
                     France_w1, France_w2, Geneva_w1, Geneva_w2, Geneva_w3,
                     Neuchatel_w1, Neuchatel_w2, Neuchatel_w3,
                     Fribourg_w1, Fribourg_w2, Fribourg_w3,
                     Ticino_w1, Ticino_w2, Ticino_w3, Ticino_w4,
                     Freiburg_w1, Freiburg_w2, Reutlingen_w1, Reutlingen_w2,
                     Osnabruck_w1, Osnabruck_w2, Aachen_w1, Aachen_w2,
                     Magdeburg_w1, Magdeburg_w2,
                     Munich_w1, Munich_w2, Denmark_w1, Denmark_w2,
                     England_deaths_w1, England_deaths_w2,
                     England_deaths_w3, England_deaths_w4,
                     England_hospital_w1, England_hospital_w2,
                     England_hospital_w3, England_hospital_w4,
                     Scotland_w1, Scotland_w2, Wales_w1, Wales_w2,
                     Ireland_w1, Ireland_w2)


write.csv(countriesDf, "../data/raw_data/collected_data.csv", row.names=FALSE)


############################
## California
############################
#
## Hospital https://beta.healthdata.gov/Hospital/COVID-19-Reported-Patient-Impact-and-Hospital-Capa/g62h-syeh
## https://github.com/datadesk/california-coronavirus-data/blob/master/cdph-age.csv
#
## Unstratified ICU https://public.tableau.com/views/COVID-19HospitalsDashboard/Hospitals?:embed=y&:showVizHome=no
## https://data.ca.gov/dataset/covid-19-hospital-data1
##   opendata@cdph.ca.gov
## https://jamanetwork.com/journals/jama/article-abstract/2765303
## https://www.bmj.com/content/369/bmj.m1923.full
## https://www.cambridge.org/core/journals/journal-of-clinical-and-translational-science/article/clinical-characteristics-associated-with-covid19-severity-in-california/B58EB9C431C6404D867BF70DBCAEBA19
## https://www.nytimes.com/interactive/2021/us/california-covid-cases.html
## https://gis.cdc.gov/grasp/covidnet/COVID19_5.html
## https://github.com/datadesk/california-coronavirus-data
#
#
## Seroprevalence & pop: SARS-CoV-2 Cumulative Incidence and Period Seroprevalence: Results From a Statewide PopulationBased Serosurvey in California
#age_California_seroprev <- c("18-34", "35-44", "45-54", "55-64", "65+")
#seroprev_California <- c(5.2, 6.2, 8.7, 1.9, 0.8)
#seroprevL_California <- c(2.8, 3, 2.5, 0.7, 0.1)
#seroprevH_California <- c(9.5, 12.3, 25.6, 5, 4.5)
#
#population_California <- c(9730987, 5282100, 4979745, 4786635, 5838115)
#
## Deaths https://data.chhs.ca.gov/dataset/covid-19-time-series-metrics-by-county-and-state
#deathDataCA <- read.csv("../data/downloaded_datasets/california/covid19casesdemographics.csv",
#                       stringsAsFactors=FALSE) %>%
#  as_tibble(.) %>%
#  dplyr::filter(., (lubridate::date(report_date)=="2020-11-02") &
#                demographic_category=="Age Group" &
#                !(demographic_value %in% c("Missing", "Total")))
#
## Redistribute deaths to match seroprevalence bins
#relIFR1 <- levinIFR(seq(18,49))/sum(levinIFR(seq(18,49)))
#relIFR1 <- c(sum(relIFR1[1:17]), sum(relIFR1[18:27]),
#                sum(relIFR1[28:32]))
#relIFR2 <- levinIFR(seq(50:64))/sum(levinIFR(seq(50:64)))
#relIFR2 <- c(sum(relIFR2[1:4]), sum(relIFR2[5:15]))
#
#caDeaths <- deathDataCA$deaths
#deaths_California <- c(round(caDeaths[2]*relIFR1[1:2]),
#                       round(caDeaths[2]*relIFR1[3]+caDeaths[3]*relIFR2[1]),
#                       round(caDeaths[3]*relIFR2[2]), caDeaths[4])
#
#California <- data.frame(Age=age_California_seroprev,
#                  Population=population_California,
#                  Prevalence=seroprev_California,
#                  PrevalenceL=seroprevL_California,
#                  PrevalenceH=seroprevH_California,
#                  Severe=NA,
#                  Critical=NA,
#                  Deaths=deaths_California,
#                  Type="Seroprevalence",
#                  Location="California",
#                  EndPointOutcome="2020-11-02",
#                  MidPointCases="2020-11-02",
#                  Period=2)
#

