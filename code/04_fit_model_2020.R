#############################################
#############################################
# 
# This script generates the fits for the main results of the
# paper. It fits the Bayesian model described in the manuscript
# and defined in 6_estimate_serology_outcome_rates.stan to the
# data collected in script 2_countries_values.R, and corrected
# in 4_countries_values_correction.R
# 
# This script also generates all the models reported in the
# supplementary, by taking subsets of the main dataset.
# 
# 
#############################################
#############################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(rstan)
library(bayesplot)
library(tidybayes)
library(rriskDistributions)
source("./functions_auxiliary.R")
source("./stan_utility.R")
set.seed(2691)

# Fitting parameters
nChains <- 4
nCores <- 2
nIter <- 4000

# compile model
outcome_reg <- rstan::stan_model("./severity_change.stan",
                                 model_name="time_change_IFR",
                                 warn_pedantic=TRUE)

# Load the data
simData <- read.csv("../data/raw_data/simulated_data.csv",
                        stringsAsFactors=FALSE) %>%
  as_tibble(.) %>%
  dplyr::mutate(., Prevalence1=Prevalence1/100, Prevalence1L=Prevalence1L/100,
                Prevalence1H=Prevalence1H/100,
                Prevalence2=Prevalence2/100, Prevalence2L=Prevalence2L/100,
                Prevalence2H=Prevalence2H/100) %>%
  dplyr::filter(., Prevalence1>0)
simData$meanAge <- mid_bin_age(simData$Age)

####################
# define function to make list of initial values, for STAN
####################
initial_values <- function(nChains, paramListName, lowerUni, upperUni,
                           paramSize, prevalenceVals1, prevalenceVals2) {
  initList <- list()
  for (ch in c(1:nChains)) {
    initList[[ch]] <- list()
    for (p in c(1:length(paramListName))) {
      initList[[ch]][[paramListName[p]]] <- runif(paramSize[p], min=lowerUni[p],
                                                  max=upperUni[p])
    }
    #randomVec <- runif(length(prevalenceVals1), min=1, max=1.05)
    #variatedPrev <-  prevalenceVals1*randomVec
    initList[[ch]][["prevalence_raw1"]] <- prevalenceVals1
    prevalenceDiffs <- pmax(prevalenceVals2-prevalenceVals1, 0.05)
    initList[[ch]][["prevalence_diff"]] <- prevalenceDiffs
  }
  return(initList)
}

# Set the ranges for the initial values of the parameters
paramListName <- c("ageSlope", "intercept", "interceptChage", "interceptSigma",
"changeSigma", "locationIntercept", "locationChange")
lowerUni <- c(2, -14, -1, 0.5, 0.1, -14, -1)
upperUni <- c(4, -12, 0, 2, 0.2, -12, 0)
paramSize <- c(1, 1, 1, 1, 1, length(unique(simData$Location)),
               length(unique(simData$Location)))

##############
# Fit beta distribution to seroprevalences
##############
fittedBetas1 <- fit_beta_ci(meanEstimate=simData$Prevalence1,
                            lower=simData$Prevalence1L,
                            upper=simData$Prevalence1H)
simData$seroprev1Shape1 <- fittedBetas1$shape1
simData$seroprev1Shape2 <- fittedBetas1$shape2
fittedBetas2 <- fit_beta_ci(meanEstimate=simData$Prevalence2,
                            lower=simData$Prevalence2L,
                            upper=simData$Prevalence2H)
simData$seroprev2Shape1 <- fittedBetas2$shape1
simData$seroprev2Shape2 <- fittedBetas2$shape2

# get values for standardization
meanAge <- 0
sdAge <- sd(simData$meanAge)
simData$stdAge <- (simData$meanAge-meanAge)/sdAge


# Make a list with the input we pass to STAN
simDataList <- list(N=nrow(simData),
                  K=length(unique(simData$Location)),
                  location=simData$Location,
                  ageVec=simData$stdAge,
                  population=simData$Population,
                  seroprev1Shape1=simData$seroprev1Shape1,
                  seroprev1Shape2=simData$seroprev1Shape2,
                  seroprev2Shape1=simData$seroprev2Shape1,
                  seroprev2Shape2=simData$seroprev2Shape2,
                  outcomes1=simData$Deaths1,
                  outcomes2=simData$Deaths2-simData$Deaths1)

########
# Fit the model
########
# Sample initial values
initList <- initial_values(nChains=nChains, paramListName=paramListName,
                           lowerUni=lowerUni, upperUni=upperUni,
                           paramSize=paramSize, prevalenceVals1=simData$Prevalence1,
                           prevalenceVals2=simData$Prevalence2)

# Fit model
model <- rstan::sampling(outcome_reg, data=simDataList,
                           chains=nChains, iter=nIter, refresh=0,
                           verbose=TRUE, cores=nCores, init=initList)

######
# Extract model data, and check diagnostics
######
posteriorTraces <- tidybayes::gather_draws(model, ageSlope, intercept,
                                           interceptChange, interceptSigma,
                                           changeSigma, locationIntercept[loc],
                                           locationChange[loc])

write.csv(posteriorTraces, "../data/analysis_results/simulated_data_fit.csv",
          row.names=FALSE)

#serologyTrace <- dplyr::mutate(posteriorTraces, .chain=factor(.chain))  %>%
#  dplyr::filter(., !(.variable %in% c("locationIntercept", "locationChange"))) %>%
#  ggplot(., aes(x=.iteration, y=.value, color=.chain)) +
#  geom_line() +
#  facet_wrap(.~.variable, scales="free") +
#  theme_bw()
#
#pairsPlot <- pairs(model, pars=c("ageSlope", "intercept", "interceptChange",
#                                 "interceptSigma", "changeSigma"))
##png("../data/figures/simulated_data_fit_pairs_plot.png")
#
#pairsPlot <- pairs(model, pars=c("prevalence_diff"))

